
# Script for analysis of single kit TMT quantification from Proteome Discover 2.5 output

library(tidyverse)
library(janitor)
library(data.table)
setDTthreads(4)
getDTthreads()

#library(ggrepel)
#library(svglite)
#library(scales)

# load in the peptide and peptide file names as .txt files
peptide_data = c("Muoio_BeckyWilson_STIM1-SKM_TMT10_FINAL_2020-05-03_SN2pt5_PeptideIsoforms")
protein_data = c("Muoio_BeckyWilson_STIM1-SKM_TMT10_FINAL_2020-05-03_SN2pt5_Proteins")

# load in mito_carta file for mitochondria data
mito_carta = fread("Mouse.MitoCarta3.0.txt", header = T)
colnames(mito_carta)[[1]] = "GeneName"
  
# read in the peptide and protein data
peptides = fread(paste0(peptide_data, ".txt"), header = T)
proteins = fread(paste0(protein_data, ".txt"), header = T)

# variables for the samples being compared and the total number of samples
ID_names = c("WT", "KO", "pool") 
total_number_of_samples = 11



# list of all description columns the user wants to pull from the raw data
# in Proteome Discover V2.5 we export the protein and peptide isoform data with all available
# columns. We then only choose a subset of them to keep in the post-analysis tables
# t

desc_cols = function(method){
  
  if (method == 1){
    
    list = list(description_cols = make_clean_names(c("Description", "Accession", "Master", "Exp. q-value",
                                     "# Peptides", "# PSMs", "# Protein Unique Peptides", "# Unique Peptides", "Entrez Gene ID",
                                     "Coverage [%]", "Abundance Ratio: (liver (tissue), input (fraction), TRF (feeding)) / (liver (tissue), input (fraction), control (feeding))")),
                
                description_cols_peptide = make_clean_names(c("Master Protein Accessions", "Protein Accessions","Sequence",
                                                              "# Missed Cleavages", "PSM Ambiguity","# PSMs", 
                                             "Modifications", "Modifications in Proteins",
                                             "Abundance Ratio: (liver (tissue), phospho (fraction), TRF (feeding)) / (liver (tissue), phospho (fraction), control (feeding))")))
    }
  
  else if(method == 2){
    
    list = list(description_cols = make_clean_names(c("Description", "Accession", "Master", "Exp. q-value",
                                     "# Peptides", "# PSMs", "# Protein Unique Peptides", "# Unique Peptides", "Entrez Gene ID",
                                      "Reactome Pathways", "WikiPathways")),
                
                description_cols_peptide = make_clean_names(c("Master Protein Accessions", "Protein Accessions", "Sequence",
                                             "# Missed Cleavages", "PSM Ambiguity", "# PSMs", 
                                             "Modifications", "Modifications in Proteins", "XCorr (by Search Engine): Sequest HT",
                                             "Deltam/z [Da] (by Search Engine): Sequest HT"
                                             )))
  }
}

description_cols = desc_cols(2)


# reducing cols function
# reduces the number of columns. Keeps the columns set in the description cols list, and abundance quantification
# values
# protein1 removes and renames the protein dataframe based on colnames and description_cols
# protein2 retains is attempt to normalize between TMT kits, it will pull the all input abundance data

# df1 
reduce_cols = function(df1, df2, method, fraction){
  if(method == "protein"){
    
    df1 = clean_names(df1)
    
    reduced_protein <- df1 %>%
      setnames(names(select(., contains("q_value"))), "exp_q_value") %>%
      select(., all_of(description_cols[[1]]), starts_with(paste0("abundance_", paste0(fraction, "_")))) 

    # splits the description column and extracts just the gene name
    split_ref = as.data.frame(str_split_fixed(reduced_protein$description, "GN=",2))
    split_ref2 = as.data.frame(str_split_fixed(split_ref$V2, "PE=", 2))
    rm_white = as.data.frame(str_replace_all(split_ref2$V1, pattern = " ", replacement = ""))
    colnames(rm_white) = c("gene_name")
    
    reduced_protein2 = cbind(rm_white, reduced_protein)
    
    reduced_protein3 = reduced_protein2 %>%
      mutate(mitocarta = gene_name %in% mito_carta$symbol)
    
    x = select(df1, starts_with(paste0("abundance_", paste0(fraction, "_"))))
    
    print("column names that have been pulled, please confirm if correct and in proper order")
    print(names(x))
    
    return(reduced_protein3)
    
    
  }
  
  else if(method == "peptide-phos"){

    df1 = clean_names(df1)
    df2 = clean_names(df2)
    
    reduced_peptides <- df1 %>%
      select(., all_of(description_cols[[2]]), starts_with(paste0("abundance_", paste0(fraction, "_"))))
      #rename_at(vars(starts_with(paste0("abundance_", paste0(fraction, "_")))), ~ paste0("abundance_phos_", column_names, sep = " ")) 

    # split the master protein accessions column into two, in order to get the first accession number
    accession = as.data.frame(str_split_fixed(df1$`master_protein_accessions`, ";", 2))
    colnames(accession)[1] = c("accession")
    
    # add it back to the df
    reduced_peptides2 = cbind(reduced_peptides, accession[1])
    
    # pull the accession and description columns from the protein df and join them with the
    # reduced_peptides2 data
    proteins2 = select(df2, accession, description, "entrez_gene_id")
    
    joining = reduced_peptides2 %>%
      left_join(., proteins2, by = "accession") %>%
      mutate(mitocarta = entrez_gene_id %in% mito_carta$GeneName)
    
    # splits the description column and extracts just the gene name
    split_ref = as.data.frame(str_split_fixed(joining$description, "GN=",2))
    split_ref2 = as.data.frame(str_split_fixed(split_ref$V2, "PE=", 2))
    rm_white = as.data.frame(str_replace_all(split_ref2$V1, pattern = " ", replacement = ""))
    colnames(rm_white) = c("gene_name")

    reduced_peptides3 = cbind(rm_white, joining)
    
    # compares the gene_name column to the mitocarta symbol column to identify mito proteins
    # also reorders the columns, this will likely break if the column names get changed in a PD update
    reduced_peptides4 = reduced_peptides3 %>%
      mutate(mitocarta = gene_name %in% mito_carta$symbol) %>%
      rename("first_master_protein_accession" = accession) %>%
      select(., gene_name, first_master_protein_accession, master_protein_accessions, protein_accessions, description, mitocarta, sequence,
             number_missed_cleavages, psm_ambiguity, number_ps_ms, modifications, modifications_in_proteins, entrez_gene_id, contains("abundance"))
    
    return(reduced_peptides4)
    
  }
  
  else if(method == "peptide-acetyl"){
    
    #starts_with(paste("Abundance:", paste(fraction, ":", sep = ''), sep = ' '))
    # reduces the columns and renames the acetyl peptide abundance column names
    reduced_peptides <- df1 %>%
      select(., all_of(description_cols[[2]]), contains(paste("Abundance:", paste(fraction, ":", sep = ''), sep = ' '))) %>%
      rename_at(vars(starts_with(paste("Abundance:", paste(fraction, ":", sep = ''), sep = ' '))), ~ paste("abundance_acetyl_", column_names, sep = " ")) %>%
      clean_names() 
    
    
    
    # split the master protein accessions column into two, in order to get the first accession number
    Accession = as.data.frame(str_split_fixed(df1$`Master Protein Accessions`, ";", 2))
    colnames(Accession)[1] = c("Accession")
    
    # add it back to the df
    reduced_peptides2 = cbind(reduced_peptides, Accession[1])
    
    # pull the accession and description columns from the protein df and join them with the
    # reduced_peptides2 data
    proteins2 = select(df2, Accession, Description)
    
    joining = reduced_peptides2 %>%
      left_join(., proteins2, by = "Accession") 
    
    # splits the description column and extracts just the gene name
    split_ref = as.data.frame(str_split_fixed(joining$Description, "GN=",2))
    split_ref2 = as.data.frame(str_split_fixed(split_ref$V2, "PE=", 2))
    colnames(split_ref2) = c("gene_name", "toss")
    
    reduced2 = cbind(joining, split_ref2[1])
    
    # creates the final version of the data frame and renames the new accession column
    final = reduced2 %>%
      rename("1st_master_protein_accession" = Accession)
    
    
    return(final)
    
  }
  
}

reduced_cols_protein = reduce_cols(proteins, proteins, method = "protein", "F2")

reduced_cols_peptide = reduce_cols(peptides, proteins, method = "peptide-phos", "F6")

# reducing rows of dataframe
# "tmt_protein" is for only keeping IsMasterProtein, and qvalue < 0.01 and removes proteins that are
# only in 1/2 the samples
# "tmt_peptide_phos" filters data that has "Phospho" in modification column

reduce_rows = function(df, method){
  if(method == "tmt_protein"){
    
    reduced_df = df %>%
      mutate(na_row_count = rowSums(is.na(select(df, starts_with("abundance"))))) %>%
      select(everything(), -starts_with("abundance"),  na_row_count, starts_with("abundance")) %>%
      filter(na_row_count <= total_number_of_samples/2) %>%
      filter(master == "IsMasterProtein" & exp_q_value < 0.01)
      
    change2 = function(df, oldname, newname){
      df %>%
        rename_at(vars(starts_with(oldname)), list(~ paste0(newname)))
    }
    
    reduced_df2 = change2(reduced_df, oldname = "abundance_ratio", newname = "pd ratio")
    
    return(reduced_df2)
    
  }
  
  else if(method == "tmt_peptide_phos"){
    # removes peptides that don't have quantifications for at least half the peptides and
    # removes peptides that don't have phospho modification listed in the modifications column

    description_cols_peptide = make_clean_names(description_cols[[2]])
    
    reduced_df <- df %>%
      mutate(na_row_count = rowSums(is.na(select(df, starts_with("abundance"))))) %>%
      select(everything(), -starts_with("abundance"), mitocarta, na_row_count, starts_with("abundance")) %>%
      filter(na_row_count <= total_number_of_samples/2)  %>%
      filter(., grepl('Phospho', modifications)) 
    
    change2 = function(df, oldname, newname){
      df %>%
        rename_at(vars(starts_with(oldname)), list(~ paste0(newname)))
    }
    
    reduced_df2 = change2(reduced_df, oldname = "abundance_ratio", newname = "pd ratio")
    
    return(reduced_df2)

  }
  
}

reduced_rows_protein = reduce_rows(reduced_cols_protein, "tmt_protein")

reduced_rows_peptide = reduce_rows(reduced_cols_peptide, "tmt_peptide_phos")


# normalization
# Normalizes the data to the ratio of the  colsums / mean(columns)
# the df1 argument is ideally the reduced_rows of either the protein or peptide data
# the df2 argument is the 

# change the fraction to the fraction in which the input is in

normalize_df = function(df1, df2, method, fraction, show_norm_boxplot) {
  
  if(method == "protein"){
    
    # normalize the protein abundance between sample by using the ratio (sum/avg) of the peptide data
    # make sure to use the original peptide data and not a df that was already reduced
    
    norm_names_df1 = grep("abundance", names(df1), value = TRUE)

    reduced_pep = df2 %>%
      select(starts_with(paste("Abundance:", paste(fraction, ":", sep = ''), sep = ' '))) %>%
      mutate_if(is.character, as.numeric)

    sum = colSums(reduced_pep, na.rm = T)
    avg = mean(sum)
    
    ratio = sum/avg
    
    show(paste("Ratio of peptide"))
    show(round(ratio, digits = 4))
    
    bar = barplot(ratio, names.arg = norm_names_df1, ylab = "ratio", ylim = c(0,ceiling(max(ratio,na.rm=T))), 
                  main = paste("Ratio of peptide input"))
    bar
    
    # normalize protein data based on peptide ratios
    
    
    adj_by_ratio = sweep(select(df1, all_of(norm_names_df1)), 2, ratio, "/")

    #norm_to_pep = cbind(df1, adj_by_ratio)
    
    
    #log2 trasnform the data
    
    logged_norm = apply(select(adj_by_ratio, all_of(norm_names_df1)), 2, log2)
    
    log_avg = rowMeans(logged_norm, na.rm = TRUE)
    
    norm_to_avg = logged_norm - log_avg
    
    final = as.data.frame(cbind(select(df1, -starts_with("abundance")), norm_to_avg))
    
    if (show_norm_boxplot == T){
      
      before = df1 %>%
        select(starts_with("abundance"))
      
      box1 = boxplot(log2(before), col = "blue", main = "Before Normalization")
      box1
      

      box2 = boxplot(logged_norm, col = "green", main = "After Loading Normalization")
      box2
      

      box3 = boxplot(norm_to_avg, col = "red", main = "After Normalizing To The Mean")
      box3
      
    }
    
    return(final)
    
  }
  
  else if(method == "peptide"){
    
    norm_names_df1 = grep("abundance", names(df1), value = TRUE)
    norm_names_df2 = grep("abundance", names(df2), value = TRUE)
    
    # normalize all the peptide abundance values based on the ratio of the sum of each column
    # divided by the avg of the sums.
    # be sure to use the unprocessed peptide data
    
    reduced_prot = peptides %>%
      select(starts_with(paste("Abundance:", paste(fraction, ":", sep = ''), sep = ' '))) %>%
      mutate_if(is.character, as.numeric)
    
    
    sum = colSums(reduced_prot, na.rm = T)
    avg = mean(sum)
    
    ratio = sum/avg
    
    show(paste("Ratio of peptide"))
    show(round(ratio, digits = 4))
    
    bar = barplot(ratio, names.arg = norm_names_df1, ylab = "ratio", ylim = c(0,ceiling(max(ratio,na.rm=T))), 
                  main = paste("Ratio of peptide samples"))
    bar
    
    # normalize protein data based on peptide ratios
    
    adj_by_ratio = sweep(select(df1, all_of(norm_names_df1)), 2, ratio, "/")
    
    
    #description_cols_peptide = make_clean_names(description_cols[[2]])
    
    norm_to_pep = cbind(select(df1, -starts_with("abundance")), adj_by_ratio)
    
    
    # log transform data, take the mean of all samples and subtract each sample from the mean
    
    logged_norm = apply(select(norm_to_pep, all_of(norm_names_df1)), 2, log2)
    logged_norm_to_pep = cbind(select(norm_to_pep, -contains("abundance")), logged_norm)
    
    log_avg = rowMeans(select(logged_norm_to_pep, all_of(norm_names_df1)), na.rm = T)
    
    norm_to_avg = logged_norm - log_avg
    
    logged_norm_to_pep2 = cbind(select(norm_to_pep, -contains("abundance")), norm_to_avg)
    
    # pulling the normalized protein input data from reduced_rows_protein into the peptide data
    
    norm_to_prot <- rename(logged_norm_to_pep2, "accession" = "first_master_protein_accession")
    
    protein_norm_df = normalized_protein %>%
      select(accession, contains(norm_names_df2))
    
    join = left_join(norm_to_prot, protein_norm_df, by = "accession")

    relative_ptm_occupancy = select(join, contains(norm_names_df1)) - select(join, contains(norm_names_df2))
    colnames(relative_ptm_occupancy) = paste0("relative_occupancy_", norm_names_df1)
    
    
    final = as.data.frame(cbind(select(df1, -contains("abundance")), norm_to_avg, relative_ptm_occupancy))
    
    return(final)
  }
  
  
}

normalized_protein = normalize_df(reduced_rows_protein, peptides, "protein", "F2", T)

normalized_peptide = normalize_df(reduced_rows_peptide, reduced_rows_protein, "peptide", "F2")


# quantification
samples_to_compare = c("KO", "WT")
number_of_reps = c(5,5)
control = c("WT")
variables = c("KO")

# type = 1 is the standard proteomics analysis, which consists of getting the mean, SD, FC, two.sided t.test,
# and the adjusted p.value (q.value) using BH method
# type = 2 uses the limma package to calculate the FC (coefficient) and the adjusted p.value
# each type is split into 2 parts, the first is if the method is for just abundance, the second is for if the method is more than 2

new_final_df = function(df, method, type){
  if (type == 1){
    
    multi_method = lapply(method, function(z){
      
      df_names = df %>%
        select(., starts_with(z)) %>%
        names()
      
      sample_names = lapply(make_clean_names(ID_names), function(x) grep(x, df_names, value = T, fixed=T))
      
      protein2 = lapply(sample_names, function(x) as.matrix(subset(df, select = x)))
      
      means1 = lapply(protein2, function(x){
        df2 = as.matrix(x)
        means = rowMeans(df2, na.rm=T)
      })
      names(means1) = c(paste0(z, "_mean_", ID_names))
      
      
      sd = lapply(protein2, function(x){
        df2 = as.matrix(x)
        sd = apply(df2, 1, sd, na.rm=T)
      })
      names(sd) = c(paste0(z, "_sd_", ID_names))
      
      control_1 = means1[paste0("mean_", control)]
      variables_1 = means1[paste0("mean_", variables)]
      
      # try this do.call(cbind, purrr::map2(BASE, FUTURE, ~ .x[, 4] - .y[, 4]))
      fc = vector('list', length(variables))
      names(fc) = c(paste(z, "_FC_", variables, "-", control, sep = ""))
      for(i in 1:length(variables)){
        fc[[i]] = unlist(variables_1[[i]] - unlist(control_1[[1]]))
      }
      
      df2 = df %>%
        dplyr::select(.,contains(z))
      
      control_sample_names = lapply(make_clean_names(control), function(x) grep(x, df_names, value = T, fixed=T))
      variables_sample_names = lapply(make_clean_names(variables), function(x) grep(x, df_names, value = T, fixed=T))
      t_test2 = vector('list', length(variables))
      names(t_test2) = c(paste("p.value_", variables, "_to_", control, sep = ""))
      for(i in 1:length(variables)){
        t_test2[[i]] = apply(df2, 1, function(y) {tryCatch(t.test(y[control_sample_names[[1]]], y[variables_sample_names[[i]]], alternative = "two.sided", var.equal = T)$p.value, error = function(err){return(NA)})})
      }
      
      q_value = as.data.frame(p.adjust(unlist(t_test2), method = "BH"))
      stat_results = as.data.frame(cbind(p.value = t_test2, q.value =  q_value))
      
      colnames(stat_results) = paste0(z, c("_p.value_", "_q.value_"), variables, "_over_", control)
      
      final = as.data.frame(cbind(means1, sd, stat_results))

    })
    
    final1 = as.data.frame(cbind(df, multi_method))
  }
    else if (type == 2){
    library(limma)
    
    multi_method = lapply(method, function(z){
      
      df_names = df %>%
        dplyr::select(., starts_with(z)) %>%
        names()
      
      sample_names = lapply(make_clean_names(ID_names), function(x) grep(x, df_names, value = T, fixed=T))
      
      protein2 = lapply(sample_names, function(x) as.matrix(subset(df, select = x)))
      
      means1 = lapply(protein2, function(x){
        df2 = as.matrix(x)
        means = rowMeans(df2, na.rm=T)
      })
      names(means1) = c(paste0(z, "_mean_", ID_names))
      
      
      sd = lapply(protein2, function(x){
        df2 = as.matrix(x)
        sd = apply(df2, 1, sd, na.rm=T)
      })
      names(sd) = c(paste0(z, "_sd_", ID_names))
      
      # this is where i start the limma analysis
      # use the normalzied data, it is already log2 transformed
      
      df1 = df %>%
        dplyr::select(., starts_with(z)) %>%
        dplyr::select(., contains(control) | contains(variables))
      
      # Design linear model with no intercept and no interaction
      # Group is the Genotype
      group_fct = factor(paste(rep(samples_to_compare, number_of_reps)))
      
      group_design = model.matrix(~0 + group_fct)
      
      # Rename column names in the design
      colnames(group_design) = samples_to_compare
      
      # Limma Step 1: Least Squares Estimates 
      group_fit = lmFit(df1, group_design)
      
      # Generate the contrast maps
      group_contMat = makeContrasts(VariableVsControl = paste0(variables, "-", control),
                                    levels=group_design)
      
      # Fit the Contrasts to the product of lmFit
      group_fit_cont = contrasts.fit(group_fit, group_contMat)
      
      # Now perform the eBayes on this new fit
      group_fit_cont_eb = eBayes(group_fit_cont)
      
      group_fit_cont_eb_decide = decideTests(group_fit_cont_eb, 
                                              method = "separate", 
                                              adjust.method = "BH", 
                                              p.value = 0.05)
      
      q_value = as.data.frame(p.adjust(group_fit_cont_eb$p.value, method = "BH"))
      
      limma_analysis = as.data.frame(cbind(group_fit_cont_eb$coefficients, group_fit_cont_eb$p.value, 
                                           q_value, group_fit_cont_eb_decide))
      
      colnames(limma_analysis) = paste0(z, c("_FC_", "_p.value_", "_q.value_", "_significance_"), variables, "over", control)
      
      final = as.data.frame(cbind(means1, sd, limma_analysis))
      
      
    })
    
    final1 = as.data.frame(cbind(df, multi_method))
      
    }
    
}

final_proteins = new_final_df(df = normalized_protein, method = list("abundance"), type = 1)

final_peptides = new_final_df(df = normalized_peptide, method = list("abundance", "relative_occupancy"), type = 1)



# exporting the data to a csv.
to_export = list(final_proteins, final_peptides_abundance, final_peptides_occupancy)
file_names = paste0(c(protein_data, paste0(peptide_data, "_abundance"), paste0(peptide_data, "_occupancy")), "_REDUCED")

output_csv = function(to_export, file_names){
  fwrite(to_export, file = paste0(file_names, ".csv"))
}

list(to_export = to_export,
     file_names = file_names) %>%
  pmap(output_csv)

###
###



# write file to excel spreadsheet
library(writexl)
write_df = function(df2, df3, df4, df5, df6, df7, df8, df9){
  library(writexl)
  data_reduced_name = paste(protein_data,"_REDUCED.xlsx", sep = "")
  write_xlsx(list(reduced_cols_peptide=df2, reduced_cols_protein = df3, 
                  reduced_rows_peptide=df4, reduced_rows_protein = df5, normalized_peptide = df6, normalized_protein = df7,
                  final_peptide = df8, final_protein = df9), data_reduced_name, format_headers = T)
  
}

write = write_df(reduced_cols_peptide, reduced_cols_protein, reduced_rows_peptide, reduced_rows_protein,
                 normalized_peptide, normalized_protein, final_peptide, final_protein)
###
###


####
# LIMMA ANALYSIS

library(limma)

samples_to_compare = c("control", "TRF")
number_of_each_samples = c(5,5)
control = c("control")
variables = c("TRF")

limma_analyis = function(df){
  
  # use the normalzied data, it is already log2 transformed
  
  df1 = df %>%
    select(starts_with("abundance"))
  
  # Design linear model with no intercept and no interaction
  # Group is the Genotype
  group_fct = factor(paste(rep(samples_to_compare, number_of_each_samples)))
  
  group_design = model.matrix(~0+group_fct)

  # Rename column names in the design
  colnames(group_design) = samples_to_compare

  # Limma Step 1: Least Squares Estimates 
  group_fit = lmFit(df1, group_design)

  # Generate the contrast maps
  group_contMat = makeContrasts(TRFvscontrol = TRF - control,
                              levels=group_design)
  

  
  # Fit the Contrasts to the product of lmFit
  group_fit_cont = contrasts.fit(group_fit, group_contMat)

  # Now perform the eBayes on this new fit
  group_fit_cont_eb = eBayes(group_fit_cont)

  
  # Perform the comparisons, with each contrast individually
  # Use Benj. Hoch. for With FDR = 0.05
  group_fit_cont_eb_decide = decideTests(group_fit_cont_eb, 
                                       method="separate", 
                                       adjust.method = "BH", 
                                       p.value = 0.05)
  
  # Export the data
  #fwrite(group_fit_cont_eb_decide, file = "test_limma.csv")
  #write.fit(group_fit_cont_eb, group_fit_cont_eb_decide, "test_limma_skm_trf.tsv")

  
}

limma = limma_analyis(normalized_protein)


####




# volcano  plot
# need to load htmlwidgets in order to save/export
# types are all_protein, all_peptide, mito_protein, mito_peptide
library(htmlwidgets)
library(plotly)
volcano_plot = function(df, type, title, exp, print, plotly_plot){
  options(ggrepel.max.overlaps = 100) 
  if(type == "all_peptide"){
    
    volcano <- df %>% 
      select(gene_name, modifications_in_proteins, mitocarta, FC = contains("FC"), t.test = starts_with("p.value"), q_value = starts_with("q_value")) %>%
      mutate(t_test = -log10(t.test), .keep = c("unused")) %>%
      filter(!is.na(FC))
    
    
    
    
    if(plotly_plot == TRUE){
      
      q = ggplot(volcano, aes(x = FC, y = t_test, q_value = q_value, modifications_in_proteins = modifications_in_proteins))+
            geom_point(aes(text = gene_name), size = 2, color = ifelse(volcano$mitocarta == TRUE, "black", "gray"))+
            scale_x_continuous(breaks=seq(min(floor(volcano$FC)), max(ceiling(volcano$FC)), 0.5), 
                               limits = c(min(floor(volcano$FC)), max(ceiling(volcano$FC))))+
            scale_y_continuous(breaks=seq(0, max(ceiling(volcano$t_test)), 0.5), limits = c(0, max(ceiling(volcano$t_test))))+
            theme_minimal()+
            geom_hline(aes(yintercept = 1.3), linetype = "dashed")+
            labs(title = paste(title, sep = ", "), x = paste("Log2 fold change", exp[1], "/", exp[2]), y = "-log10 p-value")+
            #labs(y = expression(paste("-Log"[10], " ", "p-value")), x = bquote(paste("Log"[2], " ", "fold change", " ", .(exp[1]), "/", .(exp[2]))))+
            ggtitle(title)+
            theme(axis.text=element_text(size=12, color = "black"),
                  axis.title=element_text(size=13,face="plain", color = "black"),
                  plot.title = element_text(size = 14, face = "plain", color = "black"),
                  legend.title = element_blank(),
                  axis.ticks = element_line(color = "black"),
                  panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
                  plot.background = element_rect(fill = "transparent", color = NA)) # bg of the plot
      
      ggplotly(q, tooltip = c("gene_name","FC", "t_test", "q_value","modifications_in_proteins"))
      
    }
    else{
      
      save_file_format = paste(title,".svg", sep = "")
      
      q = ggplot(volcano, aes(FC, t_test))+
            geom_point(size = 2, color = ifelse(volcano$t_test > 1.3, "red", "black"))+
            geom_text_repel(data = subset(volcano, (FC > 1.0 | FC < -1.0) & t_test > 2.0), aes(label = gene_name))+
            scale_x_continuous(breaks=seq(min(floor(volcano$FC)), max(ceiling(volcano$FC)), 0.5), 
                               limits = c(min(floor(volcano$FC)), max(ceiling(volcano$FC))))+
            scale_y_continuous(breaks=seq(0, max(ceiling(volcano$t_test)), 0.5), limits = c(0, max(ceiling(volcano$t_test))))+
            theme_minimal()+
            geom_hline(aes(yintercept = 1.3), linetype = "dashed")+
            labs(y = expression(paste("-Log"[10], " ", "p-value")), x = bquote(paste("Log"[2], " ", "fold change", " ", .(exp[1]), "/", .(exp[2]))))+
            ggtitle(title)+
            theme(axis.text=element_text(size=12, color = "black"),
                  axis.title=element_text(size=13,face="plain", color = "black"),
                  plot.title = element_text(size = 14, face = "plain", color = "black"),
                  legend.title = element_blank(),
                  axis.ticks = element_line(color = "black"),
                  panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
                  plot.background = element_rect(fill = "transparent", color = NA)) # bg of the plot
          #panel.grid.major = element_blank(), # get rid of major grid
          #panel.grid.minor = element_blank(), # get rid of minor grid
          #legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
          #legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel bg))
      if(print == T){
        ggsave(save_file_format, height = 7, width = 14 , units = "in", bg='transparent')
        while (!is.null(dev.list()))  dev.off()
      }
      else{
        print(q)
      }
      
      
    }
    
    
  }
  else if(type == "all_protein"){
    volcano <- df %>% 
      select(gene_name, mitocarta ,FC = contains("FC"), t.test = starts_with("p.value"), q_value = starts_with("q_value")) %>%
      mutate(t_test = -log10(t.test), .keep = c("unused")) %>%
      filter(!is.na(FC))
    
    
    
    
    if(plotly_plot == TRUE){
      
      q = ggplot(volcano, aes(x = FC, y = t_test, q_value = q_value))+
        geom_point(aes(text = gene_name), size = 2, color = ifelse(volcano$mitocarta == TRUE, "black", "gray"))+
        scale_x_continuous(breaks=seq(min(floor(volcano$FC)), max(ceiling(volcano$FC)), 0.5), 
                           limits = c(min(floor(volcano$FC)), max(ceiling(volcano$FC))))+
        scale_y_continuous(breaks=seq(0, max(ceiling(volcano$t_test)), 0.5), limits = c(0, max(ceiling(volcano$t_test))))+
        theme_minimal()+
        geom_hline(aes(yintercept = 1.3), linetype = "dashed")+
        labs(title = paste(title, sep = ", "), x = paste("Log2 fold change", exp[1], "/", exp[2]), y = "-log10 p-value")+
        #labs(y = expression(paste("-Log"[10], " ", "p-value")), x = bquote(paste("Log"[2], " ", "fold change", " ", .(exp[1]), "/", .(exp[2]))))+
        ggtitle(title)+
        theme(axis.text=element_text(size=12, color = "black"),
              axis.title=element_text(size=13,face="plain", color = "black"),
              plot.title = element_text(size = 14, face = "plain", color = "black"),
              legend.title = element_blank(),
              axis.ticks = element_line(color = "black"),
              panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
              plot.background = element_rect(fill = "transparent", color = NA)) # bg of the plot
      
      ggplotly(q, tooltip = c("gene_name", "FC", "t_test", "q_value"))
      
    }
    else{
      
      save_file_format = paste(title,".svg", sep = "")
      
      q = ggplot(volcano, aes(FC, t_test))+
        geom_point(size = 2, color = ifelse(volcano$t_test > 1.3, "red", "black"))+
        geom_text_repel(data = subset(volcano, (FC > 1.0 | FC < -1.0) & t_test > 2.0), aes(label = gene_name))+
        scale_x_continuous(breaks=seq(min(floor(volcano$FC)), max(ceiling(volcano$FC)), 0.5), 
                           limits = c(min(floor(volcano$FC)), max(ceiling(volcano$FC))))+
        scale_y_continuous(breaks=seq(0, max(ceiling(volcano$t_test)), 0.5), limits = c(0, max(ceiling(volcano$t_test))))+
        theme_minimal()+
        geom_hline(aes(yintercept = 1.3), linetype = "dashed")+
        labs(y = expression(paste("-Log"[10], " ", "p-value")), x = bquote(paste("Log"[2], " ", "fold change", " ", .(exp[1]), "/", .(exp[2]))))+
        ggtitle(title)+
        theme(axis.text=element_text(size=12, color = "black"),
              axis.title=element_text(size=13,face="plain", color = "black"),
              plot.title = element_text(size = 14, face = "plain", color = "black"),
              legend.title = element_blank(),
              axis.ticks = element_line(color = "black"),
              panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
              plot.background = element_rect(fill = "transparent", color = NA)) # bg of the plot
      #panel.grid.major = element_blank(), # get rid of major grid
      #panel.grid.minor = element_blank(), # get rid of minor grid
      #legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
      #legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel bg))
      if(print == T){
        ggsave(save_file_format, height = 7, width = 14 , units = "in", bg='transparent')
        while (!is.null(dev.list()))  dev.off()
      }
      else{
        print(q)
      }
      
      
    }
    
    
    
  }
  else if(type == "mito"){
    

    mito_reduced = mito_carta %>%
      select(GeneName, Symbol) %>%
      mutate(GeneName = as.character(GeneName))
    
    reduced = df %>%
      filter(!is.na(select(.,contains("t.test")))) %>%
      rename("GeneName" = entrez_gene_id) %>%
      semi_join(., mito_reduced, by = "GeneName", copy = F)
    
    
    volcano <- reduced %>% 
      select(gene_name, FC = contains("FC"), t.test = starts_with("t.test"), q_value = starts_with("q_value"))%>%
      mutate(t_test = -log10(t.test), .keep = c("unused")) %>%
      filter(!is.na(FC))
    
    
    save_file_format = paste(title,".svg", sep = "")
    
    if(plotly_plot == TRUE){
      library(plotly)
      q = ggplot(volcano, aes(FC, t_test, text = gene_name))+
            geom_point(size = 2, color = ifelse(volcano$t_test > 1.3, "red", "black"))+
            scale_x_continuous(breaks=seq(min(floor(volcano$FC)), max(ceiling(volcano$FC)), 0.5), 
                               limits = c(min(floor(volcano$FC)), max(ceiling(volcano$FC))))+
            scale_y_continuous(breaks=seq(0, max(ceiling(volcano$t_test)), 0.5), limits = c(0, max(ceiling(volcano$t_test))))+
            theme_minimal()+
            geom_hline(aes(yintercept = 1.3), linetype = "dashed")+
            labs(title = paste(title, sep = ", "), x = paste("Log2 fold change", exp[1], "/", exp[2]), y = "-log10 p-value")+
            #labs(y = expression(paste("-Log"[10], " ", "p-value")), x = bquote(paste("Log"[2], " ", "fold change", " ", .(exp[1]), "/", .(exp[2]))))+
            ggtitle(title)+
            theme(axis.text=element_text(size=12, color = "black"),
                  axis.title=element_text(size=13,face="plain", color = "black"),
                  plot.title = element_text(size = 14, face = "plain", color = "black"),
                  legend.title = element_blank(),
                  axis.ticks = element_line(color = "black"),
                  panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
                  plot.background = element_rect(fill = "transparent", color = NA)) # bg of the plot
      
      ggplotly(q)
      
      
    }
    else{
      q = ggplot(volcano, aes(FC, t_test))+
            geom_point(size = 2, color = ifelse(volcano$t_test > 1.3, "red", "black"))+
            scale_x_continuous(breaks=seq(min(floor(volcano$FC)), max(ceiling(volcano$FC)), 0.5), 
                               limits = c(min(floor(volcano$FC)), max(ceiling(volcano$FC))))+
            scale_y_continuous(breaks=seq(0, max(ceiling(volcano$t_test)), 0.5), limits = c(0, max(ceiling(volcano$t_test))))+
            theme_minimal()+
            geom_hline(aes(yintercept = 1.3), linetype = "dashed")+
            labs(y = expression(paste("-Log"[10], " ", "p-value")), x = bquote(paste("Log"[2], " ", "fold change", " ", .(exp[1]), "/", .(exp[2]))))+
            geom_text_repel(data = subset(volcano, (FC > 0.5 | FC < -0.5) | t_test > 2.0), aes(label = gene_name))+
            ggtitle(title)+
            theme(axis.text=element_text(size=12, color = "black"),
                  axis.title=element_text(size=13,face="plain", color = "black"),
                  plot.title = element_text(size = 14, face = "plain", color = "black"),
                  legend.title = element_blank(),
                  axis.ticks = element_line(color = "black"),
                  panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
                  plot.background = element_rect(fill = "transparent", color = NA)) # bg of the plot
          #panel.grid.major = element_blank(), # get rid of major grid
          #panel.grid.minor = element_blank(), # get rid of minor grid
          #legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
          #legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel bg))
      if(print == T){
        ggsave(save_file_format, height = 7, width = 14 , units = "in", bg='transparent')
        while (!is.null(dev.list()))  dev.off()
      }else{
        print(q)
      }
      
    }
    
      
    }
  
  
}

v_plot = volcano_plot(final_proteins, "all_protein", "all proteins, LIVER relative occupancy TRF/control", c("control", "TRF"), F, T)
v_plot

saveWidget(v_plot, "all peptides LIVER relative occupancy TRF over control.html", selfcontained = T, libdir = "lib")
# basic values for proteins IDs

# volcano plot for limma data

volcano_plot = function(df, type, title, exp, print, plotly_plot){
  options(ggrepel.max.overlaps = 100) 
  if(type == "all_peptide"){
    
    volcano <- df %>% 
      select(gene_name, modifications_in_proteins, mitocarta, FC = contains("FC"), t.test = starts_with("p.value"), q_value = starts_with("q_value")) %>%
      mutate(t_test = -log10(t.test), .keep = c("unused")) %>%
      filter(!is.na(FC))
    
    
    
    
    if(plotly_plot == TRUE){
      
      q = ggplot(volcano, aes(x = FC, y = t_test, q_value = q_value, modifications_in_proteins = modifications_in_proteins))+
        geom_point(aes(text = gene_name), size = 2, color = ifelse(volcano$mitocarta == TRUE, "black", "gray"))+
        scale_x_continuous(breaks=seq(min(floor(volcano$FC)), max(ceiling(volcano$FC)), 0.5), 
                           limits = c(min(floor(volcano$FC)), max(ceiling(volcano$FC))))+
        scale_y_continuous(breaks=seq(0, max(ceiling(volcano$t_test)), 0.5), limits = c(0, max(ceiling(volcano$t_test))))+
        theme_minimal()+
        geom_hline(aes(yintercept = 1.3), linetype = "dashed")+
        labs(title = paste(title, sep = ", "), x = paste("Log2 fold change", exp[1], "/", exp[2]), y = "-log10 p-value")+
        #labs(y = expression(paste("-Log"[10], " ", "p-value")), x = bquote(paste("Log"[2], " ", "fold change", " ", .(exp[1]), "/", .(exp[2]))))+
        ggtitle(title)+
        theme(axis.text=element_text(size=12, color = "black"),
              axis.title=element_text(size=13,face="plain", color = "black"),
              plot.title = element_text(size = 14, face = "plain", color = "black"),
              legend.title = element_blank(),
              axis.ticks = element_line(color = "black"),
              panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
              plot.background = element_rect(fill = "transparent", color = NA)) # bg of the plot
      
      ggplotly(q, tooltip = c("gene_name","FC", "t_test", "q_value","modifications_in_proteins"))
      
    }
    else{
      
      save_file_format = paste(title,".svg", sep = "")
      
      q = ggplot(volcano, aes(FC, t_test))+
        geom_point(size = 2, color = ifelse(volcano$t_test > 1.3, "red", "black"))+
        geom_text_repel(data = subset(volcano, (FC > 1.0 | FC < -1.0) & t_test > 2.0), aes(label = gene_name))+
        scale_x_continuous(breaks=seq(min(floor(volcano$FC)), max(ceiling(volcano$FC)), 0.5), 
                           limits = c(min(floor(volcano$FC)), max(ceiling(volcano$FC))))+
        scale_y_continuous(breaks=seq(0, max(ceiling(volcano$t_test)), 0.5), limits = c(0, max(ceiling(volcano$t_test))))+
        theme_minimal()+
        geom_hline(aes(yintercept = 1.3), linetype = "dashed")+
        labs(y = expression(paste("-Log"[10], " ", "p-value")), x = bquote(paste("Log"[2], " ", "fold change", " ", .(exp[1]), "/", .(exp[2]))))+
        ggtitle(title)+
        theme(axis.text=element_text(size=12, color = "black"),
              axis.title=element_text(size=13,face="plain", color = "black"),
              plot.title = element_text(size = 14, face = "plain", color = "black"),
              legend.title = element_blank(),
              axis.ticks = element_line(color = "black"),
              panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
              plot.background = element_rect(fill = "transparent", color = NA)) # bg of the plot
      #panel.grid.major = element_blank(), # get rid of major grid
      #panel.grid.minor = element_blank(), # get rid of minor grid
      #legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
      #legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel bg))
      if(print == T){
        ggsave(save_file_format, height = 7, width = 14 , units = "in", bg='transparent')
        while (!is.null(dev.list()))  dev.off()
      }
      else{
        print(q)
      }
      
      
    }
    
    
  }
  else if(type == "all_protein"){
    volcano <- df %>% 
      select(gene_name, mitocarta ,FC = contains("FC"), t.test = starts_with("p.value")) %>%
      mutate(t_test = -log10(t.test), .keep = c("unused")) %>%
      filter(!is.na(FC))
    
    
    
    
    if(plotly_plot == TRUE){
      
      q = ggplot(volcano, aes(x = FC, y = t_test))+
        geom_point(aes(text = gene_name), size = 2, color = ifelse(volcano$mitocarta == TRUE, "black", "gray"))+
        scale_x_continuous(breaks=seq(min(floor(volcano$FC)), max(ceiling(volcano$FC)), 0.5), 
                           limits = c(min(floor(volcano$FC)), max(ceiling(volcano$FC))))+
        scale_y_continuous(breaks=seq(0, max(ceiling(volcano$t_test)), 0.5), limits = c(0, max(ceiling(volcano$t_test))))+
        theme_minimal()+
        geom_hline(aes(yintercept = 1.3), linetype = "dashed")+
        labs(title = paste(title, sep = ", "), x = paste("Log2 fold change", exp[1], "/", exp[2]), y = "-log10 p-value")+
        #labs(y = expression(paste("-Log"[10], " ", "p-value")), x = bquote(paste("Log"[2], " ", "fold change", " ", .(exp[1]), "/", .(exp[2]))))+
        ggtitle(title)+
        theme(axis.text=element_text(size=12, color = "black"),
              axis.title=element_text(size=13,face="plain", color = "black"),
              plot.title = element_text(size = 14, face = "plain", color = "black"),
              legend.title = element_blank(),
              axis.ticks = element_line(color = "black"),
              panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
              plot.background = element_rect(fill = "transparent", color = NA)) # bg of the plot
      
      ggplotly(q, tooltip = c("gene_name", "FC", "t_test"))
      
    }
    else{
      
      save_file_format = paste(title,".svg", sep = "")
      
      q = ggplot(volcano, aes(FC, t_test))+
        geom_point(size = 2, color = ifelse(volcano$t_test > 1.3, "red", "black"))+
        geom_text_repel(data = subset(volcano, (FC > 1.0 | FC < -1.0) & t_test > 2.0), aes(label = gene_name))+
        scale_x_continuous(breaks=seq(min(floor(volcano$FC)), max(ceiling(volcano$FC)), 0.5), 
                           limits = c(min(floor(volcano$FC)), max(ceiling(volcano$FC))))+
        scale_y_continuous(breaks=seq(0, max(ceiling(volcano$t_test)), 0.5), limits = c(0, max(ceiling(volcano$t_test))))+
        theme_minimal()+
        geom_hline(aes(yintercept = 1.3), linetype = "dashed")+
        labs(y = expression(paste("-Log"[10], " ", "p-value")), x = bquote(paste("Log"[2], " ", "fold change", " ", .(exp[1]), "/", .(exp[2]))))+
        ggtitle(title)+
        theme(axis.text=element_text(size=12, color = "black"),
              axis.title=element_text(size=13,face="plain", color = "black"),
              plot.title = element_text(size = 14, face = "plain", color = "black"),
              legend.title = element_blank(),
              axis.ticks = element_line(color = "black"),
              panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
              plot.background = element_rect(fill = "transparent", color = NA)) # bg of the plot
      #panel.grid.major = element_blank(), # get rid of major grid
      #panel.grid.minor = element_blank(), # get rid of minor grid
      #legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
      #legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel bg))
      if(print == T){
        ggsave(save_file_format, height = 7, width = 14 , units = "in", bg='transparent')
        while (!is.null(dev.list()))  dev.off()
      }
      else{
        print(q)
      }
      
      
    }
    
    
    
  }
  else if(type == "mito"){
    
    
    mito_reduced = mito_carta %>%
      select(GeneName, Symbol) %>%
      mutate(GeneName = as.character(GeneName))
    
    reduced = df %>%
      filter(!is.na(select(.,contains("t.test")))) %>%
      rename("GeneName" = entrez_gene_id) %>%
      semi_join(., mito_reduced, by = "GeneName", copy = F)
    
    
    volcano <- reduced %>% 
      select(gene_name, FC = contains("FC"), t.test = starts_with("t.test"), q_value = starts_with("q_value"))%>%
      mutate(t_test = -log10(t.test), .keep = c("unused")) %>%
      filter(!is.na(FC))
    
    
    save_file_format = paste(title,".svg", sep = "")
    
    if(plotly_plot == TRUE){
      library(plotly)
      q = ggplot(volcano, aes(FC, t_test, text = gene_name))+
        geom_point(size = 2, color = ifelse(volcano$t_test > 1.3, "red", "black"))+
        scale_x_continuous(breaks=seq(min(floor(volcano$FC)), max(ceiling(volcano$FC)), 0.5), 
                           limits = c(min(floor(volcano$FC)), max(ceiling(volcano$FC))))+
        scale_y_continuous(breaks=seq(0, max(ceiling(volcano$t_test)), 0.5), limits = c(0, max(ceiling(volcano$t_test))))+
        theme_minimal()+
        geom_hline(aes(yintercept = 1.3), linetype = "dashed")+
        labs(title = paste(title, sep = ", "), x = paste("Log2 fold change", exp[1], "/", exp[2]), y = "-log10 p-value")+
        #labs(y = expression(paste("-Log"[10], " ", "p-value")), x = bquote(paste("Log"[2], " ", "fold change", " ", .(exp[1]), "/", .(exp[2]))))+
        ggtitle(title)+
        theme(axis.text=element_text(size=12, color = "black"),
              axis.title=element_text(size=13,face="plain", color = "black"),
              plot.title = element_text(size = 14, face = "plain", color = "black"),
              legend.title = element_blank(),
              axis.ticks = element_line(color = "black"),
              panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
              plot.background = element_rect(fill = "transparent", color = NA)) # bg of the plot
      
      ggplotly(q)
      
      
    }
    else{
      q = ggplot(volcano, aes(FC, t_test))+
        geom_point(size = 2, color = ifelse(volcano$t_test > 1.3, "red", "black"))+
        scale_x_continuous(breaks=seq(min(floor(volcano$FC)), max(ceiling(volcano$FC)), 0.5), 
                           limits = c(min(floor(volcano$FC)), max(ceiling(volcano$FC))))+
        scale_y_continuous(breaks=seq(0, max(ceiling(volcano$t_test)), 0.5), limits = c(0, max(ceiling(volcano$t_test))))+
        theme_minimal()+
        geom_hline(aes(yintercept = 1.3), linetype = "dashed")+
        labs(y = expression(paste("-Log"[10], " ", "p-value")), x = bquote(paste("Log"[2], " ", "fold change", " ", .(exp[1]), "/", .(exp[2]))))+
        geom_text_repel(data = subset(volcano, (FC > 0.5 | FC < -0.5) | t_test > 2.0), aes(label = gene_name))+
        ggtitle(title)+
        theme(axis.text=element_text(size=12, color = "black"),
              axis.title=element_text(size=13,face="plain", color = "black"),
              plot.title = element_text(size = 14, face = "plain", color = "black"),
              legend.title = element_blank(),
              axis.ticks = element_line(color = "black"),
              panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
              plot.background = element_rect(fill = "transparent", color = NA)) # bg of the plot
      #panel.grid.major = element_blank(), # get rid of major grid
      #panel.grid.minor = element_blank(), # get rid of minor grid
      #legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
      #legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel bg))
      if(print == T){
        ggsave(save_file_format, height = 7, width = 14 , units = "in", bg='transparent')
        while (!is.null(dev.list()))  dev.off()
      }else{
        print(q)
      }
      
    }
    
    
  }
  
  
}

v_plot = volcano_plot(final_proteins, "all_protein", "all proteins, LIVER relative occupancy TRF/control", c("control", "TRF"), F, T)
v_plot

saveWidget(v_plot, "all peptides LIVER relative occupancy TRF over control.html", selfcontained = T, libdir = "lib")



barplot_proteins = function(df, names.arg, main){
  

  df1 <- rename(df, "GeneName" = "entrez_gene_id")

  # semi_join appears to be the correct way of filtering out the mito proteins
  mito_reduced = df1 %>%
    semi_join(., mito_carta, by = "GeneName") %>%
    filter(!is.na(select(.,contains("t.test"))))

  # this will greatly reduce the number of identified peptides,if the peptide doesn't have a p.value
  df_reduced = df1 %>%
    filter(!is.na(select(.,contains("t.test"))))
  
  mito_sig = mito_reduced %>%
    filter(select(.,starts_with("t.test")) < 0.05)
  
  y = data.frame(variables = c(names.arg),
                 values = c(nrow(df_reduced), nrow(mito_reduced), nrow(mito_sig)))
  
  ggplot(y, aes(x = variables, y = values, label = values))+
    geom_col()+
    theme_minimal()+
    geom_text(aes(label = values), vjust = -0.5)+
    labs(title = main, x = "", y = "# peptides/proteins")+
    theme(axis.text=element_text(size=12, color = "black"),
          axis.title=element_text(size=13,face="plain", color = "black"),
          plot.title = element_text(size = 14, face = "plain", color = "black"),
          legend.title = element_blank(),
          axis.ticks = element_line(color = "black"),
          panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA)) # bg of the plot
  
  
  
  
  # out of all the mito specific proteins how many significantly changed
  

  
  
  
  
  
  
}

barplot_proteins(final_peptides, c("all_proteins", "mito", "mito p<0.05"), "# quantified proteins in non pyruvate set")

####

int1 = final_proteins %>%
  select(pd_ratio = "pd ratio", FC = starts_with("FC")) %>%
  mutate(pd_ratio = log2(pd_ratio)) %>%
  ggplot(., aes(x = pd_ratio, y = FC))+
  geom_point()+
  labs(title = "phos peptides occupancy FC compared to PD ratio")
int1

####

# quantified phospho sites

barplot_phos = function(df, sample, names.arg, main){
  
  x = filter(final_peptides, !is.na(select(.,contains(ID_names))))
  
  y = data.frame(variables = c("all", "quant"),
                 values = c(nrow(df), nrow(x)))
  
  ggplot(y, aes(x = variables, y = values, label = values))+
    geom_col()+
    theme_minimal()+
    geom_text(aes(label = values), vjust = -0.5)+
    labs(title = "# quantified phos sites", x = "", y = "phos_peptides")+
    theme(axis.text=element_text(size=12, color = "black"),
          axis.title=element_text(size=13,face="plain", color = "black"),
          #plot.title = element_text(size = 16, face = "plain", color = "black"),
          legend.title = element_blank(),
          axis.ticks = element_line(color = "black"),
          panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA)) # bg of the plot
    
  
}

multi_plot = function(df, sample, title, exp, volcano_xlim, volcano_ylim, print, leg_height, hist_legend){
  
  mito_carta = fread("Human.MitoCarta3.0.txt", header = T)
  
  mito_reduced = mito_carta %>%
    select(HumanGeneID, Symbol)
  
  reduced = df %>%
    filter(!is.na(select(.,contains("t.test")))) %>%
    rename("HumanGeneID" = entrez_gene_id) %>%
    inner_join(., mito_reduced, by = "HumanGeneID", copy = F)
  
  
  volcano <- reduced %>% 
    select(Symbol, FC = contains("FC"), t.test = starts_with("t.test"), q_value = starts_with("q_value"))%>%
    mutate(t_test = -log10(t.test), .keep = c("unused")) %>%
    filter(!is.na(FC))
  
  
  save_file_format = paste(title,".svg", sep = "")
  
  q = ggplot(volcano, aes(FC, t_test))+
    geom_point(size = 2, color = ifelse(volcano$t_test > 1.3, "red", "black"))+
    scale_x_continuous(breaks=seq(-8, 8, 0.5), limits = volcano_xlim)+
    scale_y_continuous(breaks=seq(-8, 8, 0.5), limits = volcano_ylim)+
    theme_minimal()+
    geom_hline(aes(yintercept = 1.3), linetype = "dashed")+
    labs(y = expression(paste("-Log"[10], " ", "p-value")), x = bquote(paste("Log"[2], " ", "fold change", " ", .(exp[1]), "/", .(exp[2]))))+
    geom_text_repel(data = subset(volcano, (FC > 0.5 | FC < -0.5) | t_test > 2.0), aes(label = Symbol))+
    ggtitle(paste(sample, title, sep = ", "))+
    theme(axis.text=element_text(size=12, color = "black"),
          axis.title=element_text(size=13,face="plain", color = "black"),
          plot.title = element_text(size = 14, face = "plain", color = "black"),
          legend.title = element_blank(),
          axis.ticks = element_line(color = "black"),
          panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA)) # bg of the plot
  #panel.grid.major = element_blank(), # get rid of major grid
  #panel.grid.minor = element_blank(), # get rid of minor grid
  #legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
  #legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel bg))
  if(print == T){
    ggsave(save_file_format, height = 7, width = 14 , units = "in", bg='transparent')
    while (!is.null(dev.list()))  dev.off()
  }else{
    print(q)
  }
  
  
  
  
  
  par(mfrow = c(1,1))
  hist(pull(select(df, starts_with("t.test"))), col = rgb(1,0,0,0.5), 
       main = paste(sample, title, sep = ", "), xlab = "p.value", breaks = 20)
  hist(pull(select(mito_reduced, starts_with("t.test"))), col = rgb(0,0,1,0.5), add=T, breaks = 20)
  legend(0.8, leg_height, hist_legend, lty = c(1,1),lwd = c(2.5, 2.5), 
         col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))
  
  
}

multi_plot(final_peptides, "proteins", "PDK4_fa/BGAL_fa, with pyr", c("PDK4_fa", "Bgal_fa"), c(-2,2), c(0,3),
           F, 1000, c("all", "mito"))


# correlation plots
library(psych)

pairs.panels(final_protein[,16:18], method = 'pearson', density = T)

plot(x=final_protein$mean_pdk4_fa_nopy, y = final_protein$mean_bgal_fa_nopy)

ggplot(final_protein, aes(x=mean_pdk4_fa_nopy, y = mean_bgal_fa_nopy))+
  geom_point()


# histogram

histo = function(df, sample, main, legend, height){
  
  mito_carta = fread("Human.MitoCarta3.0.txt", header = T)
  
  if(sample == "peptide"){
    df <- rename(df, "accession" = "1st_master_protein_accession")
  }
  
  
  mito_reduced = mito_carta %>%
    select(., accession, Symbol) %>%
    right_join(df,.,by = "accession")
  
  
  par(mfrow = c(1,1))
  hist(pull(select(df, starts_with("t.test"))), col = rgb(1,0,0,0.5), 
       main = paste(sample, main, sep = ", "), xlab = "p.value", breaks = 20)
  hist(pull(select(mito_reduced, starts_with("t.test"))), col = rgb(0,0,1,0.5), add=T, breaks = 20)
  legend(0.8, height, legend, lty = c(1,1),lwd = c(2.5, 2.5), 
         col = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)))
  
  
}

histo(final_peptide, "peptide", "p.values, PDK4_fa/Bgal_fa, with pyr", c("all", "mito"), 400)

# heatmap
library(pheatmap)

mito_carta = fread("Human.MitoCarta3.0.txt", header = T)

final_peptides <- rename(final_peptides, "HumanGeneID" = entrez_gene_id)

mito_carta = mito_carta %>%
  mutate(GeneName = as.character(GeneName))


heat = final_proteins %>%
  rename(., "GeneName" = entrez_gene_id) %>%
  #filter(!is.na(select(.,contains("t.test")))) %>%
  semi_join(., mito_carta, by = "GeneName") %>%
  select(starts_with("abundance")) %>%
  as.matrix()

pheatmap(heat, show_rownames = F, clustering_method = "ward.D2")

# dendogram

dend = dist(t(heat))
hc = hclust(dend)
plot(hc)





# below is looking at how multiple technical replicate injections effects phospho peptide abundance in the TRF liver phos dataset
df = fread("2021-03-02_AW_TRF_LIVER_inj_phos_PeptideGroups.txt")

int1 = df %>%
  select("Master Protein Accessions", starts_with("Found in Fraction: F4")) %>%
  select(1:19)

new_names = paste0(rep(c("inj1_", "inj2_", "inj3_"), each = 6), rep("input"), rep(1:6, 3))

colnames(int1) = c("accession", new_names)

int2 = as.data.frame(int1) %>%
  pivot_longer(cols = contains("inj")) %>%
  separate(name, c("injection_number", "input_number")) %>%
  mutate(value = as.factor(value)) %>%
  group_by(injection_number,input_number) %>%
  count(value) %>%
  ungroup() %>%
  filter(value == "High") %>%
  ggplot(., aes(x = input_number ,y = n, fill = injection_number))+
  geom_col(position = "dodge")+
  labs(title = "TRF liver phos dataset, number of phos proteins identified from inputs 1-6 and injections 1-3")

int2


replicate_injections = function(df){
  
  int1 = df %>%
    select(accession, starts_with("Found in Fraction")) 
  %>%
    select(1:19)
  
  new_names = paste0(rep(c("inj1_", "inj2_", "inj3_"), each = 6), rep("input"), rep(1:6, 3))
  
  colnames(int1) = c("Accession", new_names)
  
  inputs = paste0(rep("input"), rep(1:6))
  
  test = map(inputs, function(x){
    int3 = int1 %>%
      select(accession, contains(x)) %>%
      mutate_all(as.factor) 
    
    z = apply(int3, 1, function(x) length(which(x == "High")))
    
    print(table(z))
  })
  
  table(test)
  int3 = int1 %>%
    select(gene_name, contains("input1")) %>%
    mutate_all(as.factor) 
  
  z = apply(int3, 1, function(x) length(which(x == "High")))
  
  int4 = cbind(int3, z)
  
  
}


int3 = int1 %>%
  select(accession, contains("input1")) %>%
  mutate_all(as.factor) 

z = apply(int3, 1, function(x) length(which(x == "High")))

int4 = cbind(int3, z)

ggplot(int4, aes(x = gene_name, y = z)) +
  geom_point()

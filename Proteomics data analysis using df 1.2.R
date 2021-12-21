
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
# reduces the number of columns. Keeps the columns set in the description cols list, and abundance quantification values
# There are 2 different methods that can be used, the first is method = "protein" for protein data and
# the other is method = "peptide-phos" for phospho data. 


# df1 
reduce_cols = function(df1, method, fraction){
  if(method == "protein"){
    
    df1 = clean_names(df1)
    
    reduced_protein <- df1 %>%
      setnames(names(select(., contains("q_value"))), "exp_q_value") %>%
      select(., all_of(description_cols[[1]]), starts_with(paste0("abundance_", paste0(fraction, "_")))) 

   reduced_protein[reduced_protein == 0] = NA
    
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
    df2 = clean_names(proteins)
    
    reduced_peptides <- df1 %>%
      select(., all_of(description_cols[[2]]), starts_with(paste0("abundance_", paste0(fraction, "_"))))

    
    reduced_peptides[reduced_peptides == 0] = NA
    
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
}

reduced_cols_protein = reduce_cols(proteins, method = "protein", "F2")

reduced_cols_peptide = reduce_cols(peptides, method = "peptide-phos", "F6")

# reducing rows of dataframe

# df is the data.frame that will be used in the analysis
# method = "tmt_protein" is for only keeping IsMasterProtein, and qvalue < 0.01 and removes proteins that are
# only in 1/2 the samples
# method = "tmt_peptide_phos" filters data that has "Phospho" in modification column

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

reduced_rows_protein = reduce_rows(reduced_cols_protein, method = "tmt_protein")

reduced_rows_peptide = reduce_rows(reduced_cols_peptide, method = "tmt_peptide_phos")


# normalization
# Normalizes the data to the ratio of the  colsums / mean(columns)
# the df1 argument is ideally the reduced_rows of either the protein or peptide data
# the df2 argument is the 

# change the fraction to the fraction in which the input is in

normalize_df = function(df1, method, fraction, show_norm_boxplot) {
  
  if(method == "protein"){
    
    # normalize the protein abundance between sample by using the ratio (sum/avg) of the peptide data
    # make sure to use the original peptide data and not a df that was already reduced
    
    norm_names_df1 = grep("abundance", names(df1), value = TRUE)

    reduced_pep = peptides %>%
      select(starts_with(paste("Abundance:", paste(fraction, ":", sep = ''), sep = ' '))) %>%
      mutate_if(is.character, as.numeric)

    sum = colSums(reduced_pep, na.rm = T)
    avg = mean(sum)
    
    ratio = sum/avg
    
    show(paste("Ratio of peptide"))
    show(round(ratio, digits = 4))
    
    axis =  c("KO1", "KO2", "KO3", "KO4", "KO5", "Pool",
              "WT1", "WT2", "WT3", "WT4", "WT5")
    
    bar = barplot(ratio, ylab = "ratio", ylim = c(0,ceiling(max(ratio,na.rm=T))), 
                  main = paste("Ratio of peptide input"), names.arg = axis)
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
        dplyr::select(.,starts_with("abundance"))
      
      colnames(before) = axis
    
      par(mfrow = c(1,2))
      
      plotDensities(log2(before), main = "Before Normalization")
      
      box1 = boxplot(log2(before), col = "dodgerblue", main = "Before Normalization", ylab = "Protein abundance")
      
      colnames(logged_norm) = axis
      
      plotDensities(logged_norm, main = "After Loading Normalization")
      
      box2 = boxplot(logged_norm, col = "springgreen", main = "After Loading Normalization", ylab = "Protein abundance")

      colnames(norm_to_avg) = axis
      
      box3 = boxplot(norm_to_avg, col = "firebrick1", main = "After Normalizing To The Mean", ylab = "Protein abundance")

      plotDensities(logged_norm, main = "After Loading Normalization")
      
      
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

normalized_protein = normalize_df(df1 = reduced_rows_protein, "protein", "F2", T)

normalized_peptide = normalize_df(df1 = reduced_rows_peptide, "peptide", "F2")


# quantification
samples_to_compare = c("WT", "KO")
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
      
      control_1 = means1[paste0(z, "_mean_", control)]
      variables_1 = means1[paste0(z, "_mean_", variables)]
      

      fc = map2(control_1, variables_1, function(a,b){
        int1 = b-a
      })
      
      names(fc) = c(paste0(z, "_FC_", variables, "_over_", control))
      
      
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
      
      final = as.data.frame(cbind(means1, sd, as.data.frame(fc), stat_results))

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

final_proteins = new_final_df(df = normalized_protein, method = list("abundance"), type = 2)

final_peptides = new_final_df(df = normalized_peptide, method = list("abundance", "relative_occupancy"), type = 2)


# scatter plot of phosphopeptide abundance vs relative occupancy

abundance = final_peptides %>%
  dplyr::select(., -contains("relative_occupancy")) %>%
  dplyr::select(., FC_abund = contains("FC"))

rel_occ = final_peptides %>%
  dplyr::select(., -contains("abundance")) %>%
  dplyr::select(., FC_rel_occ = contains("FC"))

comb = cbind(abundance, rel_occ)

plot(comb$FC_abund, comb$FC_rel_occ, pch = 20, xlab = "Phosphopeptide abundance", ylab = "Relative occupancy", 
     xlim = c(-2,4), ylim = c(-2,4), col = "gray31",
     main = "Correlation between phosphopeptide abundance and relative occupancy")
abline(a=0, b=1, lwd = 2)





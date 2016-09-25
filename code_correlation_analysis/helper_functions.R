get_lncRNA_mRNA_pairs = function(originalCorrMatrix, corrThreshold){
  correlation_pairs = which(originalCorrMatrix > corrThreshold, arr.ind = TRUE)
  lncRNA = rownames(originalCorrMatrix)[correlation_pairs[,1]]
  mRNA = colnames(originalCorrMatrix)[correlation_pairs[,2]]
  dataframe = as.data.frame(cbind(correlation_pairs, lncRNA, mRNA))
  colnames(dataframe) = c("lncRNA_index", "mRNA_index","lncRNA", "mRNA")
  dataframe$lncRNA_index = as.numeric(as.character(dataframe$lncRNA_index))
  dataframe$mRNA_index = as.numeric(as.character(dataframe$mRNA_index))
  # correlation_vector = normal_lncRNA_mRNA_corr_matrix[dataframe$lncRNA_index,dataframe$mRNA_index]
  dataframe = cbind(dataframe, originalCorrMatrix[which(originalCorrMatrix > corrThreshold)])
  colnames(dataframe) = c("lncRNA_index", "mRNA_index","lncRNA", "mRNA", "corr")
  dataframe = dataframe[order(-dataframe$corr),]
  return(dataframe)
}

get_encRNA = function(matrix, miRNA_mRNA_corr, miRNA_lncRNA_corr, lncRNA_mRNA_corr, threshold){
  triple = which(matrix > threshold, arr.ind = TRUE)
  encRNAs = rownames(matrix)[triple[,1]]
  miRNA = colnames(matrix)[triple[,2]]
  dataframe = as.data.frame(cbind(triple, encRNAs, miRNA))
  colnames(dataframe) = c("encRNA_pair_index", "miRNA_index","encRNA_pair", "miRNA")
  dataframe$encRNA_pair_index = as.numeric(as.character(dataframe$encRNA_pair_index))
  dataframe$miRNA_index = as.numeric(as.character(dataframe$miRNA_index))
  dataframe = cbind(dataframe, matrix[which(matrix > threshold)])
  colnames(dataframe) = c("encRNA_pair_index", "miRNA_index","encRNA_pair", "miRNA", "sensitivity")
  dataframe = dataframe[order(-dataframe$sensitivity),]
  # dataframe$lncRNA = substr(dataframe$encRNA_pair, start = 1, stop = 17)
  # dataframe$mRNA = substr(dataframe$encRNA_pair, start = 19, stop = nchar(as.character(dataframe$encRNA_pair)))
  
  pos = regexpr("-",as.character(dataframe$encRNA_pair))
  dataframe$lncRNA = substr(as.character(dataframe$encRNA_pair), start = 1, stop = pos - 1)
  dataframe$mRNA = substr(as.character(dataframe$encRNA_pair), start = pos + 1, stop = nchar(as.character(dataframe$encRNA_pair))) 
  
  dataframe$lncRNA_miRNA_corr = rep(NA, nrow(dataframe))
  dataframe$mRNA_miRNA_corr = rep(NA, nrow(dataframe))
  
  dataframe$lncRNA_miRNA_corr = sapply(1:nrow(dataframe), function(index){
    miRNA_lncRNA_corr[as.character(dataframe$miRNA[index]),as.character(dataframe$lncRNA[index])]
  })
  
  dataframe$mRNA_miRNA_corr = sapply(1:nrow(dataframe), function(index){
    miRNA_mRNA_corr[as.character(dataframe$miRNA[index]),as.character(dataframe$mRNA[index])]
  })
  
  dataframe$lncRNA_mRNA_corr = sapply(1:nrow(dataframe), function(index){
    lncRNA_mRNA_corr[as.character(dataframe$lncRNA[index]),as.character(dataframe$mRNA[index])]
  })
  
  dataframe$encRNA_triple = paste(dataframe$encRNA_pair, dataframe$miRNA, sep = "-")
  dataframe$lncRNA_miRNA_pair = paste(dataframe$lncRNA, dataframe$miRNA, sep = "-")
  dataframe$mRNA_miRNA_pair = paste(dataframe$mRNA, dataframe$miRNA, sep = "-")
  dataframe = dataframe[,c("lncRNA", "mRNA", "miRNA", "sensitivity", 
                           "lncRNA_mRNA_corr", "lncRNA_miRNA_corr", "mRNA_miRNA_corr",
                           "encRNA_pair", "encRNA_triple", 
                           "lncRNA_miRNA_pair" , "mRNA_miRNA_pair", 
                           "encRNA_pair_index", "miRNA_index")]
  return(dataframe)
}



# for each lncRNA ensemble ids, find potential miRNA interaction
getMiRNAs = function(miRNA_family = NULL){
  if (is.null(miRNA_family)){
    print("must provide miRNA_family")
    break;
  }
  part1 = substr(miRNA_family, start = 1, stop = 3)
  part1 = paste("hsa",tolower(part1),sep = "-")
  
  # substring anything after miR till the end, divided by "/"
  part2 = substr(miRNA_family,start = 5, stop = nchar(miRNA_family))
  part2 = unlist(strsplit(x = part2, split = "/"))
  
  # foreach element, remove 3p and 5p parts
  part2 = gsub("-3p","",part2)
  part2 = gsub("-5p","",part2)
  
  # return individual mircRNA, 
  # example: 106abc will be disconstructed into 106a, 106b, 106c
  part2 =  sapply(part2, function(element){
    if (grepl("\\D",element)){
      digit_part = gsub(pattern = "\\D", replacement = "", x = element)
      character_parts = gsub(pattern = "\\d", replacement = "", x = element)
      character_parts = unlist(strsplit(x = character_parts,split = ""))
      returned_value = paste(digit_part, character_parts,sep = "")
    }else{
      element
    }
  })
  part2 = unname(unlist(part2))
  return(paste(part1,part2,sep="-"))
}

get_putative_lncRNA_miRNA = function(dataframe){
  require(rlist)
  l = list()
  #i = 1;
  apply(dataframe, 1, function(r){
    k = getMiRNAs(r[2])
    l <<- list.append(l, k)
    #print(i); i <<- i + 1;
  })
  names(l) = dataframe$gene_id
  df = reshape2::melt(l)
  colnames(df) = c("miRNA", "putative_lncRNAs")
  df$lncRNA_miRNA_pair = paste(df$putative_lncRNAs, df$miRNA, sep = "-")
  return(df)
}

get_putative_lncRNA_miRNA_2 = function(dataframe){
  require(rlist)
  l = list()
  i = 1;
  apply(dataframe, 1, function(r){
    k = getMiRNAs(r[2])
    l <<- list.append(l, k)
    print(i); i <<- i + 1;
  })
  names(l) = dataframe$gene_id
  df = data.frame(
    miRNA = unlist(l),
    putative_lncRNAs = rep(names(l), lapply(l, length))
  )
  colnames(df) = c("miRNA", "putative_lncRNAs")
  df$lncRNA_miRNA_pair = paste(df$putative_lncRNAs, df$miRNA, sep = "-")
  return(df)
}


get_matched_enRNA_sensitivity_with_putative_binding = function(encRNA_sensitivity){
  load("mircode_objects.rda")
  lncRNAs_overlapped = intersect(unique(mircode_lncRNA$gene_id), unique(encRNA_sensitivity$lncRNA))
  # subset the encRNA_sensivivity to include only lncRNAs matches lncRNAs in miRcode 
  encRNA_sensitivity_subset1 = encRNA_sensitivity[which(encRNA_sensitivity$lncRNA %in% lncRNAs_overlapped),] 
  # similarly, subset the mircode_lncRNA to include only lncRNAs matches lncRNAs in encRNA_sensivivity
  mircode_lncRNA_subset1 = mircode_lncRNA[which(mircode_lncRNA$gene_id %in% lncRNAs_overlapped), 
                                          c("gene_id","microrna")]
  mircode_lncRNA_subset1 = get_putative_lncRNA_miRNA(mircode_lncRNA_subset1) # divide miRNAs familily into individual miRNAs
  # now, subset encRNA_sensivivity_subset1 to include only the lncRNA-miRNA pairs which also shows up in mircode_lncRNA_subset1
  # length(intersect(unique(encRNA_sensitivity_subset1$lncRNA_miRNA_pair), unique(mircode_lncRNA_subset1$lncRNA_miRNA_pair)))
  intersected_lncRNA_miRNA_pairs = intersect(unique(encRNA_sensitivity_subset1$lncRNA_miRNA_pair), unique(mircode_lncRNA_subset1$lncRNA_miRNA_pair))
  encRNA_sensitivity_subset2 = encRNA_sensitivity_subset1[which(encRNA_sensitivity_subset1$lncRNA_miRNA_pair %in% intersected_lncRNA_miRNA_pairs),]
  # now, we have already found all lncRNA_miRNA pairs in the sensitivity matrix that are also included in miRcode, thus the duty of miRcode is done now 
  # next, we will be working on starbase. First, find all the intersected miRNAs between starbase and encRNA_sensitivity_subset2
  starbase = process_starBase()
  intersected_miRNAs = intersect(unique(starbase$miRNA), unique(encRNA_sensitivity_subset2$miRNA))
  # subset starbase to include only miRNA shown up in encRNA_sensitivity_subset2;
  # similarly, subset encRNA_sensitivity_subset2
  starbase_subset = starbase[which(starbase$miRNA %in% intersected_miRNAs),]
  encRNA_sensitivity_subset3 = encRNA_sensitivity_subset2[which(encRNA_sensitivity_subset2$miRNA %in% intersected_miRNAs),]
  # now, find all intersected miRNA_mRNA pairs between encRNA_sensitivity_subset3 and starbase_subset
  intersected_lncRNA_miRNA_pairs = intersect(unique(encRNA_sensitivity_subset3$mRNA_miRNA_pair), unique(starbase_subset$mRNA_miRNA_pair))
  encRNA_sensitivity_subset4 = encRNA_sensitivity_subset2[which(encRNA_sensitivity_subset3$mRNA_miRNA_pair %in% intersected_lncRNA_miRNA_pairs),]
  return(encRNA_sensitivity_subset4)
}


process_starBase = function(){
  load("starbase_mRNA_miRNA_interactions.rda")
  processed = starbase_mrna_mirna_interaction
  colnames(processed)[1:2] = c("miRNA", "putative_mRNA")
  #dim(processed); View(processed)
  processed$miRNA = tolower(processed$miRNA)
  # foreach element, remove 3p and 5p parts
  processed$miRNA = gsub("-3p","",processed$miRNA)
  processed$miRNA = gsub("-5p","",processed$miRNA)
  processed$mRNA_miRNA_pair = paste(processed$putative_mRNA, processed$miRNA, sep = "-")
  processed = processed[,c("miRNA", "putative_mRNA", "mRNA_miRNA_pair")]
  processed = unique(processed)
  return(processed)
}



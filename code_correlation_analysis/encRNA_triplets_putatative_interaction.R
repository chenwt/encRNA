setwd("/media/ducdo/UUI/Bioinformatics/Summer Research/Cancer_Survival/encRNA_methylation_260616")
load("data_Saved_R_Objects/corr_matrices/normal_tumor_encRNA_triplets.rda")
load("mircode_objects.rda") 
load("starbase_mRNA_miRNA_interactions.rda") 


### -------------- some stats from sensitivity matrix ----------------------------------

dim(normal_encRNA) # [1] 29169    12
unique_normal_lncRNAs = unique(normal_encRNA$lncRNA); length(unique_normal_lncRNAs) # 200
unique_normal_mRNAs = unique(normal_encRNA$mRNA); length(unique_normal_mRNAs) # 1456
unique_normal_miRNAs = unique(normal_encRNA$miRNA); length(unique_normal_miRNAs) # 33

length(unique(normal_encRNA$lncRNA)) # [1] 200

## --------------- identify triple ----------------------------------------------------

## idea: starting from overlapped lncRNAs between sensitivity matrix and miRcode, then 
# get all possible lncRNA-miRNA-mRNA from miRcode and starbase, put them into a putative matrix. 
# then compare that matrix with the triplets in the sensitivity matrix

normal_encRNA_sensitivity_matched = get_matched_enRNA_sensitivity_with_putative_binding(normal_encRNA)
dim(normal_encRNA_sensitivity_matched)
tumor_encRNA_sensitivity_matched = get_matched_enRNA_sensitivity_with_putative_binding(tumor_encRNA)
dim(tumor_encRNA_sensitivity_matched)

intersect(normal_encRNA_sensitivity_matched$mRNA, tumor_encRNA_sensitivity_matched$mRNA) # 0
intersect(normal_encRNA_sensitivity_matched$miRNA, tumor_encRNA_sensitivity_matched$miRNA) # [1] "hsa-mir-22"
intersect(normal_encRNA_sensitivity_matched$lncRNA, tumor_encRNA_sensitivity_matched$lncRNA)
# [1] "ENSG00000229645.4" "ENSG00000228639.2" "ENSG00000249042.1"
# [4] "ENSG00000245812.2"

## --------------- helper function -----------------------------------------------------


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



get_putative_lncRNA_miRNA = function(dataframe){
  require(rlist)
  l = list()
  apply(dataframe, 1, function(r){
    k = getMiRNAs(r[2])
    l <<- list.append(l, k)
  })
  names(l) = dataframe$gene_id
  df = reshape2::melt(l)
  colnames(df) = c("miRNA", "putative_lncRNAs")
  df$lncRNA_miRNA_pair = paste(df$putative_lncRNAs, df$miRNA, sep = "-")
  return(df)
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

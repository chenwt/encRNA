changeSampleNameVersion = function(sampleNames = NULL, original = ".", type = "normal"){
  if (original == "."){
    names = strsplit(sampleNames, split = "[.]")
    names = unlist(lapply(names, function(name) paste(name[3],name[4],name[5], sep="-")))
    return(names)
  }else if (type == "normal"){
    names = strsplit(sampleNames, split = "[-]")
    names = unlist(lapply(names, function(name) paste("BRCA","Normal",name[1],name[2],name[3],sep=".")))
    return(names)
  }
  else{
    names = strsplit(sampleNames, split = "[-]")
    names = unlist(lapply(names, function(name) paste("BRCA","Tumor",name[1],name[2],name[3],sep=".")))
    return(names)
  }
}
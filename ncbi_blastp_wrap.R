### set some env variable
library(Biostrings)
# Path to fasta that corresponding to each row of information csv file.
# in our case "transcript_type_info.csv"
path_to_domain_files = "data/"

####=Functions required to run STEP 1 and 2 #####
# splitting function
splitseqpair <- function(path_dfiles, info_file, query_id){
  file_raw = open_input_files(paste0(path_dfiles, info_file$name, ".fasta"))
  i <- 0
  while (TRUE) {
    i <- i + 1
    aa <- readAAStringSet(file_raw, nrec=1)
    if (length(aa) == 0L)
      break
    cat("\n[INFO] Splitting", i, "...\n")
    if (grepl(pattern = query_id, aa@ranges@NAMES, fixed = TRUE) == TRUE){
      writeXStringSet(aa, 
                      paste0(info_file$name,"/query_", aa@ranges@NAMES, ".fasta"), 
                      format="fasta")
    }
    else{
      writeXStringSet(aa, 
                      paste0(info_file$name,"/subject_", aa@ranges@NAMES, ".fasta"), 
                      format="fasta")
    }
  }
}
# saving sub meta file
save_metafile <- function(path_dfiles, info_file, query_id){
  file_raw = open_input_files(paste0(path_dfiles, infofile$name, ".fasta"))
  i <- 0
  while (TRUE) {
    i <- i + 1
    aa <- readAAStringSet(file_raw, nrec=1)
    if (length(aa) == 0L)
      break
    cat("\n[INFO] saving information for", i, "...\n")
    if (grepl(pattern = query_id, aa@ranges@NAMES, fixed = TRUE) == TRUE){
      sub_meta_df[i, 1] = aa@ranges@NAMES
      sub_meta_df[i, 2] = "query"
    }
    else{
      sub_meta_df[i, 1] = aa@ranges@NAMES
      sub_meta_df[i, 2] = "subject"
    }
  }
}

#pick file for running blastp
pick_file <- function(foldername, sub_meta_df){
  queryseqdf <- sub_meta_df[sub_meta_df$seq_type=="query",]
  subjectseqdf <- sub_meta_df[sub_meta_df$seq_type=="subject",]
  for (q in queryseqdf$seq){
    qlist = unlist(strsplit(q, "_"))[-length(unlist(strsplit(q, "_")))]
    for (s in subjectseqdf$seq){
      slist = unlist(strsplit(s, "_"))[-length(unlist(strsplit(s, "_")))]
      if (identical(qlist, slist)==TRUE){
        runblast(foldername,q,s)
      }
    } 
  }
}
#runblastp function generic function
runblast <- function(foldername, query, subject){
  qfilename<- paste0(foldername,"/query_",query, ".fasta")
  sfilename<- paste0(foldername,"/subject_",subject, ".fasta")
  cat("\n[BLASTp] Processing", qfilename, "and\n", sfilename," for blast.")
  run_cmd <- "blastp"
  cmd_arg <- paste("-query",qfilename,"-subject",sfilename, "-outfmt 0")
  rr <- base::system2(command = run_cmd, 
                      args = cmd_arg,
                      stdout=TRUE,
                      stderr=TRUE,
                      wait = TRUE)
  saveoutput(foldername=foldername,result = rr, qy=query, sb=subject)
}
# save output in the alingment folder
saveoutput <- function(foldername, result,qy,sb){
  if (!dir.exists(paste0(foldername, "/alingment"))){
    dir.create(paste0(foldername, "/alingment"))
  }
  rr <- result
  out_filename = paste0(foldername, "/alingment/Alingment_",qy,"_",sb,"_out.txt")
  cat("\n[OUTFILE] saving result for",out_filename)
  write(rr, file = out_filename)
}

##### Step1: get the information from the file about pairs for domain name and ids ######
infofile <- read.csv("transcript_type_info.csv")
infofile <- infofile[!apply(is.na(infofile) | infofile == "", 1, all),]

####=====step 2 : RUN splitting===####
# split the file and save sub meta file for each domain from infofile
for(j in 1:nrow(infofile)){
  query_id = unlist(strsplit(infofile$prin[j], "_"))[1]
  subject_id = unlist(strsplit(infofile$alt[j], "_"))[1]
  if (!dir.exists(infofile$name[j])){
    dir.create(infofile$name[j])
    if(file.exists(paste0(path_to_domain_files, infofile$name[j], ".fasta"))){
      splitseqpair(path_dfiles=path_to_domain_files ,info_file = infofile[j], query_id = query_id)
      sub_meta_df <- setNames(data.frame(matrix(ncol = 2)), c("seq", "seq_type"))
      save_metafile(path_dfiles= path_to_domain_files, info_file = infofile[j], query_id = query_id)
    }
    else{
      cat("\n[Warning] Domain file:", infofile$name[j],".fasta Not found")
    }
    
  }
}
#### ---Step 3: RUN BLAST----#####
for (f in infofile$name){
  cat("\n [INFO]: Processing domain", f, "....\n")
  pick_file(foldername = f, sub_meta_df = sub_meta_df)
}

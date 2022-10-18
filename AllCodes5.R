

#' @title 
#' A helper function to runDeploid()that creates panel and plaf dataframes 
#' 
#' @description 
#' This function takes a long dataframe (containing differences in coding compared to the rCRS, 
#' in a select set of individuals) and returns panel and plaf dataframes 

#' @param long Table with vectors Sample.name (sample ID), Pos (nucleotide postion) and Nuc (nucleotide)

#' 
#' @examples 
#' createPlafPan(long)

createPlafPan<-function(EmpopLong) {
  Wide<-tidyr::spread(EmpopLong, Sample.name, Nuc)#go from long format to wide again
  Wide$Pos<-as.numeric(Wide$Pos)#coerce position to number to be able to sort
  Wide<-dplyr::arrange(Wide, Pos)#sort by position
  
  # final panel formating
  Wide[,2:ncol(Wide)]<-ifelse(is.na(Wide[ ,-1])==TRUE,0,1)#replace NA and AGTC with 0 and  1 respectively 
  panelFile<-Wide%>%dplyr::mutate(CHROM = "CHRM")#rename object, add mandatory CHRM column
  panelFile<-panelFile[,c(ncol(panelFile),1:ncol(panelFile)-1)]#and reorder
  
  
  EmpopLong %>%dplyr::group_by(Pos)%>%
    dplyr::summarise(PLAC=length(Nuc), #to get po level allele count or numerator
                     u=length(unique(Nuc))) %>% # pick out the number of unique alleles 
    dplyr::filter(u<=2) ->Numerator # filter out multi-allelic positions
  Numerator$Pos<-as.numeric(Numerator$Pos)# convert to numeric to allow sorting
  plafFile<-dplyr::arrange(Numerator, Pos) # sort by position
  plafFile<-dplyr::mutate(plafFile,PLAF = ifelse(PLAC==0 | PLAC==1, 
                                                       (PLAC + 1)/ (length(unique(EmpopLong$Sample.name))+2), # correction formula
                                                       PLAC/length(unique(EmpopLong$Sample.name))))
  
  #final PLAF formating 
  plafFile[-c(2,3)] ->plafFile # remove unnecesary columns
  plafFile<-plafFile %>%dplyr::mutate(CHROM = "CHRM")#rename object, add mandatory CHRM column
  plafFile<-plafFile[,c(ncol(plafFile),1:ncol(plafFile)-1)] #and reorder
  return(list(plafFile, panelFile))# return both panel and plaf dataframes as a list
  }


#' @title 
#' A helper function to runDeploid()that harmonizes the panel  
#' 
#' @description 
#' This function takes a panel (createPlafPan object),
#' and a vector q of sites from am AltRef (modified createAltRef object)
#' and matches or "harmonizes" the number of sites between both dataframes

#' @param panel a dataframe representing a reference panel for a set of individuals (createPlafPan object)
#' @param q vector of sites from an AltRef dataframe with reference and alternate allele counts (modified createAltRef object)
#' 
#' @examples 
#' harmonize(plaf,q)

harmonizePa<-function(panel,q) {
  
  pan<-dplyr::full_join(panel,q)%>%dplyr::arrange(Pos)# join panel and POS from alt ref keeping all rows
  pan<-dplyr::mutate(pan, CHROM=ifelse(is.na(CHROM), "CHRM", "CHRM"))# convert all NAs in the CHROM column to CHRM
  pan[is.na(pan)]<-0# convert all remaing NAs in the table to 0s
  CVPa<-pan # save final panel dataframe to enviroment
  
  return(CVPa) # return results of function
}


#' @title 
#' A helper function to runDeploid()that harmonizes the plaf  
#' 
#' @description 
#' This function takes a plaf (createPlafPan object),
#' and a vector q of sites from am AltRef (modified createAltRef object)
#' and matches or "harmonizes" the number of sites between both dataframes

#' @param plaf a dataframe of population level allele frequencies (createPlanPlaf object) 
#' @param q vector of sites from an AltRef dataframe with reference and alternate allele counts (modified createAltRef object)
#' 
#' @examples 
#' harmonize(plaf,q)

harmonizePl<-function(Plaf,q) {
  plaf<-dplyr::full_join(Plaf,q)%>% dplyr::arrange(Pos)# join plaf and POS from alt ref keeping all rows
  
  plaf<-dplyr::mutate(plaf, CHROM=ifelse(is.na(CHROM), "CHRM", "CHRM"))# convert all NAs in the CHROM column to CHRM
  plaf[is.na(plaf)]<-1e-6# convert all remaing NAs in the table to 1s
  CVPl<-plaf # save final plaf table to enviroment
  
  return(CVPl)# return results of function
}


#' @title 
#' A helper function to runDeploid()that harmonizes the AltRef  
#' 
#' @description 
#' This function takes an AltRef (modified createAltRef object),
#' and a vector p of sites from a plaf/pan (createPlafPan object)
#' and matches or "harmonizes" the number of sites between both dataframes

#' @param AltRef a dataframe of reference and alternate allele counts (modified createAltRef object) 
#' @param p vector of sites from a plaf/pan dataframe (createPlafPan object)
#' 
#' @examples 
#' harmonize(AltRef,p)

harmonizeAR<-function(AltRef, p){
  
  readCount <- AltRef$AltCount[1] + AltRef$RefCount[1]# create object representing total count
  
  AR<-dplyr::full_join(AltRef, p, by="Pos")%>% dplyr::arrange(Pos)# join alt ref tabel and POS from panel keeping all rows
  AR<-dplyr::mutate(AR, CHROM= ifelse(is.na(CHROM), "CHRM", "CHRM"))# convert all NAs in the CHROM column to CHRM
  AR$AltCount%>%tidyr::replace_na(0)->AR$AltCount # convert all NAs in ALT column to 0s
  AR$RefCount%>%tidyr::replace_na(readCount)->AR$RefCount # convert all NAs in REF column to total counts
  
  return(AR) # return results of function
}

#' @title 
#' A helper function to getHapDis()that calculates the hamming distance between two haplotypes 
#' 
#' @description 
#' This function takes numeric sample ID from a single phased haplotype (dEploid object) and a single true haplotype provided by the user,
#' and returns the hamming distance between the corresponding haplotype sequences

#' @param dEpHap a single sample ID associated with a phased haplotype (type numeric)
#' @param TruHap a single sample ID associated with a true haplotype (type numeric)
#' 
#' @examples 
#' HamDis(1,2)

HamDis<-function(dEpHap, TruHap, Dhaps, Thaps){ 
  dPos=dEpHap # create tibble of dEploid haplotype sampleIDs and allelic positions 
  tPos=TruHap # create tibble of true haplotype sampleIDs and allelic positions 
  DPos<-dplyr::filter(Dhaps, Samples==dPos)%>%dplyr::select(Pos) # select only those dEploid haplotypes/positions(rows) that have the sampleID specified in the arguments  
  TPos<-dplyr::filter(Thaps, Samples==tPos)%>%dplyr::select(Pos) # select only those true haplotypes/positions(rows) that have the sampleID specified in the arguments  
  haplos<-dplyr::bind_rows(DPos, TPos)%>%dplyr::add_count(Pos) # bind the selected haplotypes together by row
  haplos1<-dplyr::mutate(haplos, Dist=if_else(n<2, 1, 0)) # add a column to the combined haplotype dataframe 
  #where repeated positions are represented by a 0 and unique postions are represented by 0
  HamDis<-sum(haplos1$Dist) # sum the values of the new column to give hamming distance between dEploid and true haplotypes
  return(HamDis) # return the hamming distance value calculate in previous line
}


#' @title 
#' Get variants from SNP data
#' 
#' @description 
#' This function takes three vectors Position, Allele and Type from the SNPs 
#' and returns a vector of variants (e.g. "73G", "73+G", 73-" )

#' @param Pos A vector of genomic positions for the SNPs (type numeric)
#' @param Allele A vector of nucleotide bases present in the SNPs (character strings). Deletions are represented by ""
#' @param Type A vector of the type of mutation represented by the SNP. Possible values include "Substitution", "Insertion" 
#' and "Deletion" (charachter strings)
#' 
#' @examples 
#' Snp2variant("73", "G", "Substitution")
#' Snp2variant(c(73, 95, 100, 146), c("G", "", "c", "C"), c("Substitution", "Deletion", "Deletion", "Insertion"))

snp2variant<-function(Pos, Allele, Type){
  tib<-tibble::tibble(Pos=Pos, Allele=Allele, Type=Type) # make a tibble with three variables Pos, Allele and Type
  #and assign them the same values as are in the input vectors
  
  if (!(is.numeric(tib$Pos))){ # if the column Pos does not contain numeric data
    stop("Pos vector needs to be of the type 'numeric'") # give an error message to user
  } 
  v1<-if_else(tib$Allele=="A"|tib$Allele=="T"|tib$Allele=="G"|tib$Allele=="C"|tib$Allele=="a"|tib$Allele=="t"|tib$Allele=="g"|tib$Allele=="c"|tib$Allele=="", 0, 1) # if the column Allele contains ACTorG record a 0, otherwise a 1
  if (sum(tibble(v1))!=0) { # if the sum of the above vector is not 0 
    stop("Allele vector needs to be either A or T or G or C (case sensitive) or an empty string") # give an error message to user
  } 
  if (!(is.character(tib$Type))){ # if the column Type does not contain numeric data
    stop("Type vector needs to be of the type 'charachter'") # give an error message to user
  }
  v2<-if_else(tib$Type=="Substitution"|tib$Type=="Insertion" |tib$Type=="Deletion", 0, 1) # if the column Type contains Insertion/Deletion/Substiution record a 0, otherwise a 1
  if (sum(tibble(v2)!=0)) { # if the sum of the above vector is not 0 
    stop("Type vector needs to be either 'Substitution', 'Insertion' or 'Deletion' (case senstitive)") # give an error message to user
  } 
  
  
  Variant<-dplyr::case_when( # create and keep only one vector - Variant
    Type == "Substitution" ~ paste(Pos, Allele, sep=""), # if the Type is Substitution, concatenate Pos and Allele
    Type == "Insertion" ~ paste(Pos, Allele, sep = "+"), # if the Type is Insertion, concatenate Pos and Allele
    #seperated by +
    Type == "Deletion" & Allele == "" ~ paste(Pos, "-", sep = ""), # if the Type is Deletion, concatenate Pos and -
    Type == "Deletion" & (Allele == "a" | Allele =="t" | Allele =="g" | Allele =="c") ~ paste(Pos, Allele, sep = ""),
    TRUE               ~ "?"
  )
  
  return(Variant)
}


#' @title 
#' Write inforatuon on sample ID and variants in the EMPOP format 
#' 
#' @description 
#' This function takes two vectors SampleID and Variant
#' and returns a tibble in the EMPOP format.
#' Each row in the tibble will correspond to 1 individual (SampleID)
#' and their empop string (tab separated) will be the value in the second column
#' 
#' @param SampleID A vector of sample IDs (character strings)
#' @param Variant A vector of variants (character strings). See Snp2variant for description
#' 
#' @examples 
#' write_mbop("NA12871", c("73G", "95-"))
#' write_mbop(c("NA12871", "NA12871", "NA12872"), c("73G", "95-", "73G"))

write_mbop<-function(SampleID, Variant){
  tib<-tibble::tibble(SampleID=SampleID, Variant=Variant) #create tibble with two variables SampleID and Variant
  
  if (!(is.character(tib$SampleID))){ # if the column SampleID does not contain string data
    stop("SampleID vector needs to be of the type 'charachter'") # give an error message to user
  }
  
  if (!(is.character(tib$Variant))){ # if the column Variant does not contain string data
    stop("Variant vector needs to be of the type 'charachter'") # give an error message to user
  }
  
  # if (is.unsorted(tib$Variant, strictly = FALSE)){
  #   stop("Variant vector needs to be in ascending order")
  # }
  
  Empop<-tib%>%dplyr::group_by(SampleID)%>%dplyr::summarise(empopstring = paste(unique(Variant),collapse = "\t")) 
}


#' @title 
#' Get variants from an EMPOP file
#' 
#' @description 
#' This function takes an empop file containing at least four (mandatory) columns
#' (refer to https://empop.online/downloads for the emp file format) 
#' and returns a tibble with two columns - SampleID and Variant 
#' 
#' @param EMPOP An EMPOP file (tab seperated)
#' @param s a numeric argument that tells the function how many rows of data to skip while reading in EMPOP file
#' If no value is provided, the first line is ommited by default
#' 
#' @examples
#' Empop2variant("EMPOP.emp")
#' Empop2variant("EMPOP.emp", s = 3)

empop2variant<-function(EMPOP, s = 1){
  EMPOP<- readr::read_delim(EMPOP, delim = "\t", skip = s, col_names = FALSE) %>% # read the empop files in 
    dplyr::rename(SampleID = X1, Haplogroup = X2, Frequencies = X3) #rename the columns 
  LongF<- tidyr::gather(EMPOP,"Key","Variant",-SampleID, -Haplogroup, -Frequencies, factor_key = FALSE) #gather data 
  #to go from wided to long format 
  if(ncol(LongF)==3) { # if their is no info on variants in the EMPOP (i.e. sama as rCRS)... 
    dplyr::mutate(LongF,Variant = NA) %>% select(-Haplogroup, -Frequencies)-> LongF # add NA in the variant column
    #then remove all columns except Variant
  } 
  LongF%>%dplyr::select(SampleID, Variant) -> LongF #%>% drop_na() -> LongF1# Otherwise proceed with columns SampleID and Variant
  #longF<-dplyr::arrange(longF, SampleID )#sort by individuals
  
  return(LongF)
}


#' @title 
#' Get SNPs from variant data
#' 
#' @description 
#' This function takes a vector of variants (e.g. "73G", "73+G", 73-" ) and returns a tibble
#' with three columns : Pos, Allele and Type that contain information on the genomic position,
#' nucleotide base and type of mutation (i.e. "Substitution", "Insertion" and "Deletion")resecptively for each variant
#'  
#' @param Variant A vector of variants (character strings). See Snp2variant for description
#' 
#' @examples 
#' Variant2snp("73G")
#' Variant2snp(c("73G", "95-", "100c", "146+C"))

variant2snp<-function(Variant){
  
  if (!(is.character(Variant))){
    stop("Variant vector needs to be of the type 'charachter'")
  }
  
  # if (is.unsorted(Variant, strictly = FALSE)){
  #  stop("Variant vector needs to be in ascending order")
  # } 
  
  df<-(stringr::str_replace(Variant,"\\+", "\\."))# replace all occurences of "+" with a "."
  df1<-stringr::str_replace_all(df, "^[^0-9]","")%>% enframe()# replace any non numeric element at the start of the string with an empty string and make the result into a dataframe
  df2<-dplyr::rename(df1,Var = "value")%>% dplyr::select(-(1))# rename column to "Var" to be able to manipulate it downstream 
  df3<-df2%>%
    tidyr::extract(Var, into=c("neg", "Pos", "Ins", "Allele"), "^([A-Z-])?(\\d+)(\\.\\d*)?([^.]*)") %>% # split Var into four columns
    dplyr::mutate(Len=ifelse(Ins==".", 1, as.integer(sub(".", "", Ins))), # create a column "Len"
                  Allele = dplyr::case_when(                                     # populate the columen "Allele"...
                    grepl("^del|-", Allele, ignore.case = T)~ "",         # with an empty string if vector contains "del" or "-"
                    is.na(Len) | Len==1 ~ Allele,                         # with existing string if the Len is 1 or NA
                    !is.na(Len) | Len>1 ~ stringr::str_dup(Allele, Len),           # with existing string mutliplied by the number in column Len
                    TRUE ~ "?"                                            # with a "?" if none of these cases occur 
                  ),
                  Type = dplyr::case_when(                                                         # create a column "Type" and populate it...
                    grepl("\\.\\d|\\.", Ins)~ "Insertion",                                  # with "Insertion" if the column "Ins" contains a "." or digit
                    grepl("^[ACGTRYSWKMBDHVN]$", Allele, ignore.case = FALSE)~"Substitution",# with "Substitution" if the column "Allele" contains an upper-case nucleotide or the IUPAC code
                    grepl("^[acgtryswkmbdhvn]$", Allele, ignore.case = FALSE)~"Deletion",    # with "Deletion" if the column "Allele" contains a lower-case nucleotide or the IUPAC code 
                    Allele == "" ~ "Deletion",                                              # with "Deletion" if the column "Allele" has an empty string
                    TRUE ~ "?"                                                              # with a "?" if none of these cases occur 
                  )
    ) %>% 
    dplyr::select(-(neg), -(Ins), -(Len)) -> snpTab                               # drop all unnecessary columns 
  snpTab$Pos<-as.integer(snpTab$Pos)
  return(snpTab)  # return table
}


#' @title 
#' Join mutliple EMPOP files into one large EMPOP file using a look-up table 
#' 
#' @description 
#' This function takes a lookup-table file and a path to EMPOP files
#' and returns a bigger EMPOP file (i.e. all EMPOP files in the path bound together)
#' 
#'@param LUT A lookup-table file (tab seperated)
#'@param Pa The path to folder where EMPOP files are stored (charachter string)
#'
#'@examples
#'Hmtdb2Empop(LUT, Pa = R.home())

hmtdb2Empop<-function(LUT,Pa) {
  
  #function to add the amp column to the long format  
  add_amp<-function(e) {
    tib<-Empop2variant(e) # run Empop2Variant 
    tib$Amp<-stringr::str_extract(e, "[\\d]+") # add a column called "Amp" 
    #and populate with only the digits extracted from filenames
    return(tib) #return the dataframe
  }
  
  l<-list.files(Pa, pattern = ".emp")# read in the name of empop files in current folder
  empop<-lapply(l, add_amp)%>% #run through the list and for each file in the list add an "Amp" column 
    dplyr::bind_rows()%>%tidyr::drop_na()#  join all elements of the list into a dataframe and drop missing values
  
  
  t<-readr::read_tsv(LUT, col_names = FALSE )# read lookup table into R
  t%>%dplyr::select_if(~!(all(is.na(.)) | all(. == "")))-> t # remove empty columns
  t<-dplyr::rename(t, "Amp" = X5)#change the column name to "Amp" to be able to join with empop data later
  t$Amp<-as.character(t$Amp)# coerce number to charachter to be able to join with empop data later
  
  #empop1<- left_join(empop, t, by = "Amp" ) # join long dataframe with lookuptable 
  empop1<- dplyr::inner_join(t, empop, by = "Amp" ) # join long dataframe with lookuptable
  dplyr::anti_join(t, empop1, by = "X6")%>%dplyr::pull(X6) -> ids # pull out all individuals that have all amps same as rCRS
  empop2<-dplyr::mutate(empop1, Sort =(stringr::str_extract(Variant, "\\d+"))) # create a column "Sort" and populate it with just the numbers 
  #from the Variant column 
  empop2$Sort<-as.numeric(empop2$Sort)#convert "Sort" into type numeric
  empop3<-empop2%>%dplyr::arrange(Sort) # sort "Sort" in ascending order
  empop4<-dplyr::select(empop3, -(X1:X4), -(Amp), -(SampleID), -(Sort)) # delete unwanted columns
  #empop5<-rename(empop4, SampleID = ncol(empop4))# rename the last column to "SampleID"
  empop5<-dplyr::rename(empop4, SampleID = X6)# rename the last column to "SampleID"
  Empop<-empop5%>%dplyr::group_by(SampleID)%>%dplyr::summarise(empopstring = paste(unique(Variant),collapse = "\t"))#go from long to wide format
  Empop<-dplyr::select(Empop, SampleID, everything())%>% # move "SampleID" to the beginging of dataframe
    dplyr::mutate(Haplogroup = "", Frequencies = "") %>% # add the two other needed columns to dataframe
    dplyr::select(SampleID, Haplogroup, Frequencies, everything())# change order of the columns
  Empop <-tibble::add_row(Empop, .before = 1)%>% # add the first row for title of study and author details
    tibble::add_row(.before = 2) %>% # add the second row for geo background
    tibble::add_row(.before = 3, SampleID = "!#" ) # add third row for sequence range
  
  #create a tibble of indvidials with same amps as rCRS
  tibble::tibble(SampleID = unique(ids),
         Haplogroup = "",
         Frequencies = "",
         empopstring = "")-> foo
  
  Empop <- dplyr::bind_rows(Empop, foo)# bind the above tibble to the final EMPOP data  
  
  write.table(Empop, "Empop.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE, na="") # write tab-delimited file to file
  #without column or row names
}


#' @title 
#' Wrapper function for dEploid
#' 
#' @description 
#' This function takes information on alternate and reference allele counts,   
#' number of individuals and the recombination rate and runs 
#' the dEploid function on them. It retuns the dEploid output 
#' with includes estimated haplotypes, mixture proprotions and other phasing statistics. 
#' Additionaly this function also returns a list of nucleotide positions that violate the infinte sites model
#' and those that have been excluded by the user. 
#' 
#' @param l A list of sites excluded by the user
#' @param mPos A vector of genomic positions for the SNPs in the mixture (type numeric)
#' @param mAllele A vector of nucleotide bases present in the SNPs in the mixture (character strings)
#' @param Count A vector of read counts for each base (type numeric)
#' @param IsRef A vector containing either "Y" (Yes) or "N" (No), indicating whether allele is the same as reference or not
#' @param SampleID A vector of sample IDs in the reference data (character strings)
#' @param rPos A vector of genomic positions for the SNPs in the reference data (type numeric)
#' @param rAllele A vector of nucleotide bases present in the SNPs in the reference data (character strings)
#' @param NumMCMC A numeric argument that tells the function the number of MCMC samples (default value 800)
#' @param exportPostProb Save the posterior probabilities of the final iteration of all individuals
#' @param recomb numeric argument that gives the function the constant recombination probability (default value of 0.0)
#' @param k A numeric argument that tells the function the number of individuals in the mixture (default value 2, maximum of 5)
#' 
#' @examples
#  rundEploid(l, mPos, mAllele, Count, IsRef, SampleID, rPos, rAllele, NumMCMC=3000, exportPostProb, recomb=0.0, k=5)

rundEploid<-function(l=l, mPos, mAllele, Count, IsRef, SampleID, rPos, rAllele, NumMCMC=800, exportPostProb=TRUE, recomb=0.0, k=2){
  
  AR<-tibble::tibble(Pos=mPos, Nuc=mAllele, Count=Count, IsRef = IsRef)
  long<-tibble::tibble(Sample.name = SampleID, Pos=rPos, Nuc=rAllele)
  
  if (!(is.numeric(AR$Pos))){ # if the column Pos does not contain numeric data
    stop("mPos vector needs to be of the type 'numeric'") # give user an error message
  } 
  v1<-(if_else (AR$Nuc=="A"|AR$Nuc=="T"|AR$Nuc=="G"|AR$Nuc=="C", 0, 1)) # if the column Nuc contains ACTorG record a 0, otherwise a 1
  if (sum(tibble(v1))!=0) { # if the sum of the above vector is not 0 
    stop("mAllele vector needs to be either A or T or G or C") # give user an error message
  } 
  if (!(is.numeric(AR$Count))) { # if the column Count does not contain numeric data
    stop("Count vector needs to be of the type 'numeric'") # give user an error message
  }
  if (!(is.character(AR$IsRef))){ # if the column IsRef does not contain string data
    stop("IsRef vector needs to be of the type 'charachter'") # give user an error message
  }
  v2<-as.vector(if_else (AR$IsRef=="Y"|AR$IsRef== "N", 0, 1)) # if the column IsRef contains string data that is Y/N record a 0, otherwise a 1
  if (sum(tibble(v2))!=0) { # if the sum of the above vector is not 0 
    stop("IsRef vector needs to be either be Y (Yes) or N (No)") # give user an error message
  } 
  if (!(is.character(long$Sample.name))){ # if the column Sample.name does not contain string data
    stop ("SamplID vector needs to be of the type 'character'") # give user an error message
  }
  if (!(is.numeric(long$Pos))){ # if the column Pos does not contain numeric data
    stop("mPos vector needs to be of the type 'numeric'") # give user an error message
  }
  v3<-as.vector(if_else (long$Nuc=="A"|long$Nuc=="T"|long$Nuc=="G"|long$Nuc=="C", 0, 1)) # if the column Nuc contains ACTorG record a 0, otherwise a 1
  if (sum(tibble(v3)!=0)) { # if the sum of the above vector is not 0 
    stop("rAllele vector needs to be either A or T or G or C") # give user an error message
  } 
  
  empopLong<-dplyr::ungroup(long) # ungroup long dataframe 
  #identifying positions with infinite sites violations
  AR1<-subset(AR, IsRef=="N") #remove reference alleles
  V<-dplyr::bind_rows(AR1, empopLong)#bind AR to long table by appending rows 
  V%>%dplyr::group_by(Pos)%>% # group by position/site
    dplyr::mutate(UniqueAltAll=length(unique(Nuc)))-> V # calculate the number of unique alleles for each site
  V%>%dplyr::filter(UniqueAltAll >=2)%>%dplyr::select(Pos)%>%unlist()%>%unname()-> ISV # filter out infinte site violators
  exSites<-c(l, ISV)%>%sort()%>%as.data.frame()# merge user excluded sites and ISV 
  names(exSites)<-"ExSites" #name the vector
  V%>%dplyr::filter(!(Pos%in%exSites$ExSites))-> V2 #select the remaining sites to make  long dataframe
  
  
  AR3<-AR[which(!(AR$Pos %in% exSites$ExSites)),] # remove excluded sites from altref
  AR4<-AR3%>%dplyr::mutate(RefCount =dplyr::if_else(IsRef=="Y", Count, as.double(0)), AltCount =dplyr::if_else(IsRef == "N", Count, as.double(0)))
  AltRef <- dplyr::select(AR4, -c(Nuc, Count, IsRef))%>% mutate(CHROM = "CHRM")%>%dplyr::select(CHROM, everything()) # remove unnecessary columns 
  AltRef%>%dplyr::group_by(CHROM, Pos)%>%dplyr::summarise(RefCount = sum(RefCount), AltCount = sum(AltCount))->AltRef
  #convert remaining sites into long format
  dplyr::select(V2, -c(Count, IsRef, UniqueAltAll))%>% dplyr::select(Sample.name, Pos, Nuc)%>%dplyr::filter(Sample.name!= is.na(Sample.name)) -> EmpopLong #parse out updated long table
  
  #Create panel and plaf dataframes
  Plafpanel<-createPlafPan(EmpopLong)# run function that creates panel and plaf tables without candidate haps
  Pan<-Plafpanel[[2]] # split Pan from Plafpanel
  Plaf<-Plafpanel[[1]] # split Plaf from Plafpanel
  
  #Create vectors of positions for the harmonize functions
  p<-tibble::tibble(Pos=as.numeric(Pan$Pos))# create a table to accomodate just the POS column from panel 
  q<-tibble::tibble(Pos=as.numeric(AltRef$Pos))# create a table to accomodate just t he POS column from AltRef   #
  
  #harmonize Panel, Plaf and AltRef so that they all have same number of sites
  pan<-harmonizePa(Pan, q) # harmonize panel  
  plaf<-harmonizePl(Plaf, q) # harmonize plaf 
  altRef1<-harmonizeAR(AltRef, p) # harmonize AltRef 
  
  #write dataframes to file
  readr::write_tsv(pan,"Pan.tsv")#export data.frame as panelfile.tsv
  readr::write_tsv(plaf, "Plaf.tsv") #export data.frame as panelfile.tsv
  readr::write_tsv(altRef1[,c("CHROM", "Pos", "AltCount")],"Alt.tsv" )#write position and count of alternate alleles to tsv
  readr::write_tsv(altRef1[,c("CHROM", "Pos", "RefCount")], "Ref.tsv")#write postion and count of reference alleles to tsv
  
  #assigning objects to dEploid input files
  Pan<-"Pan.tsv" #assinging the panel file
  Plaf<-"Plaf.tsv" #assigning the plaf file
  Alt<-"Alt.tsv" #assigning the alt file
  Ref<-"Ref.tsv" # assinging the ref file
  
  #core dEploid function
  Com<-(paste("-ref", Ref, "-alt", Alt, "-plaf", Plaf, "-panel", Pan, "-nSample", NumMCMC,
                    ifelse(exportPostProb==TRUE, "-exportPostProb", ""), "-recomb", recomb, "-k", k, sep = " "))
  dEploid.run<-DEploid::dEploid(Com)
  
  
  colnames(dEploid.run$Haps)<-altRef1$Pos # assign column names to the dEploid matrix to associate sites to alleles
  return(append(dEploid.run, unique(exSites))) #return modified dEploid object with a list of positions that violate infinite sites model
  } 


#' @title 
#' Function to get nearest neighbors of the mixed haplotypes to construct reference panel
#' 
#' @description 
#' This function takes information on alternate and reference allele counts of the mixed haplotypes,  
#' the rCRS, and sequences from mtDB to identify nearest neighbors of the mixed haplotyoes based on an edit distance selected by the user. 
#' It retuns a vector of haplotypes that is then used to make the reference panel.
#' 
#' 
#' @param mPos A vector of genomic positions for the SNPs in the mixture (type numeric)
#' @param mAllele A vector of nucleotide bases present in the SNPs in the mixture (character strings)
#' @param Count A vector of read counts for each base (type numeric)
#' @param IsRef A vector containing either "Y" (Yes) or "N" (No), indicating whether allele is the same as reference or not
#' @param ref A vector representing the rCRS (character strings)
#' @param mtGenomes A list of mitogenomic sequences from hmtDB (characther strings)
#' @param db Haplotype sequences curated from mtDB (S4 object of the class SQLiteConnection)
#' @param ED A numeric value that tells the function what edit distance to use to capture nearest neighbours (default value of 4)
#' @param l A list of sites excluded by the user (character strings)
#' 
#  getNeNe(mPos, mAllele, Count, IsRef, ref=ref, mtGenomes=mtGenomes, db=db, ED=4, l=l )


#Funtion to provide a list of nearest neighbors for any two individuals in database
getNeNe <- function(mPos, mAllele, Count, IsRef, ref, mtGenomes=mtGenomes, db=db, ED=4, l=l) {
  
  AR<-tibble::tibble(Pos=mPos, Nuc=mAllele, Count=Count, IsRef = IsRef)
  
  mix<-dplyr::rename(AR, "position" = "Pos", "basecall" = "Nuc")
  mixy<-dplyr::filter(mix,!(position%in%l)) #remove user selected sites
  mixy$position<-as.integer(mixy$position)
  
  # We need to emualate the basecalls as per converge
  # namely-- if the two people have the same allele, we need to have *1 row*; in the current form there's two rows with the same INFO
  # (other than who has the allele)
  
  # if the two people have different alleles and their both not the rCRS, we need two rows (with both alleles)
  # (deploid will filter that out anyways)
  # AND
  # if the two people are different and one == the rCRS at that position
  # the current representation has 1 row
  #
  # for converge processing there would be two alleles at this location instead
  # (a "heterozygote" call)
  # which would get split into two rows, one for each allele (with the value of each row either being the allele call or the reference)
  
  # get all of the reference alleles
  dplyr::pull(mixy, position) %>% unique() %>% sort -> mixPos
  refAlleles <- substring(ref, mixPos, mixPos)
  refTib <- tibble::tibble(position=mixPos, RefAllele=refAlleles)
  
  # add a RefAllele column; inner_join is equivalent here
  dplyr::left_join(mixy, refTib, by="position") ->mixy1
  
  dplyr::group_by(mixy1, position) %>%
    dplyr::summarize(N=dplyr::n(), # number of allele calls at site. 1 or 2
                     A1=basecall[[1]], # no matter what, first allele is what we have
                     A2=ifelse(N==1, RefAllele[[1]], basecall[[2]]), # second allele is the ref if there's no other info,
                     basecall=ifelse(A1==A2, A1, paste(A1, A2, sep=";"))  # concatenate the two records as a ; separated string
    ) -> mixSplit
  
  
  dplyr::select(mixSplit, position, basecall) %>%
    tidyr::separate_rows(basecall, sep=";") %>% # split the Col string; two rows for het calls; 1 row for hom calls
    dplyr::mutate(event="X") %>% # define the event to be a substitution
    dplyr::arrange(position) -> mixAsPerConverge
  
  
  # make a sequence graph
  # supports indels...
  makeDeploidSeqGraph(db, mixAsPerConverge$position, mixAsPerConverge$basecall, mixAsPerConverge$event, ignoreIndels=TRUE) -> gr
  
  # all edit distances; capped at 4 (0-3 are meaningful distances)
  mtGenomes$edDist <- Haplotypical::fastBoundedHammingGraphDist(gr, mtGenomes$sequence, ED)
  dplyr::filter(mtGenomes, edDist<4) %>%dplyr::pull(sampleid) %>% as.character() -> nearNeighbors
  # old code (levenshtein distance)
  # mtGenomes$edDist <- -1
  
return(nearNeighbors)
  }


#' @title 
#' Get proportion of mixtures from dEploid output
#' 
#' @description 
#' This function takes the output of dEploid and returns mixture proportions for major and minor contributors
#' 
#' 
#' @param dEploid.run dEploid object 
#'
#'@examples
#'getMixProps(dEploid.run)
#'

getMixProps<-function(dEploid.run){
  MixProps<-tail(dEploid.run$Proportions, n=1)%>%sort(decreasing = TRUE)
}


#' @title 
#' Get estimated haplotypes based on the dEploid output
#' 
#' @description 
#' This function takes the output of dEploid and an Altref table and returns estimated haplotypes of the 
#' contributors and their genomic positions
#' 
#' 
#' @param dEploid.run dEploid output 
#' @param mPos A vector of genomic positions for the SNPs in the mixture (type numeric)
#' @param mAllele A vector of nucleotide bases present in the SNPs in the mixture (character strings)
#' @param Count A vector of read counts for each base (type numeric)
#' @param IsRef A vector containing either "Y" (Yes) or "N" (No), indicating whether allele is the same as reference or not
#' 
#' @examples getdEploidHaps(dEploid.run, mPos, mAllele, Count, IsRef)

getdEploidHaps <-function(dEploid.run, mPos, mAllele, Count, IsRef) {
  
  dEploid.run<-dEploid.run
  AR<-tibble::tibble(Pos=mPos, Nuc=mAllele, Count=Count, IsRef = IsRef)
  AR <- dplyr::filter(AR, IsRef == "N")
  if (!(is.numeric(AR$Pos))){ # if the column Pos does not contain numeric data
    stop("mPos vector needs to be of the type 'numeric'")# give user an error message
  } 
  v1<-(if_else (AR$Nuc=="A"|AR$Nuc=="T"|AR$Nuc=="G"|AR$Nuc=="C", 0, 1)) # if the column Nuc contains ACTorG record a 0, otherwise a 1
  if (sum(tibble(v1))!=0) { # if the sum of the above vector is not 0 
    stop("mAllele vector needs to be either A or T or G or C") # give user an error message
  } 
  if (!(is.numeric(AR$Count))) { # if the column Count does not contain numeric data
    stop("Count vector needs to be of the type 'numeric'") # give user an error message
  }
  if (!(is.character(AR$IsRef))){ # if the column IsRef does not contain numeric data
    stop("IsRef vector needs to be of the type 'charachter'") # give user an error message
  }
  v2<-as.vector(if_else (AR$IsRef=="Y"|AR$IsRef== "N", 0, 1)) # if the column IsRef contains string data that is Y/N record a 0, otherwise a 1
  if (sum(tibble(v2))!=0) { # if the sum of the above vector is not 0 
    stop("IsRef vector needs to be either be Y (Yes) or N (No)") # give user an error message
  } 
  
  Haps<-dplyr::as_tibble(dEploid.run$Haps)%>% tibble::rownames_to_column("Samples") # convert Haps from matrix to tibble
  tidyr::gather(Haps,"Key", "Value", -Samples)-> Long # go from wide to long
  Long[order(Long$Samples), c(1,2,3)]->HapsTab # reorder columns
  HapsTab1<- dplyr::filter(HapsTab, Value!=0)%>% dplyr::rename("Pos"=Key) # remove all rows that have a reference allele (i.e. 0)
  HapsTab1$Pos<-as.numeric(HapsTab1$Pos) # convert the the vector POS to numeric to be able to join
  HapsTab2<-dplyr::left_join(HapsTab1, AR, by = "Pos")%>% dplyr::select(Samples, Pos, Nuc) # join dEploid haps with AltRef data to
  #to associate the site to a unique alternate allele and then remove unecessary columns
  }


#'@title
#'Get Hamming distances between the haplotypes estimated by dEploid and the true haplotypes
#' #'
#'@description
#'This function takes the haplotypes estimated by dEploid and the true haplotypes and returns the Hamming distance between them
#' #'
#'@param dSampleID A vector of sample IDs for the phased haplotypes provided by the dEploid object (type numeric)
#'@param dPos A vector of genomic positions for the phased haplotypes provided by the dEploid object (type numeric)
#'@param tSampleID A vector of sample IDs for the true haplotypes provided by the user (charachter strings or type numeric)
#'Note that the number of true haplotypes can not be larger than the number of phased haplotypes (type numeric)
#'@param tPos A vector of genomic positions for the true haplotypes provided by the user (type numeric)

#'@examples getHapDis(c(1,1,1,2,2,2), c(73, 95, 146, 95, 146, 750), c(1,1,2,2,2), c(95, 146, 95, 146, 750))

getHapDis<-function(dSampleID, dPos, tSampleID, tPos){
  
  Dhaps<-tibble::tibble(Samples=dSampleID, Pos=dPos) #make a tibble with vectors of dEploid sampleID and allelic positions
  Thaps<-tibble::tibble(Samples=tSampleID, Pos=tPos) #make a tibble with vectors of true sampleID and allelic positions
  
  if (!(is.numeric(Dhaps$Pos))){ # if the column Pos in  df Dhaps does not contain numeric data
    stop("dPos vector needs to be of the type 'numeric'") # give user an error message
  }
  
  if (!(is.numeric(Thaps$Samples))){ # if the column Samples in df Thaps does not contain numeric data
    stop("tSampleID vector needs to be of the type 'numeric'") # give user an error message
  }
  
  if (!(is.numeric(Thaps$Pos))){ # if the column Pos in df Thaps does not contain numeric data
    stop("tPos vector needs to be of the type 'numeric'") # give user an error message
  }
  
  n<-(nrow(subset(Dhaps,!duplicated(Samples)))) # get the number of individuals haplotypes provided by dEploid
  perms<-gtools::permutations(n, n, 1:n) # calculate all possible permutations of odering the haplotypes
  perms1<-dplyr::rowwise(as.data.frame(perms))%>% dplyr::group_split()#split the permutations by groups (group==1 specific permutation)
  names(perms1)<-1:length(perms1) # provide each group in the list an unique ID (to be able to group_by later) 
  perms2<-lapply(perms1, t) # transpose each group in the list
  perms3<-do.call(rbind, lapply(perms2, as.data.frame)) # bind all groups in list into one dataframe
  perms4<-tibble::rownames_to_column(perms3, var = "Groups") # use the rownames to make a column of unique IDs for each group
  perms4$Groups<-as.numeric(gsub("\\..*","", perms4$Groups)) # remove everything after the "." to make cleaner IDs
  V2<-unique(Thaps$Samples) # create a vector of true sampleIDs
  length(V2)<-n # modify the new vectors length to be the same as the lenght of the dEploid sampleID vector 
  perms4$TrueHap<-V2 # add the vector as a column to the main dataframe
  perms5<-dplyr::rename(perms4, PhasedHap=V1) # rename the column of dEploid sampleIDs
  perms6<-stats::na.omit(perms5) # remove rows with NA
  perms6$HD<-purrr::map2(perms6$PhasedHap, perms6$TrueHap, HamDis, Dhaps, Thaps) # run the hamming distance function on each row of the dataframa
  #and add the resulting values as a column to the main dataframe
  Sumhd<-perms6%>%dplyr::group_by(Groups)%>%dplyr::summarise(SumHD=sum(as.numeric(HD))) # sum the hamming distance values by group
  MinHD<-dplyr::filter(Sumhd, SumHD==min(SumHD)) # select the group with the lowest hamming distance
  g<-as.numeric(MinHD$Groups[1]) # create a variable to store the group ID with the smallest hamming distance
  return(dplyr::filter(perms6, Groups==g)%>%dplyr::select(-Groups)) # return the sampleIDs and hamming distances of the group with the smallest hamming distance
  }


#'@title 
#'Get the haplotype of the minor contributor
#' 
#'@description 
#'This function interpolates the minor haplotype if the major haplotype has been estimated accurately in a two person mixture. 
#' 
#'@param major A numeric identifier for the major haplotype
#'@param dEploid.run dEploid output 
#'@param mPos A vector of genomic positions for the SNPs in the mixture (type numeric)
#'@param mAllele A vector of nucleotide bases present in the SNPs in the mixture (character strings)
#'@param Count A vector of read counts for each base (type numeric)
#'@param IsRef A vector containing either "Y" (Yes) or "N" (No), indicating whether allele is the same as reference or not

#'@examples 
#'getMinor(major=1, dEploid.run, mPos, mAllele, mCount, IsRef)

getMinorHap <-function( major, dEploid.run, mPos, mAllele, mCount, IsRef) {
  
  dEploid.run<-dEploid.run
  AR<-tibble::tibble(Pos=mPos, Nuc=mAllele, Count=mCount, IsRef = IsRef)
  
  if (!(is.numeric(major))){ # if the value of the argument "major" is not numeric
    stop("major needs to be of the type 'numeric'")# give user an error message
  }
  if (!(is.numeric(AR$Pos))){ # if the column Pos does not contain numeric data
    stop("mPos vector needs to be of the type 'numeric'")# give user an error message
  }
  v1<-(if_else (AR$Nuc=="A"|AR$Nuc=="T"|AR$Nuc=="G"|AR$Nuc=="C", 0, 1)) # if the column Nuc contains ACTorG record a 0, otherwise a 1
  if (sum(tibble(v1))!=0) { # if the sum of the above vector is not 0
    stop("mAllele vector needs to be either A or T or G or C") # give user an error message
  }
  if (!(is.numeric(AR$Count))) { # if the column Count does not contain numeric data
    stop("Count vector needs to be of the type 'numeric'") # give user an error message
  }
  if (!(is.character(AR$IsRef))){ # if the column IsRef does not contain numeric data
    stop("IsRef vector needs to be of the type 'charachter'") # give user an error message
  }
  v2<-as.vector(if_else (AR$IsRef=="Y"|AR$IsRef== "N", 0, 1)) # if the column IsRef contains string data that is Y/N record a 0, otherwise a 1
  if (sum(tibble(v2))!=0) { # if the sum of the above vector is not 0
    stop("IsRef vector needs to be either be Y (Yes) or N (No)") # give user an error message
  }
  eMajor<-dplyr::filter(Haps,Samples==major) #extract the alleles of the correctly estimated sample
  eMajor$Samples<-as.integer(eMajor$Samples)
  AR1<-dplyr::left_join(AR,eMajor) # join with the altref file to figure out which alleles belong to the second individual
  AR2<-dplyr::filter(AR1,IsRef=="N") #remove all reference alleles
  Max<-max(AR2$Count)# store total read count
  AR3<-dplyr::filter(AR2, is.na(Samples)|Count==Max) # only extract sites that are either NA (only found in minor) or homozygous alternate (found in both)
  minor1<-dplyr::select(AR3, Pos, Nuc)%>%dplyr::filter(!Pos%in%dEploid.run$ExSites)# clean up dataframe and remove sites that violate the infinite sites model
  return(minor1) # return the haplotype of the minor
} 


#' Title
#' A function to decode IUPAC coded SNPs and split indels per position in long format.
#' 
#' @description 
#' This function takes three vectors of genomic positions, nucleotide base, and type of mutation respectively for each variant and decodes 
#' any IUPAC code to its respective alleles in a long format. This also converts any mixture indels into a long format.
#' 
#' @param Pos: A vector of genomic positions of alleles
#' @param Allele: A vector of nucleotide bases 
#' @param Type: A vector of type of mutation  
#'
#' @return: Returns a tibble of genomic positions, nucleotide bases, and type of mutation for all variants in a long format
#' @export
#'
#' @examples
#' \dontrun{
#' # UnfoldSNP(Pos = c(73L,152L,217L,264L,512L),
#' #           Allele=c("G","C","R","ac",""),
#' #           Type=c("Substitution","Substitution","Substitution","Insertion","Deletion"))
#' }
#' 
UnfoldSNP<-function(Pos, Allele, Type){
  if (!(is.integer(Pos))){
    stop("Pos vector needs to be of the type 'integer'")
  }
  
  if (!(is.character(Allele))){
    stop("Allele vector needs to be of the type 'character'")
  }
  
  if (!(is.character(Type))){
    stop("Type vector needs to be of the type 'character'")
  }
  
  
  mydata <- tibble(Pos, Allele, Type)
  IUPAC_amb_codes <- tribble(
    ~IUPACcode, ~DecodAllele,
    "M", "A", "m", "a",
    "M", "C", "m", "c",
    "R", "A", "r", "a",
    "R", "G", "r", "g",
    "W", "A", "w", "a",
    "W", "T", "w", "t",
    "S", "C", "s", "c",
    "S", "G", "s", "g",
    "Y", "C", "y", "c",
    "Y", "T", "y", "t",
    "K", "G", "k", "g",
    "K", "T", "k", "t",
    "V", "A", "v", "a",
    "V", "C", "v", "c",
    "V", "G", "v", "g",
    "H", "A", "h", "a",
    "H", "C", "h", "c",
    "H", "T", "h", "t",
    "D", "A", "d", "a",
    "D", "G", "d", "g",
    "D", "T", "d", "t",
    "B", "C", "b", "c",
    "B", "G", "b", "g",
    "B", "T", "b", "t",
    "N", "G", "n", "g",
    "N", "A", "n", "a",
    "N", "T", "n", "t",
    "N", "C", "n", "c",
    "A", "A", "a", "a",
    "T", "T", "t", "t",
    "G", "G", "g", "g",
    "C", "C", "c", "c"
  )
  
  mytib = mydata %>% 
    dplyr::left_join(IUPAC_amb_codes, by = c("Allele" = "IUPACcode" )) %>%
    dplyr::mutate(Allele = if_else(is.na(DecodAllele),
                                   Allele,
                                   DecodAllele),
                  DecodAllele=NULL) %>%
    dplyr::mutate(Allele=ifelse(Type=="Deletion" & Allele=="","",Allele),
                  Pos=as.integer(Pos))
  rbind(mytib, mytib %>%
          filter(grepl("[acgtryswkmbdhvn]+", Allele)) %>% 
          dplyr::mutate(Allele ="")) %>% 
    dplyr::mutate(Allele=str_to_upper(Allele)) %>%
    dplyr::arrange(Pos) -> myfinaldata  
  
  
  return(myfinaldata)
}

library(readxl)
#@param: unfoldVar - output of UnfoldSNP
#@param: quantXL - excel file from converge that have allele counts

# quantify unfolded variants Empop2AltRef
Empop2AltRef <- function(unfoldVar, quantXL, normVal = 100){
  
  mydata_xl <- read_excel(quantXL,skip = 7)
  mydata_xlLong <-  mydata_xl %>% 
    filter(Type == "SNP") %>%
    dplyr::select(Position:Variant,Type,Polymorphism,contains("%")) %>%
    tidyr::pivot_longer(-c(Position,Ref,Sample,Variant,Type,Polymorphism), names_to="Bases_indels", values_to="Percentage") %>%
    tidyr::separate(Bases_indels, into="Bases_indels",extra="drop") %>%
    dplyr::mutate(Position=as.integer(Position),
                  Percentage=as.numeric(Percentage)) %>%
    dplyr::mutate(Type = ifelse(Type == "SNP", "Substitution", Type))
  
  #browser()
  mydata_2gether <- unfoldVar %>% filter(Type == "Substitution") %>% 
    dplyr::left_join(mydata_xlLong, 
                     by=c("Pos"="Position",
                          "Allele"="Bases_indels",
                          "Type"="Type")) %>%
    group_by(Pos) %>%
    mutate(sumn = sum(Percentage))
  
  mydata_final <- mydata_2gether %>% 
    mutate(#browser(),
      Count = round((Percentage/sumn)* normVal)
    ) %>%
    mutate(isRef = ifelse(Allele == Ref, "Y","N")) %>% 
    dplyr::rename(NormalizedCount = "Count") %>%
    select(Pos, Allele, NormalizedCount, isRef)
  
  # trycatch to display the pos with 0 readcount. This is only valid if the
  # normalization is set to 100(default) and readcounts are between 0-0.5% in 
  # the excel file (quantitative file from converge)
  tryCatch(
    if(nrow(mydata_final %>% filter(NormalizedCount == 0)) >= 1){
      stop("There are alleles with zero read count.")
    },
    error = function(e){
      cat("There [is an]/ are alleles with read count '0'. Please recheck the empop",
          "file. Read count of such allele[s] will be converted to read count of '1'\n")
      mydata_final %>% filter(NormalizedCount == 0) %>% pull(Pos) -> allWithzero
      cat("The position[s] with read count of zero is/are: ", allWithzero, "\n")
    }
  )
  
  mydata_final <- mydata_final %>% 
    mutate(NormalizedCount = 
             if_else(NormalizedCount == 0, 
                     1, 
                     NormalizedCount))
  
  return(mydata_final)
  
}





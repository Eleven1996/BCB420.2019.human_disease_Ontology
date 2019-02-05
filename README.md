# `BCB420.2019.Human_disease_Ontology`

#### (**HDO** data annotation of human genes)

&nbsp;

###### Liwen Zhuang &lt;zhuangyedda@gmail.com&gt;

----

**If any of this information is ambiguous, inaccurate, outdated, or incomplete, that is because I wrote this with fever,please just slap me.**

----

<!-- TOCbelow -->
1. About this package
2. Human_disease_Ontology Data
3. Data Download and clearnup
4. Mapping HDO to HGNC
5. Annotation
6. Annotation of example gene set
7. Reference and citation
8. Acknowledgements
<!-- TOCabove -->

----


# 1 About this package:
This package describe the workflow to download disease ontology data from [**Human Disease Ontology Database**](http://www.obofoundry.org/ontology/doid.html),how to map the IDs to HGNC symbols when avaliable,how to analyze the database on different aspects.<br/>

Package checks pass without errors,warnings or notes.
structure:
```text
 --BCB420.2019.human_disease_Ontology/
   |__.gitignore
   |__.Rbuildignore
   |__BCB420.2019.human_disease_Ontology.Rproj
   |__DESCRIPTION
   |__dev/
      |__functionTemplate.R
      |__rptTwee.R
   |__inst/
      |__extdata/
         |__HDOsys.RData
      |__scripts/
         |__scriptTemplate.R
   |__LICENSE
   |__NAMESPACE
   |__R/
      |__toBrowser.R
      |__zzz.R
   |__README.md
```
&nbsp;

# 2 Human Disease Ontology Data
HDO(Human Disease Ontology) is a sub database of OLS(Ontology Lookup Service)[https://www.ebi.ac.uk/ols/index],which is a repository for biomedical ontologies that aims to provide a single point of access to the latest ontology versions. HDO has been developed as a standardized ontology for human disease with the purpose of providing the biomedical community with consistent, reusable and sustainable descriptions of human disease terms, phenotype characteristics and related medical vocabulary disease concepts.<br/>
The database has one [Term] for each disease, each [Term] has following infomration:<br/>

* a unique DOID to indicate the [Term]
* the parent of this [Term]: is this a sub-disease of some disease
* xref:cross reference to other database of this [Term]

**Example**<br/>
[Term]<br/>
id: DOID:0002116<br/>
name: pterygium<br/>
def: "A corneal disease that is characterized by a triangular tissue growth located_in cornea of the eye that is the result of collagen degeneration and fibrovascular proliferation." <br/>
synonym: "surfer's eye" EXACT []<br/>
xref: MESH:D011625<br/>
xref: UMLS_CUI:C0033999<br/>
is_a: DOID:10124 ! corneal disease<br/>
created_by: laronhughes<br/>
creation_date: 2010-06-30T02:44:30Z<br/>

&nbsp;

#3 Data Download and clearnup
To download the source data from Human Disease Ontoogy:

1. Navigate to the [**Human Disease Ontology Database**](http://www.obofoundry.org/ontology/doid.html)
2. choose "HumanDO.obo"(4.7Mb),copy the url
3. in R,use the following code , OR just download to your filepath
```R
file<-download.file("https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/master/src/ontology/HumanDO.obo",
                    destfile = "./inst/extdata/HumanDO.obo",method = "curl")
```
we also need data from OMIM database later, similarly:

1. Navigate to [**OMIM Database**](http://omim.org/downloads/)
2. choose "mim2gene.txt"(857.6KB),copy the url
3. download

#4 Mapping DOID to HGNC symbols
HDO are basicly disease names, some of them may have xref to OMIM ID, some of the OMIM ID have Ensemble IDs,and Ensemble IDs can used to map to HGNC symbols. HDO->OMIM->Ensemble->HGNC

#### Preparations: packages, functions, files

**`biomaRt`** (optional)biomaRt is a Bioconductor package that implements the RESTful API of biomart,
the annotation framwork for model organism genomes at the EBI. It is a Bioconductor package, and as such it needs to be loaded via the **`BiocManager`**,
&nbsp;

```R
if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}
```
We need to fetch HGNC reference data from Github.

1. Navigate to [https://github.com/hyginn/BCB420-2019-resources]
2. copy HGNC.RData
3. paste to your own working direcotry

&nbsp;

#### 4.1 Step one: how to read obo file?
[**Obo files**](http://www.geneontology.org/faq/what-obo-file-format) are "The OBO file format is one of the formats that the Gene Ontology is made available in." One of the very important concept of ontology is the ancestor/descendant relationship,there is a R package called "ontologyIndex" used to do this kind of analysis, but for the purpose of this project, we only want to map the ids in different database, so we want all information presented in a flat-form (dataframe) way. There isn't any avaliable packages for reading obo files into dataframe, but we can treat it as a simple txt file, then parse with regular expression to extract the inforamtion we need. I found this great and useful [**github resource**](https://github.com/Molecular-Proteomics/roboparse/blob/master/main.R) for doing this.<br/>
By reading the code, we can see that it use "startWith" function to extract information we required. So when we are using the function, the "col"" attribute need to be "xref: OMIM".(Instead of "xref", since the result will only show the first xref.We are not interested in other cross-references, so we should specificly ask it to find OMIM id). Also,note that the "obofile" attribute is a character vector from readLines, not filepath.

```R
#
# Edward Lau 2017
# lau1@stanford.edu
#

readOboFile <- function(obofile, cols) {
    "
    Read in an obo file and return a flat table. 
    obofile:    character vector from readLines() of an obo file
    cols:       list of values to read, e.g., c('name', 'id')
    "
  
    require("dplyr")

    # Make an empty dataframe to host all entries, with the column names matching the supplied column lists
    dataframe <- matrix(ncol=length(cols), nrow=0) %>% as.data.frame()
    colnames(dataframe) <- cols
    
    current_section = ""
    
    # Iterate each line
    for (i in 1:length(obofile)) {
        print(i)
      
        # Skip empty line
        if (obofile[i] == "") { 
            next
        
        # If this is a [Term] line, start a new mini dataframe "current_entry" and grab the entries specified by the cols arg
        } else if (obofile[i] == "[Term]") {
          
              current_section <- "Term"
              # Check if "current entry" already exists, if it does, then bind to master dataframe and remove it
              if (exists("current_entry")) {
                  #print(current_entry)
                  dataframe <- bind_rows(dataframe, current_entry)
                  rm(current_entry)
              }
                
              # If "current entry" doesn't already exist or is just removed, create the new mini dataframe for the next entry
                current_entry <- matrix(ncol=length(cols), nrow=1) %>% as.data.frame()
                colnames(current_entry) <- cols
        
        # Skipping over TypeDef section so we are only pulling things from [Term]
                # Note this only works if the TypeDef is at the
        } else if (obofile[i] == "[Typedef]"){
            if(exists("current_entry")){
               current_section <- "Typedef"
            }
              
        # For any other line, check if starts with the same as the specified column ID  
        } else {
            if (exists("current_entry") & current_section == "Term") {
                for (j in 1:length(cols)){
                    if (startsWith(obofile[i], cols[j])){
                        current_entry[1,j] = gsub(paste0("^",cols[j],": "),"", obofile[i])
                    }
                }
            }
        }
    }
    
    # Store the last entry since we are not looping back to 
    if (exists("current_entry")) {
        dataframe <- bind_rows(dataframe, current_entry)
        rm(current_entry)
    }
    
    return(dataframe)
}
```
define the file and col we want to retrieve
```R
HDOfile<- readLines("./inst/extdata/HumanDO.obo")
HDOcols<- c('name', 'id',"xref: OMIM","def")
HDOdf<-readOboFile(HDOfile,HDOcols) #11514 obs. of 3 variables
```
&nbsp;

#### 4.2 Step two: adjust the dataframe
After getting the result, we will see that a lot of the rows with missing OMIM ids, we will not need these since OMIM id is the only reference possible in the dataset to map to HGNC symbols. We need to remove these from our dataframe. Also,in the third column(the OMIM id column), each element are with a prefix "xref: OMIM:", this is due to the readOboFile function,also in the DOID column, we have a prefix "DOID:" with similar reason, we want to get a more concise and informative dataframe.

```R
#disease without OMIM id will not be able to map to HGNC symbol, thus we need to filter these rows
NEWdf<-HDOdf[!is.na(HDOdf$`xref: OMIM`), ] #3389 obs. of 3 variables
#remove prefix "xref: OMIM"
for (index in 1:nrow(NEWdf)){
  NEWdf[index,3]<- gsub("xref: OMIM:","",NEWdf[index,3])
}
#remove prefix "DOID"
for (index in 1:nrow(NEWdf)){
  NEWdf[index,2]<- gsub("DOID:","",NEWdf[index,2])
}
#give better column names
colnames(NEWdf)<-c("disease name","DOID","OMIM","Description")
```

#### 4.3 Step three: we need more data
If we want to map with biomaRt, we need to get the Ensembl ID for each OMIM id we have in the database, to do this, we need to download "mim2gene.txt" from OMIM website as I described above. note that the file are seperate by tab.
```R
omim_to_ensemble <- readr::read_delim(file.path("./inst/extdata/", "mim2gene.txt"),
                         delim = "\t",
                         skip = 5,
                         col_names = c("MIM","Type","NCBI",	"HGNC"	,"Ensembl"))
head(omim_to_ensemble)
# A tibble: 6 x 5
#MIM Type                          NCBI HGNC  Ensembl
#<dbl> <chr>                        <dbl> <chr> <chr>  
#  1 100050 predominantly phenotypes        NA NA    NA     
#2 100070 phenotype                100329167 NA    NA     
#3 100100 phenotype                       NA NA    NA     
#4 100200 predominantly phenotypes        NA NA    NA     
#5 100300 phenotype                       NA NA    NA     
#6 100500 moved/removed                   NA NA    NA  
```
#### 4.4 Step four: finally let's map
Since we only to map those OMIM id that are in our NEWdf, we first want to merge omim_to_ensemble and NEWdf on OMIM id, then get a list of the corresponding Ensembl ID, then use this list and biomaRt to get the map. However, after the merge, there are 155 elememts have a HGNC symbol,but only 3 of them have Ensemble ID. We may not need biomaRt to do the mapping since there are only 3 entries, but I'll show it anyway.
```R
#merge data on OMIM id, note two table have different column names on that
mergedata<-merge.data.frame(omim_to_ensemble,NEWdf,by.x = "MIM",by.y = "OMIM")

#how many of them have a ensemble id?
nrow(mergedata[!is.na(mergedata$Ensembl), ])#3 rows

#How many of them have HGNC sym?
nrow(mergedata[!is.na(mergedata$HGNC), ])#155 
HGNCmap<-mergedata[!is.na(mergedata$HGNC), ]
#how may out in HGNCmap have Ensembl ID?
nrow(HGNCmap[!is.na(HGNCmap$Ensembl), ]) #still 3

usefuldata<-mergedata[!is.na(mergedata$Ensembl), ]
#MIM           Type  NCBI HGNC.x         Ensembl
#273  151430 gene/phenotype   596   BCL2 ENSG00000171791
#1427 602715           gene 25802  LMOD1 ENSG00000163431
#1633 606501           gene 54545 MTMR12 ENSG00000150712...not complete

myvalue<-unique(usefuldata$Ensembl)
myMart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
mymap <- biomaRt::getBM(filters = "ensembl_gene_id",
                      attributes = c("ensembl_gene_id",
                                     "hgnc_symbol"),
                      values = myvalue,
                      mart = myMart)
#> mymap
#ensembl_gene_id hgnc_symbol
#1 ENSG00000150712      MTMR12
#2 ENSG00000163431       LMOD1
#3 ENSG00000171791        BCL2
```
#### 4.5 Validation
Since there are only three DOID that can be mapped to HGNC symbol, only by inspection, we can see the mapping does not have duplicate, and it is a one-to-one relationship.But we still want to compare it with our reference table.
```R
#check with our reference table
load("./inst/extdata/HGNC.RData")
for (index in 1:nrow(NEWdf)){
  if (NEWdf[index,3]%in% HGNCtoOMIM$OMIMIDc){
    print(sprintf("OMIM: %s",NEWdf[index,3]))
  }
}
#[1] "OMIM: 151430"
#[1] "OMIM: 602715"
#[1] "OMIM: 606501"
#give usefuldata a better name
HDOsys<-usefuldata
save(HDOsys,file = file.path("inst","extdata","HDOsys.RData"))
# can be loaded with
load(file = file.path("inst", "extdata", "HDOsys.RData"))
```

# 5 Annotating 
Our database is basicly a DAG object, each node is a disease, each node has a unique DOID number, and it has a "is_a" realationship with their parents.
first we want to know how many nodes are there in the databse.<br/>
From the result of readOboFile we know that there are total 11514 node, and among these, only 3389(29.7%) of them have OMIM id.<br/>
Within the OMIM databse, 15614 out of 26135(59.7%) of them have Ensemble ID. 62.4% of them have HGNC symbol,so some of the OMIM id have a HGNC symbol but no Ensemble ID.<br/>
After merge DOID with OMIM id and OMIM database, we have 3284 entries, 3389-3284=105,which means,weiredly, some of the OMIM id shown in HDO database does not shown in OMIM dabatase. let's have a look what are they:
```R
> setdiff(mergedata$MIM,NEWdf$OMIM)
numeric(0)
> setdiff(NEWdf$OMIM,omim_to_ensemble$MIM)
  [1] "PS267700" "PS275200" "PS107970" "PS276900" "PS310500" "PS162400" "PS124900" "PS220290" "PS304500" "PS210600"
 [11] "PS608594" "PS266600" "PS208500" "PS604348" "PS308350" "PS208085" "PS168000" "PS213300" "300000"   "PS257920"
 [21] "PS600118" "PS268310" "PS272430" "PS156200" "PS309510" "PS250950" "PS253310" "PS220210" "PS605552" "PS242300"
 [31] "PS604772" "PS600513" "PS145980" "PS308240" "PS608808" "PS214700" "PS312080" "PS214450" "PS169150" "PS602014"
 [41] "PS254130" "PS211600" "PS243300" "PS251200" "500000"   "PS610805" "PS601419" "PS604004" "PS192600" "PS109730"
 [51] "PS119580" "PS219000" "PS242860" "PS183600" "PS120100" "PS147950" "PS601198" "PS614039" "PS604931" "PS159000"
 [61] "PS253600" "PS234200" "PS113900" "PS252150" "PS118220" "PS133100" "PS603075" "PS151660" "PS122470" "PS310300"
 [71] "PS166200" "PS256100" "PS115200" "PS603278" "PS130000" "PS191100" "PS607634" "PS227650" "PS125310" "PS256730"
 [81] "PS204000" "PS108800" "PS135900" "PS209900" "PS305100" "PS104500" "PS235200" "PS303350" "PS215100" "PS256300"
 [91] "PS127550" "PS137800" "PS161800" "PS603165" "PS105400" "PS163950" "PS601462" "PS601678" "PS605389" "PS212065"
[101] "PS145600" "PS300751" "PS173900" "PS193500" "PS244400"

```
Most of them start with PS,I searched but didn't find any formal explaination about this,since I do not know what happened here,I think it will be better just leave it as it is.(I don't want to remove the prefix "PS" then do the mapping, since PS605389 may not mean 605389, and the result will be incorrect)<br/>
So,only 3 out of 11514 (0.03%) are successfully mapped to HGNC symbol, which is a very low rate, however, give consideration of the nature of my database, which is just disease name, it make sense. This result is saying only very small proportion(0.03%) of diseases can be directly linked to single gene/protein, most of the disease are a result of combination influences of many genes. <br/>

# 6 Annotation of the example gene set
Sadly, none of the element in the test gene set was able to mapped to our DOID, we cannot perform any kind of annotation around it.
```R
xSet <- c("AMBRA1", "ATG14", "ATP2A1", "ATP2A2", "ATP2A3", "BECN1", "BECN2",
          "BIRC6", "BLOC1S1", "BLOC1S2", "BORCS5", "BORCS6", "BORCS7",
          "BORCS8", "CACNA1A", "CALCOCO2", "CTTN", "DCTN1", "EPG5", "GABARAP",
          "GABARAPL1", "GABARAPL2", "HDAC6", "HSPB8", "INPP5E", "IRGM",
          "KXD1", "LAMP1", "LAMP2", "LAMP3", "LAMP5", "MAP1LC3A", "MAP1LC3B",
          "MAP1LC3C", "MGRN1", "MYO1C", "MYO6", "NAPA", "NSF", "OPTN",
          "OSBPL1A", "PI4K2A", "PIK3C3", "PLEKHM1", "PSEN1", "RAB20", "RAB21",
          "RAB29", "RAB34", "RAB39A", "RAB7A", "RAB7B", "RPTOR", "RUBCN",
          "RUBCNL", "SNAP29", "SNAP47", "SNAPIN", "SPG11", "STX17", "STX6",
          "SYT7", "TARDBP", "TFEB", "TGM2", "TIFA", "TMEM175", "TOM1",
          "TPCN1", "TPCN2", "TPPP", "TXNIP", "UVRAG", "VAMP3", "VAMP7",
          "VAMP8", "VAPA", "VPS11", "VPS16", "VPS18", "VPS33A", "VPS39",
          "VPS41", "VTI1B", "YKT6")
any(xSet%in% usefuldata$HGNC) #FALSE
```
# 7 Reference and citation
* course website
* example STRING Package:https://github.com/hyginn/BCB420.2019.STRING
* parse function by ed-lau [https://github.com/Molecular-Proteomics/roboparse/blob/master/main.R]
* Lynn M Schriml,Elvira Mitraka, James Munro, Becky Tauber etc. Human Disease Ontology 2018 update: classification, content and workflow expansion, Nucleic Aicds Reaseach,Volumn 47,Issue D1,8 Jan 2019,pages D955-d962[https://academic.oup.com/nar/article/47/D1/D955/5165342]

 
# 8 Acknowlegments
Thanks to my classmates who answered my questions when I am confusing.<br/>
Thanks ed-lau for whoever you are, you saved me a lot of time.


&nbsp;

&nbsp;

<!-- END -->

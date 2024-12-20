# monkeypox_ins
a collection of codes for monkeypox phylodinamics in a global context

## the R code ## 
```r
setwd("/home/vjimenez/Documentos/monkey_research/monkey_230123")
getwd()
dir()

library(ape)
library(seqinr)
library(lubridate)
library(tidyr)

############################################
##         change peruvian headers        ##
############################################

f <- read.fasta("peru_1.fasta")
f1 <- names(f)
f1[1:6]
f2 <- data.frame(NETLAB=f1,order=1:length(f1))
head(f2)

n <- read.csv("newnames.tsv", header=TRUE, sep="\t")
head(n)
d <- merge(f2,n,by="NETLAB",all.x=TRUE)
dim(d)
head(d)

e <- d[order(d$order),]
head(e)

write.fasta(f,"peru_2.fasta",names=e$gisaid)

############################################
## merge NETLAB metadata - GISAID metadata ##
############################################

setwd("/home/vjimenez/Documentos/monkey_research/monkey_270123")
getwd()
dir()

## merge NETLAB metadata - GISAID metadta # 

a <- read.csv("gisaid_pox_2023_01_25_14.tsv",header=TRUE, sep="\t")

a$Location <- gsub(" ","",a$Location)
head(a)
b <- separate(a,"Location",c("continent","country","province","city","locality"),sep="/")
b <- data.frame(strain=b$Virus.name,Accesion=b$Accession.ID,date=ymd(b$Collection.date), 
                continent=b$continent,country=b$country,province=b$province,
                city=b$city,host=b$Host,sex=b$Gender,age=b$Patient.age)
dim(b)

pm <- read.csv("peru_1.metadata.tsv",header=TRUE, sep="\t")
unique(pm$FECHA.TOMA.DE.MUESTRA)
pm2 <- data.frame(strain=pm$gisaid,Accesion=pm$CODIGO.NETLAB,date=dmy(pm$FECHA.TOMA.DE.MUESTRA), 
                  continent=rep("South America",dim(pm)[1]),country=rep("Peru",dim(pm)[1]),province=pm$departamento,
                  city=rep("",dim(pm)[1]),host=rep("Human",dim(pm)[1]),sex=pm$SEXO,age=pm$EDAD)

c <- rbind(pm2,b)
dim(c)

write.table(c,"metadata_1.tsv", sep="\t", row.names = FALSE)

###############################################
## merge GENERAL metadata and NEXTCLADE results## 
###############################################
a <- read.csv("metadata_1.tsv", header=TRUE, sep="\t")
dim(a)
names(a)

b <- read.csv("nextclade.tsv", header=TRUE, sep="\t")
dim(b)
names(b)
b1 <- data.frame(strain=b$seqName, lineage=b$lineage,clade=b$clade, coverage=b$coverage, missing=b$totalMissing)
head(b1)

c <- merge(a,b1, by="strain", all.x=FALSE)
dim(c)
head(c)

plot(c$coverage, c$missing)
c$missing <- as.numeric(c$missing)
c$coverage <- as.numeric(c$coverage)

d <- c[c$missing < 20000 & c$coverage > 0.90 , ]
dim(d)
plot(d$coverage, d$missing)

as.data.frame(table(d$country))
names(d)

d$country <- gsub("DemocraticRepublicoftheCongo","Democratic Republic of the Congo",d$country)
d$country <- gsub("UnitedKingdom","United Kingdom",d$country)
d$continent <- gsub("SouthAmerica","South America",d$continent)
d$continent <- gsub("NorthAmerica","North America",d$continent)
d$continent <- gsub("NorthAmeric","North America",d$continent)
d$continent <- gsub("Asian","Asia",d$continent)
as.data.frame(table(d$country))
as.data.frame(table(d$continent))
names(d)

e <- d[!is.na(d$strain),]
dim(d)
dim(e)

write.table(e,"metadata_2.tsv", row.names=FALSE, sep="\t")
write.csv(d$strain,"ext.txt",row.names=FALSE)

a <- read.csv("metadata_2.tsv", header=TRUE, sep="\t")
dim(a)
names(a)

l <- read.csv("lat_longs.tsv", sep="\t", header=FALSE)
dim(l)
names(l)
head(l)

setdiff(unique(a$country),l$V2)
setdiff(l$V2,unique(a$country))

###########################
### generate final data 1 ###
###########################

f <- read.csv("list.txt", header=FALSE)
head(f)
a <- read.csv("metadata_2.tsv", header=TRUE, sep="\t")
head(a)
dim(a)

a1 <- a[a$strain %in% f$V1 , ]
dim(a1)
head(a1)
unique(a1$continent)

write.table(a1,"metadata_final.tsv", row.names=FALSE, sep="\t")

###########################
### generate final data 2 ###
###########################
dir()
a2 <- read.csv("metadata_final.tsv", header = TRUE, sep="\t")
dim(a2)
a3 <- d[!is.na(a2$date),]
dim(a3)

f <- read.csv("list.txt", header=FALSE)
head(f)
dim(f)

setdiff(f$V1,a3$strain)
length(intersect(f$V1,a3$strain))

a4 <- a3[a3$strain %in% f$V1 , ]
dim(a4)
head(a4)
write.table(a4,"metadata_final_2.tsv", sep="\t", row.names=FALSE)
write.table(a4$strain,"list2.txt", sep="\t", row.names=FALSE)
```

## The bash code ##
```r
# edit tip names ##
sed 's/|.*//g' input.fasta > input_unaligned.fasta ;
grep ">" input_unaligned.fasta ; 
grep ">" input_unaligned.fasta | wc -l ;
ls -lh ; 

# run nextclade #

nextclade dataset get --name hMPXV --output-dir DB
nextclade run --input-dataset DB/ --output-all out input_unaligned.fasta
mv out/nextclade.aligned.fasta out/nextclade.tsv .
mv nextclade.aligned.fasta sequences_total.fasta

## extract selected sequences ##
seqtk subseq sequences_total.fasta ext.txt > sequences.fasta
seqtk subseq sequences_total.fasta list2.txt > sequences_1.fasta
seqtk subseq snps2.fasta list2.txt > snps3.fasta

## run snp-sites ## 

snp-sites -o snps.fasta sequences.fasta ; 
aliview snps.fasta ; 

# run raxml-nextsrain #

raxmlHPC-PTHREADS -p 123568 -m GTRCAT -s snps.phy -T 31 -# 10 -n nwk ; 
mv RAxML_bestTree.nwk raw_tree.nwk ;
mkdir raxml ; 
mv RAx* raxml/ ; 
conda activate nextstrain ; 
augur refine --alignment sequences.fasta --tree raw_tree.nwk --metadata metadata.tsv --output-tree refine_tree.nwk --output-node-data node_Data.json --timetree --coalescent opt --gen-per-year 52 --root least-squares --covariance --date-confidence --date-inference marginal --branch-length-inference marginal --year-bounds 2017 2022 --divergence-units mutations-per-site --seed 12548746 ; 
augur ancestral --tree refine_tree.nwk --alignment sequences.fasta --output-node-data ancestral.json --inference marginal --keep-overhangs ; 
augur translate --tree refine_tree.nwk --ancestral-sequences ancestral.json --reference-sequence ref_monkey.gb --output-node-data translate.json ; 
augur traits --tree refine_tree.nwk --metadata metadata.tsv --columns continent country province lineage age sex clade --confidence --output-node-data traits.json ; 
AUGUR_RECURSION_LIMIT=10000 augur export v2 --tree refine_tree.nwk --node-data ancestral.json node_Data.json traits.json translate.json --output final.json --auspice-config auspice_config.json --geo-resolutions country --color-by-metadata continent country province lineage age sex clade --panels tree map entropy frequencies --metadata metadata.tsv --lat-longs lat_longs.tsv ; 
ls -lh ; 

# only phyplogeny # 

raxmlHPC-PTHREADS -p 123568 -m GTRCAT -s snps.phy -T 31 -# 10 -n nwk ; 
mv RAxML_bestTree.nwk raw_tree.nwk ;
mkdir raxml ; 
mv RAx* raxml/ ; 
ls -lh ; 
```

# RNAseq analysis using WDL on Terra
This documents introduced all the files or scripts needed to run the WDL on Terra.



# metadata

The metadata contains

- one inputs.json
- one outputs.json
- one or more entity_tables
- one workspace_table
- other index or reference files



## The inputs.json  and the outputs.json

The inputs.json  and the outputs.json are uploaded and downloaded in the WORKFLOWS panel. The inputs.json  and the outputs.json can also be produced by a R script. The R script help distinguish entity variable and workspace variable. Entity variables have the prefix "`this`", while Workspace variables have the prefix "`workspace`".



## entity_tables

The entity_tables are written by R script to make the information manipulation easier. entity_tables are uploaded and downloaded in the DATA panel.



## Workspace _table

The workspace table stores the key-value information that are shared by all the samples or workflows within a workspace. The workspace table are usually compiled in MS Excel spreadsheet, and then transpose into a table and then uploaded to the `Terra/WorkSpace/Data/Other Data/Workspace Data`.



## Other index or reference files

A other file that is used to load index and/or reference files is usually a text file that stores the locations or paths of the index or reference files. This text file will usually read in by the WDL script so the locations or paths mentioned in the text file are loaded as `File type` in WDL.



### Hisat Index file



### What are the contents in the `genome_index.txt` file?

```shell
gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c/grch38_snp_tran/genome_snp_tran.1.ht2
gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c/grch38_snp_tran/genome_snp_tran.2.ht2
gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c/grch38_snp_tran/genome_snp_tran.3.ht2
gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c/grch38_snp_tran/genome_snp_tran.4.ht2
gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c/grch38_snp_tran/genome_snp_tran.5.ht2
gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c/grch38_snp_tran/genome_snp_tran.6.ht2
gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c/grch38_snp_tran/genome_snp_tran.7.ht2
gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c/grch38_snp_tran/genome_snp_tran.8.ht2
```

### How to generate the `genome_index.txt` a file like this?

Go to [HiSat2 website](https://ccb.jhu.edu/software/hisat2/index.shtml) to download the index files of your desired version, say, [grch38_snp_tran](https://cloud.biohpc.swmed.edu/index.php/s/grch38_snp_tran/download) for RNAseq reads alignment. 

Upload the index files to the google drive bucket in the Terra Workspace. Then use Google Cloud SDK to list the paths of index files and copy the paths of the files into a text file called `genome_index.txt`. If you use MicroSoft Windows system, ensure that you convert the newline symbols in Windows System into according newline symbols in Linux, perhaps using notepad++. If you don't convert the newline symbols, the regex pattern of Terra Cromwell will break.

```shell
D:\Program Files\Google\Cloud SDK>gsutil ls gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c/grch38_snp_tran

gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c/grch38_snp_tran/genome_index.txt
gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c/grch38_snp_tran/genome_snp_tran.1.ht2
gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c/grch38_snp_tran/genome_snp_tran.2.ht2
gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c/grch38_snp_tran/genome_snp_tran.3.ht2
gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c/grch38_snp_tran/genome_snp_tran.4.ht2
gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c/grch38_snp_tran/genome_snp_tran.5.ht2
gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c/grch38_snp_tran/genome_snp_tran.6.ht2
gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c/grch38_snp_tran/genome_snp_tran.7.ht2
gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c/grch38_snp_tran/genome_snp_tran.8.ht2
gs://fc-6ca14431-c7c8-4f2b-91d1-20aed069eb8c/grch38_snp_tran/make_grch38_snp_tran.sh
```



### GTF annotation files

We use the latest version for GTF files [Homo_sapiens.GRCh38.99.gtf.gz](ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz) from [ensembl.org](https://uswest.ensembl.org/downloads.html).



# Workflow Overview

Measures gene level expression in HT-Seq raw read count。 Following the [documentation of HTSeq-count](https://htseq.readthedocs.io/en/release_0.11.1/count.html),  we do run the following commands:

```python
htseq-count [options] <alignment_files> <gff_file>
```


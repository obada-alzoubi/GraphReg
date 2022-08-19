#!/usr/bin/env Rscript
#######################################
# Script to run HiCDCPlus.
#######################################
args = commandArgs(trailingOnly=TRUE)
library(HiCDCPlus)
organism        = 'human'         # human/mouse
genome          = 'hg38'           # hg19/hg38/mm10
assay_type      = 'HiC'            # HiC/HiChIP/MicroC/HiCAR
cell_line       = args[1]          # GM12878/K562/hESC/mESC,..
binsize         = as.integer(args[2])
qval            = as.double(args[3])   
data_path       = args[4] 
hicfile_path    = args[5]
ncore           = 24 # number of cores

# fdr as string 
if (qval==0.01)
{
  fdr='01'
}else if (qval==0.001)
{
  fdr='001'
}else if (qval==0.1)
{
  fdr='1'
}
# bin Size 
if (binsize==5000){
  res = '5kb'            # 5kb/10kb
}else if(binsize==10000)
{
  res = '10kb'            # 5kb/10kb
}

outdir          = paste0(data_path,'/data/', cell_line, '/hic/HiC/')
# Create folder for output
dir.create(outdir,recursive = T)
# List of chroms 
chrs = paste0('chr', seq(22))

for (chr in chrs){
  # 0.1/0.01/0.001
  cat('Chr ', chr, '\n')
  out_hic_bed = paste0(data_path, '/data/', cell_line, '/hic/',
                       assay_type,'/',cell_line,'_',assay_type,'_FDR_',fdr,'_',chr)
  print(out_hic_bed)
  construct_features(output_path=paste0(outdir, "/", "hg38_", cell_line,
                                        "_5kb_GATC_", chr),
                     gen="Hsapiens",gen_ver=genome,
                     sig="GATC",
                     bin_type="Bins-uniform",
                     binsize=binsize,
                     chrs=c(chr))
  
  #generate gi_list instance
  gi_list<-generate_bintolen_gi_list(gen="Hsapiens",gen_ver=genome,
                                     bintolen_path=paste0(outdir,"/hg38_",
                                                          cell_line,
                                                          "_5kb_GATC_", chr,
                                                          "_bintolen.txt.gz"))
  
  #add .hic counts
  gi_list<-add_hic_counts(gi_list,hic_path = hicfile_path)
  
  #expand features for modeling
  gi_list<-expand_1D_features(gi_list)
  #run HiC-DC+ on 2 cores
  set.seed(1010) #HiC-DC downsamples rows for modeling
  gi_list<-HiCDCPlus_parallel(gi_list,ncore=ncore)
  head(gi_list)
  
  gi_list_write(gi_list,fname=out_hic_bed,
                significance_threshold=qval,rows='significant')
  
}
# ouput
# GraphReg - Expanded  
This repo is an extention for GraphReg. It has detialed explnations and improvment over the orignal repo.I added demo exmaples and fixed several bugs and hard coded values. It is still an ongoing process. 
### Top Configurations
Before moving on with running. You need to pay attention to the following paramters. 
* `bin size`: the bin size at which you will sample genemone. It also depends on HiC resoltion. This parameter will determine many parameters.  
* `Genome Build`: hg19, hg38 or mouse.
* `model`: sequence or epigenomic.
* `cell lines`: cell line(s) for extracting epi features. 
Each script should be examined for each of the previous parameters. There are many hard coded paramters in scripts!
### Compuational Time of data preprocessing 
Some scripts are memory greedy. You may need a computer with 72 GB at minimum.

<img
  src="assets/distal_enhancers.png"
  alt="Alt text"
  title=""
  style="display: inline-block; margin: 0 auto; max-width: 300px">

**GraphReg** ([Chromatin interaction aware gene regulatory modeling with graph attention networks](https://genome.cshlp.org/content/32/5/930.short)) is a graph neural network based gene regulation model which integrates DNA sequence, 1D epigenomic data (such as chromatin accessibility and histone modifications), and 3D chromatin conformation data (such as Hi-C, HiChIP, Micro-C, HiCAR) to predict gene expression in an informative way. **GraphReg** is a versatile model which can be used to answer interesting questions in regulatory genomics such as:

- How well we can predict expression of a gene by using the epigenomic features of its promoter and candidate enhancers and enhancer-promoter interactions? Can this model be used in unseen cell types to predict gene expression?

- What are the cis regulatory elements of the genes in each cell type? Which candidate enhancers are functional and play a role in gene regulation?

- Which transcription factor (TF) motifs are important for gene regulation? How do distal TF motifs regulate their target genes?

This repository contains all the codes for training **GraphReg** models and all the downstream analyses for gene expression prediction, enhancer validation, and discovering regulating TF motifs.

## Data preparation
Probably, this will be the most time consuming approach. I chose datasets for elaboration purposes only. 

### Reference Genome 
We need the reference genome. In this example, I'm using hg38. If you want to use hg19, please downlaod the hg19 reference genome from [Gencode](https://www.gencodegenes.org/human/)

```
mkdir -p data\genome 
!wget -O data/genome/GRCh38.primary_assembly.genome.fa.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz
!gunzip data/genome/GRCh38.primary_assembly.genome.fa.gz

```
### 1D data (epigenomic and CAGE)
We need a coverage file `bigwig` for each epigenomic track. We have used some useful functions from [Basenji](https://github.com/calico/basenji) for reading and writing the `bigwig` files, which can be found in [utils](https://github.com/karbalayghareh/GraphReg/tree/master/utils). 

We can use two different approaches to generate `bigwig` files from alignment `BAM` files:

- [`bam_cov.py`](https://github.com/karbalayghareh/GraphReg/blob/master/utils/bam_cov.py) from Basenji. This works best when we want to work with each cell type individually. The coverage tracks from different cell types are not normalized by this method. In **Epi-GraphReg** if we are interested in cross-cell-type generalization, the coverage tracks should be normalized by other techniques such as DESeq, otherwise there would be batch effect between cell types due to sequencing depths, which would hurt the generalization performance. 

- [`bamCoverage`](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) from [deepTools](https://deeptools.readthedocs.io/en/develop/index.html). This is more suitable for cross-cell-type analyses, as they offer some normalization methods for `bigwig` files. In particular, we use 1x normalization or reads per genome coverage (RPGC), which normalizes the coverage in each bin by sequencing depth. We run `bamCoverage` with bin size 100 for epigenomic tracks and 5000 for CAGE-seq.

After generating the `bigwig` files, we use [data_read.py](https://github.com/karbalayghareh/GraphReg/blob/master/utils/data_read.py) to read the `bigwig` files and save the coverage signals in `hdf5` format. We use `pool_width = 100` (to get the coverage in 100bp bins) for epigenomic tracks and `pool_width = 5000` (to get the coverage in 5Kb bins) for CAGE. The reason of using 5Kb bins for CAGE is that we use 5Kb resolution of 3D assays and want to have corresponding bins. If we use `bam_cov.py` to generate `bigwig` files, we set `sum_stat = 'sum'` to sum all the base-pair coverage in each bin; otherwise, if we use `bamCoverage` to generate `bigwig` files, we set `sum_stat = 'max'` as the coverage per bin has already been computed per bin. 
We will use K562 cell line to get CAGE data: https://www.encodeproject.org/experiments/ENCSR000CJN/. Then, we need to index the file using [samtools](http://www.htslib.org/). Then, we will be able to run `bam_cov.py` or `bamCoverage`. 
I'm going to use `bamCoverage` in this case. Thus, we need to install deeptools. It is very important to pay attention to the bin size `5000 bp`
```
conda install -c bioconda samtools # install samtools
conda install -c bioconda deeptools # install deepTools

```
Now, we can run indexing and coverage. 
```
wget -O data/ENCFF366MWI.bam https://encode-public.s3.amazonaws.com/2016/08/04/1c730238-8266-4ace-ba49-47fcd56838a8/ENCFF366MWI.bam
samtools index data/ENCFF366MWI.bam
bamCoverage -b data/ENCFF366MWI.bam \
-bs 5000 \ # Bin size 
--normalizeUsing RPGC \
-o data/K562_CAGE_binsize_5000bp.bigWig \
--effectiveGenomeSize 2913022398 # hg38
# Data should be orginized in data\cell\bam\nam.bigWig 
mkdir -p  data/K562/bam/
cp  data/K562_CAGE_binsize_5000bp.bigWig data/K562/bam/
rm data/K562_CAGE_binsize_5000bp.bigWig 
```
Now, you can run `data_read.py`. Before you run the code, make sure to change the following parts in the code. 
```
  organism = 'human'
  cell_line = 'K562'
  res = '5kb'
  genome='hg38'
  data_path = 'parent_of_data' # the parent of data directory. For example, if data is located in /home/codes/GraphReg/data, then data_path='/home/codes/GraphReg'

```
```
python GraphReg/utils/data_read.py 
```

### 3D data (chromatin conformation: Hi-C/HiChIP/Micro-C/HiCAR)
*Note: You need R for this step*
The chromatin conformation `fastq` data from various 3D assays such as Hi-C, HiChIP, Micro-C, HiCAR could be aligned to any genome (using packages like [Juicer](https://github.com/aidenlab/juicer) or [HiC-Pro](https://github.com/nservant/HiC-Pro)) to get `.hic` files. **GraphReg** needs connecivity graphs for each chromosome. As these 3D data are very noisy, we need some statistical tools to get the significant interactions for the graphs, otherwise it would be very noisy. To this end, we use [HiCDCPlus](https://github.com/mervesa/HiCDCPlus) which gets the `.hic` files and returns the significance level (FDR) for each genomic interaction (of resolution 5Kb) based on a Negative Binomial model. We filter the interactions and keep the ones with `FDR <= alpha` to form the graphs and adjacency matrices. We have worked with three different values of `alpha = 0.1, 0.01, 0.001` and noticed that its ideal value depends on the 3D data. But, we recommend `alpha = 0.1` as a default and less stringent cutoff. 
First, install `HiCDCPlus` pckage.  

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("HiCDCPlus")
```
Download data and create directories

```
!wget -O data/4DNFIW6H9U3S.hic https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/ac58fc15-48c2-4eec-a689-23b677b4b6e7/4DNFIW6H9U3S.hic'
mkdir -p data/IMR90/hic/HiC/
```
Then, use the following Rscript code to generate the significance level (FDR) . 
```
library(HiCDCPlus)
cell_line  = 'IMR90'           # GM12878/K562/hESC/mESC/IMR90
organism   = 'human'           # human/mouse
res        = '5kb'              # 5kb/10kb
binsize    = 5000
genome     = 'hg38'                # hg19/hg38/mm10
assay_type = 'HiC'        # HiC/HiChIP/MicroC/HiCAR
qval       = 0.01   
fdr        = '01'
data_path  = 'GraphReg'
hicfile_path    = 'data/4DNFIW6H9U3S.hic'
outdir = 'data/IMR90/hic/HiC/'

chrs = paste0('chr', seq(22))
chr='chr21'
for (chr in chrs){
  # 0.1/0.01/0.001
  cat('Chr ', chr, '\n')
  out_hic_bed = paste0(data_path, '/data/', cell_line, '/hic/',
                       assay_type,'/',cell_line,'_',assay_type,'_FDR_',fdr,'_',chr)
  print(out_hic_bed)
  construct_features(output_path=paste0(outdir, "/", "hg38_IMR90_5kb_GATC_", chr),
                     gen="Hsapiens",gen_ver=genome,
                     sig="GATC",
                     bin_type="Bins-uniform",
                     binsize=binsize,
                     chrs=c(chr))
  
  #generate gi_list instance
  gi_list<-generate_bintolen_gi_list(gen="Hsapiens",gen_ver=genome,
    bintolen_path=paste0(outdir,"/hg38_IMR90_5kb_GATC_", chr, "_bintolen.txt.gz"))
  
  #add .hic counts
  gi_list<-add_hic_counts(gi_list,hic_path = hicfile_path)
  
  #expand features for modeling
  gi_list<-expand_1D_features(gi_list)
  #run HiC-DC+ on 2 cores
  set.seed(1010) #HiC-DC downsamples rows for modeling
  gi_list<-HiCDCPlus_parallel(gi_list,ncore=24)
  head(gi_list)
  

  gi_list_write(gi_list,fname=out_hic_bed,
               significance_threshold=qval,rows='significant')

}
```
We also need to generate the bed files of ranges.
```
import numpy as np
import pandas as pd
import time
import os
import scipy.sparse


##### write seqs #####
'''
T = 400 
TT = T + T//2
organism = 'human'
genome='hg38'
data_path = 'parent_of_data' # the parent of data directory. For example, if data is located in /home/codes/GraphReg/data, then data_path='/home/codes/GraphReg'
for i in range(1,22+1):
    chr = 'chr'+str(i)
    filename_seqs = data_path+'/data/csv/seqs_bed/'+organism+'/'+genome+'/5kb/sequences_'+chr+'.bed'
    seq_dataframe = pd.DataFrame(data = [], columns = ["chr", "start", "end"])
    chrom = i

    if organism=='human' and genome=='hg19':
       chr_len = [249235000, 243185000, 197960000, 191040000, 180905000, 171050000, 159125000, 146300000, 141150000, 135520000, 134945000,
                       133840000, 115105000, 107285000, 102520000, 90290000, 81190000, 78015000, 59115000, 62965000, 48115000, 51240000]
    elif organism=='human' and genome=='hg38':
       chr_len = [248950000, 242185000, 198290000, 190205000, 181530000, 170800000, 159340000, 145130000, 138385000, 133790000, 135080000,
                       133270000, 114355000, 107035000, 101985000, 90330000, 83250000, 80365000, 58610000, 64435000, 46700000, 50810000]
    elif organism=='mouse':
       chr_len = [195465000, 182105000, 160030000, 156500000, 151825000, 149730000, 145435000, 129395000, 124590000, 130685000, 122075000, 
               120120000, 120415000, 124895000, 104035000, 98200000, 94980000, 90695000, 61425000]

    nodes_list = []
    for i in range(0, chr_len[chrom-1]+5000, 5000):
        nodes_list.append(i)
    nodes_list = np.array(nodes_list)
    left_padding = np.zeros(TT).astype(int)
    right_padding = np.zeros(TT).astype(int)
    nodes_list = np.append(left_padding, nodes_list)
    nodes_list = np.append(nodes_list, right_padding)
    seq_dataframe['start'] = nodes_list
    seq_dataframe['end'] = nodes_list + 5000
    seq_dataframe['chr'] = chr
    print(seq_dataframe)
    seq_dataframe.to_csv(filename_seqs, index = False, header = False, sep = '\t')
```
The outputs of HiCDCPlus and ranges are given to [hic_to_graph.py](https://github.com/karbalayghareh/GraphReg/blob/master/utils/hic_to_graph.py) to generate the adjacency matrices for each chromosome, which are saved as sparce matrices. 
You need again to edit `hic_to_graph.y`
```

cell_line  = 'K562'           # GM12878/K562/hESC/mESC
organism   = 'human'           # human/mouse
res        = '5kb'                  # 5kb/10kb
genome     = 'hg38'                # hg19/hg38/mm10
assay_type = 'HiC'        # HiC/HiChIP/MicroC/HiCAR
qval       = 0.01                    # 0.1/0.01/0.001 - please make sure the q-value here matches wiht Rscript
data_path  = 'parent_of_data' # the parent of data directory. For example, if data is located in /home/codes/GraphReg/data, then data_path='/home/codes/GraphReg'
```
Also depending on the version of `HiCDCPlus`, you might need to fix this line:
```
hic_dataframe.columns = ["chr", "start_i","end_i","chrj", "start_j","end_j", "D",  "count","pvalue", "qval", "mu", "sdev"]
```
 Instead of 
``` 
hic_dataframe.columns = ["chr", "start_i", "start_j", "qval", "count"]
```

### TSS bins and positions
We need to have a `BED` file for TSS annotations. This file could be extracted from any gene annotation `GTF` files for any genome build. We have used GENCODE annotations which can be found [here](https://www.gencodegenes.org/). The TSS annotation `BED` file is given to [find_tss.py](https://github.com/karbalayghareh/GraphReg/blob/master/utils/find_tss.py) to compute the number of TSS's in each 5Kb bin. `find_tss.py` saves four outputs as numpy files: start position of each bin, number of TSS's in each bin, and the gene names (if existent) and their TSS positions in each bin. With 5Kb bins, the majority of them would have one TSS. However, there is a chance that a bin has 2 or 3 TSS's, in which case we save the first TSS position and all the genes (in the format `gene_name_1+gene_name_2`), because we want to keep track of all the genes appearing in each bin. 

```
wget -O data/gencode.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz
gunzip data/gencode.gtf.gz
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' data/gencode.gtf | gtf2bed -> data/gencode.bed  
mkdir -p data/tss/human/hg38/
mv  data/gencode.bed  data/tss/human/hg38/gencode.v38.annotation.gtf.tss.bed
```
Then, we can run `find_tss.py` after editing the following line:
```
resolution = '5kb'   # 5kb/10kb
organism   = 'human' # human/mouse
genome     = 'hg38   # hg38/hg19/mm10
thr        = 0       # only keep the tss bins whose distance from bin borders are more than "thr" 
                     # (only applicable when want to consider the bins with 1 tss, otherwise thr = 0)
data_path = 'parent_of_data' # the parent of data directory. For example, if data is located in /home/codes/GraphReg/data, then data_path='/home/codes/GraphReg'
```

### Writing all data (1D, 3D, TSS) to TFRecords 
**GraphReg** has been implemented in TensorFlow. [TFRecord](https://www.tensorflow.org/tutorials/load_data/tfrecord) is an efficient format to store and read data in TensorFlow. We use [data_write.py](https://github.com/karbalayghareh/GraphReg/blob/master/utils/data_write.py) to read (1) epigenomic coverage files (saved in `h5` format), (2) sparse adjacency matrices (saved in `npz` format), and (3) TSS files (saved in `np` format) and save them sequentially in TFRecords in each chromosome. We start from beginning of each chromosome and write the epigenomic and CAGE coverages, adjacency matrices, and TSS annotations for the regions of 6Mb. Then we sweep the entire chromosome by steps of 2Mb. This way, there is no overlap for the middle 2Mb regions where we predict gene expression values. For each batch of 6Mb, the dimensions of data would be: `60,000` for each epigenomic track, `1200` for CAGE, and `1200 x 1200` for adjacency matrices. The predicted CAGE values in the middle `400` bins would appear in the loss function so that all the genes could see their distal enhancers up to 2Mb up- and downstream of their TSS. 
Again, we need to edit the following lines in `data_write.py`. You need also to determine the model `model=epi or model=seq`.
```
  organism = 'human'          # human/mouse
  res = '5kb'                 # 5kb/10kb
  cell_line = 'K562'          # K562/GM12878/hESC/mESC
  genome='hg19'               # hg19/hg38/mm10
  model = 'epi'               # seq/epi
  assay_type = 'HiC'        # HiC/HiChIP/MicroC/HiCAR
  qval = 0.01                    # 0.1/0.01/0.001
data_path = 'parent_of_data' # the parent of data directory. For example, if data is located in /home/codes/GraphReg/data, then data_path='/home/codes/GraphReg'
```
You need also to edit the following lines based on input data:
```
seqs_cov_files = [data_path+'/data/'+cell_line+'/seqs_cov/CAGE_cov_RPGC_'+chr_temp+'.h5',
                  #data_path+'/data/'+cell_line+'/seqs_cov/H3K4me3_cov_RPGC_'+chr_temp+'.h5',
                  #data_path+'/data/'+cell_line+'/seqs_cov/H3K27ac_cov_RPGC_'+chr_temp+'.h5',
                  #data_path+'/data/'+cell_line+'/seqs_cov/H3K4me1_cov_FC_'+chr_temp+'.h5',
                  #data_path+'/data/'+cell_line+'/seqs_cov/H3K27me3_cov_FC_'+chr_temp+'.h5',
                  #data_path+'/data/'+cell_line+'/seqs_cov/CTCF_cov_FC_'+chr_temp+'.h5',
                  data_path+'/data/'+cell_line+'/seqs_cov/DNase_cov_RPGC_'+chr_temp+'.h5'
                  ]
```
The last step in data preperation is writing data as tfrecorde. Please not this process will take long time (~ 40 min per chromosome).
```
mkdir -p data/tfrecords
python data_write.py
```
The TFRecord files are slightly different for **Epi-GraphReg** and **Seq-GraphReg** models: (1) TFRecords for **Seq-GraphReg** also contain one-hot-coded DNA sequences of the size `6,000,000 x 4`, as the DNA sequence is an input for these models, (2) The epigenomic signals for **Epi-GraphReg** undergo an extra log-normalization, via function `log2(x+1)`, to reduce their dynamic ranges, as they are inputs in  **Epi-GraphReg** models.

Now that we have generated TFRecord files, we are ready to train the models.

## Training GraphReg models

### Epi-GraphReg

Use [`Epi-GraphReg.py`](https://github.com/karbalayghareh/GraphReg/blob/master/train/Epi-GraphReg.py) to train the **Epi-GraphReg** models. You should specify the validation and test chromosomes. The remaining autosomal chromosomes are used for training. For example:
```
python Epi-GraphReg.py -c K562 -p $data_path -a HiChIP -g 1 -q 0.1 -v 3,13 -t 4,14
```
trains **Epi-GraphReg** on cell line K562, using graphs extracted from HiChIP with FDR (q-value) cutoff 0.1, in generalizable mode `-g 1`, with Chrs 3 and 13 as the  validation and Chrs 4 and 14 as the test chromosomes. `$data_path` is the directory where TFRecords have been stored. Training on generalizable mode means that the model uses the normalized epigenomic coverage tracks (for example the ones obtained from RPGC normalization) so that the trained model can be used in other cell types as well.

### Seq-GraphReg

The **Seq-GraphReg** models can be trained in two ways: (1) end-to-end, (2) separate. End-to-end training means that both epigenomic and CAGE data are predicted in a multi-task learning fashion. However, separate training means that we train two tasks (epigenomic and GACE) separately: we first use CNN layers to predict the epigenomic tracks from DNA sequence (similar to Basenji) and then feed the bottleneck representations to the graph attention layers to predict the CAGE values. So, which one should you use? It depends on the amount of GPU memory the users have access to. End-to-end training requires high GPU memory as it needs to load the entire 6Mb DNA sequence to the GPU memory. One advantage of an end-to-end model is the ability to do gradient-based feature attribution from output of any gene back to the base pair level. However, if a high GPU memory is not available, we can employ separate training, where we can use smaller genomic regions of length 100Kb (instead of 6Mb) to predict the epigenomic data from DNA sequences as these are local features and no graph is used for this task. Then after predicting the entire 6Mb (60 mini batches of 100Kb), we concatenate their corresponding bottleneck representations (with the size `60,000 x 64`, where 64 is the dimension of bottleneck representations) and feed that to the graph attention layers along with the corresponding graph of 6Mb region. 

#### End-to-end training

Use [Seq-GraphReg_e2e.py](https://github.com/karbalayghareh/GraphReg/blob/master/train/Seq-GraphReg_e2e.py) to train end-to-end **Seq-GraphReg** models. You should specify the validation and test chromosomes. The remaining autosomal chromosomes are used for training. For example:
```
python Seq-GraphReg_e2e.py -c K562 -p $data_path -a HiChIP -q 0.1 -v 3,13 -t 4,14
```
trains end-to-end **Seq-GraphReg** on cell line K562, using graphs extracted from HiChIP with FDR (q-value) cutoff 0.1 with Chrs 3 and 13 as the validation and Chrs 4 and 14 as the test chromosomes.

#### Separate training

1- Use [Seq-CNN_base.py](https://github.com/karbalayghareh/GraphReg/blob/master/train/Seq-CNN_base.py) to first train a CNN model to predict the epigenomic data from DNA sequence. For example:
```
python Seq-CNN_base.py -c K562 -p $data_path -v 3,13 -t 4,14
```
trains a CNN model on cell line K562 with Chrs 3 and 13 as the validation and Chrs 4 and 14 as the test chromosomes. We have implemented four flavors of such epigenomic CNN models: with/without dilation layers, and with/without FFT (Fast Fourier Transform) loss (so overall four models), which can be found [here](https://github.com/karbalayghareh/GraphReg/tree/master/train). The idea of adding FFT attribution prior is borrowed from [here](https://proceedings.neurips.cc/paper/2020/file/1487987e862c44b91a0296cf3866387e-Paper.pdf) to improve the interpretability of the deep learning models using DNA sequences.

2- Use [Seq-GraphReg.py](https://github.com/karbalayghareh/GraphReg/blob/master/train/Seq-GraphReg.py) to train the graph attention networks of **Seq-GraphReg** to predict the CAGE values. For example:
```
python Seq-GraphReg.py -c K562 -p $data_path -a HiChIP -q 0.1 -f 1 -d 0 -v 3,13 -t 4,14
```
first loads the epigenomic CNNs trained in step (1) on cell line K562 without dilation layers `-d 0` and with FFT loss `-f 1`. Then it uses the bottleneck representation as the input for the second task: predicting CAGE using the graphs extracted from HiChIP with FDR cutoff 0.1. Note that the validation (Chrs 3 and 13) and test (Chr 4 and 14) chromosomes are the same as step (1).

## Testing GraphReg models

After training the **GraphReg** models, it is time to use them to predict gene expression (CAGE) in held-out test chromosomes or cell types. All the scripts for saving the predictions and plotting the results are in [test](https://github.com/karbalayghareh/GraphReg/tree/master/test). We first run the prediction scripts (explained below) to save the CAGE predictions and all meta data (such as number of enhancer-promoter interactions, name of the genes, TSS positions, etc.) for the TSS bins in the test chromosomes in a `csv` file. Then we run the plotting scripts (which call the saved `csv` files) to plot the results. 

- For **Epi-GraphReg**, run [Epi-GraphReg_test_multiple_runs.py](https://github.com/karbalayghareh/GraphReg/blob/master/test/Epi-GraphReg_test_multiple_runs.py) to save predictions in the test chromosomes, and run [plot_results_epi_models.py](https://github.com/karbalayghareh/GraphReg/blob/master/test/plot_results_epi_models.py) to plot the results. 

- For generalizable (cross cell type) **Epi-GraphReg**, run [Epi-GraphReg_generalizable_test_multiple_runs.py](https://github.com/karbalayghareh/GraphReg/blob/master/test/Epi-GraphReg_generalizable_test_multiple_runs.py) to save predictions in the test chromosomes, and run [plot_results_epi_models_generalizable.py](https://github.com/karbalayghareh/GraphReg/blob/master/test/plot_results_epi_models_generalizable.py) to plot the results. Note that we recommend cross cell type and cross chromosome predictions, meaning that the chromosomes that we use in the test cell type for predictions are the same chromosomes that are held out as the test chromosomes in the training cell type. This makes sure that the model does not just use the memorized values of the training cell types in the test cell type for the genes that do not vary a lot between the cell types. 

- For end-to-end **Seq-GraphReg**, run [Seq-GraphReg_e2e_test_multiple_runs.py](https://github.com/karbalayghareh/GraphReg/blob/master/test/Seq-GraphReg_e2e_test_multiple_runs.py) to save predictions in the test chromosomes, and run [plot_results_seq_models.py](https://github.com/karbalayghareh/GraphReg/blob/master/test/plot_results_seq_models.py) to plot the results.

- For separate **Seq-GraphReg**, run [Seq-GraphReg_test_multiple_runs.py](https://github.com/karbalayghareh/GraphReg/blob/master/test/Seq-GraphReg_test_multiple_runs.py) to save predictions in the test chromosomes, and run [plot_results_seq_models.py](https://github.com/karbalayghareh/GraphReg/blob/master/test/plot_results_seq_models.py) to plot the results.

## Feature attributions of GraphReg models

### Enhancer validation

We can use feature attributions of **GraphReg** models to find out which regions are important for gene expression, meaning that they could be considered as the enhancers of any target gene. To validate enhancer predictions, we need experimental enhancer perturbation data. We have used two such datasets in the K562 cell line: (1) CRISPRi FlowFISH from [this](https://www.nature.com/articles/s41588-019-0538-0) paper and (2) Targeted Perturb-seq (TAP-seq) from [this](https://www.nature.com/articles/s41592-020-0837-5) paper. We can use the feature attributions of both **Epi-GraphReg** and **Seq-GraphReg** models. We have tried [DeepSHAP](https://github.com/slundberg/shap) and gradient-by-input for feature attribution. Note that in **Epi-GraphReg** we do the gradients of output with respect to the input epigenomic bins (of size 100bp), while in **Seq-GraphReg** we do the gradients of outputs with respect to the bottleneck representation bins (of size 100bp). This approach has been suggested by Basenji. All the scripts for these kind of analyses are in [here](https://github.com/karbalayghareh/GraphReg/tree/master/feature_attribution).

### In-silico TF motif knockout and ISM

We can use **Seq-GraphRe** models to get insights about how TF motifs regulate their target genes. To have experimental validation data, we downloaded CRISPRi TF knockout experiments in K562 cells from [ENCODE](https://www.encodeproject.org/search/?searchTerm=ENCSR016WFQ&limit=all). These data show which genes are down- (up) regulated after knockout of a TF. To see if we can predict such effects, we delete the motifs of that TF in the genome and then observe the difference in the predictions of each gene. If a TF motif is important for the expression of a gene, the predicted value of gene expression would go down. By doing so, we can find the direct effects of each TF on their target genes. However, this prediction is not perfect as there is no way for the model to consider the indirect effects of the TFs on their target genes. [Seq-GraphReg_TF_KO.py](https://github.com/karbalayghareh/GraphReg/blob/master/feature_attribution/Seq-GraphReg_TF_KO.py) does these kinds of analyses. 

We can also perform an in-silico saturation mutagenesis (ISM) in some distal enhancers of the genes to get an idea about how single mutations get affect the expression predictions of the genes. [Seq-GraphReg_ISM.py](https://github.com/karbalayghareh/GraphReg/blob/master/feature_attribution/Seq-GraphReg_ISM.py) can do ISM using **Seq-GraphReg** models. 



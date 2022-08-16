# Attention - GraphReg - Expanded  
This repo is an extention for GraphReg. It has detialed explnations and improvment over the orignal repo.I added deom exmaples and fixing several bugs and hard coded values. It is still an ongoing process. 

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

### 3D data (chromatin conformation: Hi-C/HiChIP/Micro-C/HiCAR)
The chromatin conformation `fastq` data from various 3D assays such as Hi-C, HiChIP, Micro-C, HiCAR could be aligned to any genome (using packages like [Juicer](https://github.com/aidenlab/juicer) or [HiC-Pro](https://github.com/nservant/HiC-Pro)) to get `.hic` files. **GraphReg** needs connecivity graphs for each chromosome. As these 3D data are very noisy, we need some statistical tools to get the significant interactions for the graphs, otherwise it would be very noisy. To this end, we use [HiCDCPlus](https://github.com/mervesa/HiCDCPlus) which gets the `.hic` files and returns the significance level (FDR) for each genomic interaction (of resolution 5Kb) based on a Negative Binomial model. We filter the interactions and keep the ones with `FDR <= alpha` to form the graphs and adjacency matrices. We have worked with three different values of `alpha = 0.1, 0.01, 0.001` and noticed that its ideal value depends on the 3D data. But, we recommend `alpha = 0.1` as a default and less stringent cutoff. 

The outputs of HiCDCPlus is given to [hic_to_graph.py](https://github.com/karbalayghareh/GraphReg/blob/master/utils/hic_to_graph.py) to generate the adjacency matrices for each chromosome, which are saved as sparce matrices. 

### TSS bins and positions
We need to have a `BED` file for TSS annotations. This file could be extracted from any gene annotation `GTF` files for any genome build. We have used GENCODE annotations which can be found [here](https://www.gencodegenes.org/). The TSS annotation `BED` file is given to [find_tss.py](https://github.com/karbalayghareh/GraphReg/blob/master/utils/find_tss.py) to compute the number of TSS's in each 5Kb bin. `find_tss.py` saves four outputs as numpy files: start position of each bin, number of TSS's in each bin, and the gene names (if existent) and their TSS positions in each bin. With 5Kb bins, the majority of them would have one TSS. However, there is a chance that a bin has 2 or 3 TSS's, in which case we save the first TSS position and all the genes (in the format `gene_name_1+gene_name_2`), because we want to keep track of all the genes appearing in each bin. 

### Writing all data (1D, 3D, TSS) to TFRecords 
**GraphReg** has been implemented in TensorFlow. [TFRecord](https://www.tensorflow.org/tutorials/load_data/tfrecord) is an efficient format to store and read data in TensorFlow. We use [data_write.py](https://github.com/karbalayghareh/GraphReg/blob/master/utils/data_write.py) to read (1) epigenomic coverage files (saved in `h5` format), (2) sparse adjacency matrices (saved in `npz` format), and (3) TSS files (saved in `np` format) and save them sequentially in TFRecords in each chromosome. We start from beginning of each chromosome and write the epigenomic and CAGE coverages, adjacency matrices, and TSS annotations for the regions of 6Mb. Then we sweep the entire chromosome by steps of 2Mb. This way, there is no overlap for the middle 2Mb regions where we predict gene expression values. For each batch of 6Mb, the dimensions of data would be: `60,000` for each epigenomic track, `1200` for CAGE, and `1200 x 1200` for adjacency matrices. The predicted CAGE values in the middle `400` bins would appear in the loss function so that all the genes could see their distal enhancers up to 2Mb up- and downstream of their TSS. 

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



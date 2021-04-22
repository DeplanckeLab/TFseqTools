![](https://img.shields.io/badge/build-passing-green.svg)
![](https://img.shields.io/badge/version-1.0.0-blue.svg)
![](https://img.shields.io/badge/picard-2.9.0-blue.svg)
![](https://img.shields.io/badge/java-1.8-red.svg)

# TF-seq Tools
A suite of tools for the pre-processing of TF-seq data

## 1. Preprocessing
TF-seq generates an enriched library containing one to one mapping between cell barcodes and TF barcodes.
The first step is to create a read/UMI count matrix containing TF x cell to properly assign TF to cells.

### 1.1. Download software
TF-seq preprocessing tool is provided as a [single executable jar file](../master/releases/TFseqTools-1.0.jar?raw=true).
The .jar file contains all required materials and can be run on any terminal.

#### Dependencies
- **Java version:**
For the tools to run properly, you must have Java >=1.8 installed. 
To check your java version, open your terminal application and run the following command:

```bash
java -version
```

If the output looks something like java version "1.8.x" (or above), you are good to go. 
If not, you may need to update your version; see the [Oracle Java website](http://www.oracle.com/technetwork/java/javase/downloads/) to download the latest JRE (for users) or JDK (for developers).

- **Picard:**
The software relies on [Picard Tools](http://broadinstitute.github.io/picard/), but the Picard JAR is already embedded in the released JAR, so no need to install it yourself.

### 1.2 Sequencing output, and mapping
After sequencing your enriched TF-seq libraries, you should obtain two fastq files: 
- **R1 fastq file:** This file should contain the barcodes and UMIs (for e.g. from 10x or Dropseq). Usually, the barcode comes before the UMI sequence.
- **R2 fastq file:** This file should contain the exact same number of reads (and read names) than the R1 file. Except that the sequence of the reads are the sequences representing the RNA fragments.

#### Alignment to vector
You need to align this enriched library to the default construct containing the TF barcodes. Let's call this assembly "Vector". We provide the TF-seq vector fasta file in the "genome" folder.
You can align the R2 fastq file on the " Vector" assembly using any aligner you are used to (STAR, bwa, ...).
Example using STAR:

```bash
STAR --runMode alignReads --outSAMmapqUnique 60 --runThreadN 8 --genomeDir example/genome/Vector/STAR_Index --outFilterMultimapNmax 1 --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix example/bam/ --readFilesIn example/fastq/TF_enrich_example_R2.fastq.gz
```

### 1.3 Running TF-seq preprocessing tool to generate the TF-Cell mapping matrix
To check that the preprocessing tool is working properly, run the following command:

```bash
java -jar TFseqTools.jar
```
As shown in the output of this command, this tool allows the user to run one possible tool:

```
TFseqTools vx.x

-- Options --
	Counter		Count occurence of each TF in Drop-seq/10x R2 aligned BAM and R1 fastq
```

Now, you can run the "Counter" tool to see its parameters:

```bash
java -jar TFseqTools.jar Counter
```

Which will output something like that:

```
TFseqTools vx.x

-- 'Counter' options --
	--r1 %s 	[Required] Path of R1 FastQ file.
	--r2 %s 	[Required] Path of R2 aligned BAM file [do not need to be sorted or indexed].
	--tf %s 	[Required] File containing known TF barcodes
	-o %s 		Output folder [default = folder of BAM file]

-- Additional options --
	--log %i 	Detailed log file [default: None]
	--nu %i 	Number of allowed difference (hamming distance) for two UMIs to be counted only once [default = 0].
	-p %s 		Cell barcode pattern/order found in the reads of the R1 FastQ file. Barcode names should match the barcode file [default = 'BU', i.e. barcode followed by the UMI].
				'B' [Required] is used for specifying the barcode position.
				'U' can be used for specifying a UMI value position.
				'?' can be used to ignore specific nucleotides.
	--UMI %i 	If your barcode pattern contains UMI ('U'), you should specify this parameter as the length of the UMI [e.g. 10x run is 10]
	--BC %i 	If your barcode pattern contains Barcode ('B'), you should specify this parameter as the length of the barcode [e.g. 10x run is 16]
```

You can see that the main options to use are --r1, --r2 and --tf
So the final command will be:
```
java -jar TFseqTools.jar Counter --r1 example/fastq/TF_enrich_example_R1.fastq.gz --r2 example/bam/Aligned.out.bam --tf example/barcodes/tf_barcodes.txt
```

> **Note:** You can download/edit this **[example of tf barcodes file](../master/example/barcodes/tf_barcodes.txt)**

## Directory content
* **src**: all source files required for compilation
* **lib**: all JAR dependencies required for compilation / execution
* **releases**: final, all embedded, released .jar files
* **example**: list of example files and outputs, to set up the pipeline

## Author
Vincent Gardeux - vincent.gardeux@epfl.ch

# genomePileup

## Description
This program generates training data for [deepore](https://github.com/kishwarshafin/deePore/) a neural network based variant caller. This program generates chromosome wise training data taking the VCF file as ground truth.

## Dependencies
python3.5+, pysam, numpy, scipy, PIL

## Usage
### Whole genome run
```./pileupProcessor.sh [bam_file_path] [ref_fasta_file_path] [vcf_file_path] [output_directory] [number_of_threads]```

#### Parameters: <br/>
-- `bam_file_path`: Path to alignment bam file <br/>
-- `ref_fasta_file_path`: Path to reference fasta file <br/>
-- `vcf_file_path`: Path to vcf file <br/>
-- `output_direcotry`: Path to a directory where output will be saved <br/>
-- `number_of_threads`: Most threads it can use to generate images <br/>

#### Example run:
```./pileupProcessor.sh ~/illumina/chr3.bam ~/illumina/chr3.fa  ~/illumina/chr3_whole_chr.vcf.gz ~/pileup_output/test/ 8```

---


### Specific site run using the python script
```python3 main.py --bam [bam_file_path] --ref [ref_fasta_file_path] --vcf [vcf_file_path] --output_dir [output_directory]  --parallel [bool] --max_threads [int] --contig [string] --site [string]```

#### Parameters <br/>
-- `bam_file_path`: Path to alignment bam file <br/>
-- `ref_fasta_file_path`: Path to reference fasta file <br/>
-- `vcf_file_path`: Path to vcf file <br/>
-- `output_direcotry`: Path to a directory where output will be saved <br/>
-- `parallel`: If true it will use mutiprocessing (Default is False) <br/>
-- `max_threads`: If parallel is true then will use max_thread number of threads (Default is 5) <br/>
-- `contig`: Contig to focus (Default is "chr3") [Example: chr3, chr2 etc.]<br/>
-- `site`: Site of the contig to focus (Default is empty) [Example: ":100000-200000"] <br/>

#### Example run:
```python3 main.py --bam ~/illumina/chr3.bam --ref ~/illumina/chr3.fa --vcf ~/illumina/chr3_whole_chr.vcf.gz --output_dir ~/pileup_output/test/ --contig chr3 --site :100000-200000 --parallel True```

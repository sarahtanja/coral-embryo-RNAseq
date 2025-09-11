# Step 2: QAQC RNA sequences
Sarah Tanja
October 10, 2024

- [<span class="toc-section-number">1</span> Background](#background)
  - [<span class="toc-section-number">1.1</span> fastqc](#fastqc)
  - [<span class="toc-section-number">1.2</span> multiqc](#multiqc)
  - [<span class="toc-section-number">1.3</span> fastp](#fastp)
- [<span class="toc-section-number">2</span> Create a Bash variables
  file](#create-a-bash-variables-file)
  - [<span class="toc-section-number">2.1</span> Name directory
    paths](#name-directory-paths)
  - [<span class="toc-section-number">2.2</span> Program
    paths](#program-paths)
  - [<span class="toc-section-number">2.3</span> Filename
    patterns](#filename-patterns)
  - [<span class="toc-section-number">2.4</span> Thread
    count](#thread-count)
  - [<span class="toc-section-number">2.5</span> File name
    arrays](#file-name-arrays)
  - [<span class="toc-section-number">2.6</span> Program
    array](#program-array)
- [<span class="toc-section-number">3</span> Quality check raw reads
  with `FastQC`](#quality-check-raw-reads-with-fastqc)
- [<span class="toc-section-number">4</span> Compile a report of raw
  reads with `MultiQC`](#compile-a-report-of-raw-reads-with-multiqc)
- [<span class="toc-section-number">5</span> Interpretation of `MultiQC`
  report for raw reads](#interpretation-of-multiqc-report-for-raw-reads)
- [<span class="toc-section-number">6</span> Trim & Clean reads using
  `fastp`](#trim--clean-reads-using-fastp)
- [<span class="toc-section-number">7</span> Quality check trimmed reads
  with `FastQC`](#quality-check-trimmed-reads-with-fastqc)
- [<span class="toc-section-number">8</span> Compile report of trimmed
  reads with `MultiQC`](#compile-report-of-trimmed-reads-with-multiqc)
- [<span class="toc-section-number">9</span> Summary & Next
  Steps](#summary--next-steps)

# Background

In this script, we will generate FastQC/MultiQC for raw sequences,
conduct trimming and cleaning, then generate reports for cleaned
sequences.

## [fastqc](https://anaconda.org/bioconda/fastqc)

*FastQC generates sequence quality information of your reads*

``` bash
/home/shared/FastQC-0.12.1/fastqc -h
```


                FastQC - A high throughput sequence QC analysis tool

    SYNOPSIS

        fastqc seqfile1 seqfile2 .. seqfileN

        fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam] 
               [-c contaminant file] seqfile1 .. seqfileN

    DESCRIPTION

        FastQC reads a set of sequence files and produces from each one a quality
        control report consisting of a number of different modules, each one of 
        which will help to identify a different potential type of problem in your
        data.
        
        If no files to process are specified on the command line then the program
        will start as an interactive graphical application.  If files are provided
        on the command line then the program will run with no user interaction
        required.  In this mode it is suitable for inclusion into a standardised
        analysis pipeline.
        
        The options for the program as as follows:
        
        -h --help       Print this help file and exit
        
        -v --version    Print the version of the program and exit
        
        -o --outdir     Create all output files in the specified output directory.
                        Please note that this directory must exist as the program
                        will not create it.  If this option is not set then the 
                        output file for each sequence file is created in the same
                        directory as the sequence file which was processed.
                        
        --casava        Files come from raw casava output. Files in the same sample
                        group (differing only by the group number) will be analysed
                        as a set rather than individually. Sequences with the filter
                        flag set in the header will be excluded from the analysis.
                        Files must have the same names given to them by casava
                        (including being gzipped and ending with .gz) otherwise they
                        won't be grouped together correctly.
                        
        --nano          Files come from nanopore sequences and are in fast5 format. In
                        this mode you can pass in directories to process and the program
                        will take in all fast5 files within those directories and produce
                        a single output file from the sequences found in all files.                    
                        
        --nofilter      If running with --casava then don't remove read flagged by
                        casava as poor quality when performing the QC analysis.
                       
        --extract       If set then the zipped output file will be uncompressed in
                        the same directory after it has been created. If --delete is 
                        also specified then the zip file will be removed after the 
                        contents are unzipped. 
                        
        -j --java       Provides the full path to the java binary you want to use to
                        launch fastqc. If not supplied then java is assumed to be in
                        your path.
                       
        --noextract     Do not uncompress the output file after creating it.  You
                        should set this option if you do not wish to uncompress
                        the output when running in non-interactive mode.
                        
        --nogroup       Disable grouping of bases for reads >50bp. All reports will
                        show data for every base in the read.  WARNING: Using this
                        option will cause fastqc to crash and burn if you use it on
                        really long reads, and your plots may end up a ridiculous size.
                        You have been warned!
                        
        --min_length    Sets an artificial lower limit on the length of the sequence
                        to be shown in the report.  As long as you set this to a value
                        greater or equal to your longest read length then this will be
                        the sequence length used to create your read groups.  This can
                        be useful for making directly comaparable statistics from 
                        datasets with somewhat variable read lengths.

        --dup_length    Sets a length to which the sequences will be truncated when 
                        defining them to be duplicates, affecting the duplication and
                        overrepresented sequences plot.  This can be useful if you have
                        long reads with higher levels of miscalls, or contamination with
                        adapter dimers containing UMI sequences.

                        
        -f --format     Bypasses the normal sequence file format detection and
                        forces the program to use the specified format.  Valid
                        formats are bam,sam,bam_mapped,sam_mapped and fastq
                        

        --memory        Sets the base amount of memory, in Megabytes, used to process 
                        each file.  Defaults to 512MB.  You may need to increase this if
                        you have a file with very long sequences in it.
                    
        --svg           Save the graphs in the report in SVG format.

        -t --threads    Specifies the number of files which can be processed
                        simultaneously.  Each thread will be allocated 250MB of
                        memory so you shouldn't run more threads than your
                        available memory will cope with, and not more than
                        6 threads on a 32 bit machine
                      
        -c              Specifies a non-default file which contains the list of
        --contaminants  contaminants to screen overrepresented sequences against.
                        The file must contain sets of named contaminants in the
                        form name[tab]sequence.  Lines prefixed with a hash will
                        be ignored.

        -a              Specifies a non-default file which contains the list of
        --adapters      adapter sequences which will be explicity searched against
                        the library. The file must contain sets of named adapters
                        in the form name[tab]sequence.  Lines prefixed with a hash
                        will be ignored.
                        
        -l              Specifies a non-default file which contains a set of criteria
        --limits        which will be used to determine the warn/error limits for the
                        various modules.  This file can also be used to selectively 
                        remove some modules from the output all together.  The format
                        needs to mirror the default limits.txt file found in the
                        Configuration folder.
                        
       -k --kmers       Specifies the length of Kmer to look for in the Kmer content
                        module. Specified Kmer length must be between 2 and 10. Default
                        length is 7 if not specified.
                        
       -q --quiet       Suppress all progress messages on stdout and only report errors.
       
       -d --dir         Selects a directory to be used for temporary files written when
                        generating report images. Defaults to system temp directory if
                        not specified.
                        
    BUGS

        Any bugs in fastqc should be reported either to simon.andrews@babraham.ac.uk
        or in www.bioinformatics.babraham.ac.uk/bugzilla/
                       
        

## [multiqc](https://anaconda.org/bioconda/multiqc)

*Multiqc summarizes FastQC analysis logs and summarizes results in an
html report. Checkout the* [git developer version](Installation)

## [fastp](https://anaconda.org/bioconda/fastp)

*FastP provides fast all-in-one preprocessing for FastQ files*

Reference bash variables from Roberts Lab [Code
Snippets](https://robertslab.github.io/resources/code_Snippets/)

How to do a [for loop in
bash](https://stackoverflow.com/questions/9612090/how-to-loop-through-file-names-returned-by-find)

> [!TIP]
>
> Instead of \``git add –all` execute
> `find . -type f -size -104M -exec git add {} +` in terminal before
> making a git commit to only track files that are under the git size
> limit of 10MB.

# Create a Bash variables file

This allows usage of Bash variables across R Markdown chunks.

## Name directory paths

``` bash
{
echo "#### Assign Variables ####"
echo ""

echo "# Data directories"
echo 'export project_dir=/media/4TB_JPG_ext/stanja/gitprojects/coral-embryo-ecotox'
echo 'export output_dir_top=${project_dir}/output'
echo 'export output_dir_qaqc=${project_dir}/output/03_qaqc'
echo 'export output_dir_raw=${project_dir}/output/03_qaqc/raw_reads'
echo 'export output_dir_trim=${project_dir}/output/03_qaqc/trimmed_reads'
echo 'export output_dir_fastqc=${output_dir_trim}/fastqc'
echo 'export raw_reads_dir=${project_dir}/rawfastq'
echo ""
} > .bashvars

cat .bashvars
```

## Program paths

``` bash
{
echo "# Paths to programs"
echo 'export programs_dir="/home/shared"'
echo 'export fastp="${programs_dir}/fastp"'
echo 'export fastqc=${programs_dir}/FastQC-0.12.1/fastqc'
echo 'export multiqc=/home/sam/programs/mambaforge/bin/multiqc'
echo ""
} >> .bashvars

cat .bashvars
```

## Filename patterns

``` bash
{
echo "# Set FastQ filename patterns"
echo "export fastq_pattern='*.fastq.gz'"
echo "export R1_fastq_pattern='*_R1_*.fastq.gz'"
echo "export R2_fastq_pattern='*_R2_*.fastq.gz'"
echo "export trimmed_fastq_pattern='*fastp-trim.fq.gz'"
echo ""
} >> .bashvars

cat .bashvars
```

## Thread count

``` bash
{
echo "# Set number of CPUs to use"
echo 'export threads=40'
echo ""
} >> .bashvars

cat .bashvars
```

## File name arrays

``` bash
{
echo "## Inititalize arrays"
echo 'export fastq_array_R1=()'
echo 'export fastq_array_R2=()'
echo 'export raw_fastqs_array=()'
echo 'export R1_names_array=()'
echo 'export R2_names_array=()'
echo ""
} >> .bashvars

cat .bashvars
```

## Program array

``` bash
{
echo "# Programs associative array"
echo "declare -A programs_array"
echo "programs_array=("
echo '[fastp]="${fastp}" \'
echo '[fastqc]="${fastqc}" \'
echo '[multiqc]="${multiqc}" \'
echo ")"
echo ""
} >> .bashvars

cat .bashvars
```

# Quality check raw reads with `FastQC`

FastQC on Raven lives: `/home/shared/FastQC/fastqc`

``` bash
# Load bash variables into memory
source .bashvars

${fastqc} -h
```

``` bash
# Load bash variables into memory
source .bashvars

# Populate array with FastQ files
fastq_array=(${raw_reads_dir}/*.fastq.gz)

# Pass array contents to new variable
fastqc_list=$(echo "${fastq_array[*]}")

# Run FastQC
# NOTE: Do NOT quote ${fastqc_list}
${fastqc} \
--threads ${threads} \
--outdir ${output_dir_raw} \
${fastqc_list}
```

# Compile a report of raw reads with `MultiQC`

``` bash
# Load bash variables into memory
source .bashvars

cd ../output/02_qaqc/raw_reads
${multiqc} ./
```

# Interpretation of `MultiQC` report for raw reads

Watch a quick [6-min
tutorial](https://www.youtube.com/watch?v=qPbIlO_KWN0) on how to
navigate in the MultiQC Report

# Trim & Clean reads using `fastp`

`fastp` all
[options](https://github.com/OpenGene/fastp?tab=readme-ov-file#:~:text=of%201%20~%206.-,all%20options,-usage%3A%20fastp%20%2Di)
can be found in the fastp git README file.

- `--in1` \| Path to forward read input

- `--in2` \| Path to reverse read input

- `--qualified_quality_phred` \| The quality value that a base is
  qualified. Default 15 means phred quality \>=Q15 is qualified. Chosen
  value is 30.

- `--trim_poly_x 6` \| Enable polyX trimming in 3’ ends. Trims any
  repeating sequence (AAAAAA, TTTTTT, CCCCCC, GGGGGG), the default is
  10… here we’ve chosen if it has 6 or more repeating bases to trim them
  off. This is more conservative than 10 or more repeating bases, and
  will result in cleaner sequences that have fewer polyA tails leftover.

- `--detect_adapter_for_pe` \| This automatically detects Illumina
  adapter sequences for paired end (pe) sequences and trims them off.
  Adapters can be trimmed by overlap analysis, however,
  `--detect_adapter_for_pe` will usually result in slightly cleaner
  output than overlap detection alone. This results in a slightly slower
  run time.

- `--thread 16` \| fastp uses up to 16 threads (even if there are more
  available). If this is not specified, it uses 3 as a default.

- `--html` \| The html format report file name

- `--json` \| The json format report file name

- `-out1` \| Path to forward read output

- `--out2` \| Path to reverse read output

``` bash
# Load bash variables into memory
source .bashvars

# Change to raw reads directory
cd "${raw_reads_dir}"

# Create arrays of fastq R1 files and sample names
for fastq in ${R1_fastq_pattern}
do
  fastq_array_R1+=("${fastq}")
  R1_names_array+=("$(echo "${fastq}" | awk -F"_" '{print $1}')")
done

# Create array of fastq R2 files
for fastq in ${R2_fastq_pattern}
do
  fastq_array_R2+=("${fastq}")
  R2_names_array+=("$(echo "${fastq}" | awk -F"_" '{print $1}')")
done

# Create list of fastq files used in analysis
# Create MD5 checksum for reference
if [ ! -f "${raw_reads_dir}"/*raw-fastq-checksums*fastq.gz.md5 ]; then
for fastq in *.gz
  do
    md5sum ${fastq} >> "${output_dir_qaqc}"/checksums/raw-fastq-checksums.md5
  done
fi

# Run fastp on files
# Adds JSON report output for downstream usage by MultiQC
for index in "${!fastq_array_R1[@]}"
do
  R1_sample_name=$(echo "${R1_names_array[index]}")
  R2_sample_name=$(echo "${R2_names_array[index]}")
  ${fastp} \
  --in1 ${fastq_array_R1[index]} \
  --in2 ${fastq_array_R2[index]} \
  --detect_adapter_for_pe \
  --qualified_quality_phred 30 \
  --thread 16 \
  --trim_poly_x 6 \
  --trim_front1 10 \
  --trim_front2 10 \
  --html "${output_dir_trim}"/"${R1_sample_name}".fastp-trim.report.html \
  --json "${output_dir_trim}"/"${R1_sample_name}".fastp-trim.report.json \
  --out1 "${output_dir_trim}"/"${R1_sample_name}"_R1_001.fastp-trim.fq.gz \
  --out2 "${output_dir_trim}"/"${R2_sample_name}"_R2_001.fastp-trim.fq.gz \
  2>> "${output_dir_trim}"/fastp.stderr


  # Generate md5 checksums for newly trimmed files
  cd "${output_dir_trim}"
  md5sum "${R1_sample_name}"_R1_001.fastp-trim.fq.gz > "${R1_sample_name}"_R1_001.fastp-trim.fq.gz.md5
  md5sum "${R2_sample_name}"_R2_001.fastp-trim.fq.gz > "${R2_sample_name}"_R2_001.fastp-trim.fq.gz.md5
  cd -
done
```

# Quality check trimmed reads with `FastQC`

Using a for loop to unzip and run fastqc sequentially… this uses less
CPU. When I try to run them all at once the CPU usage gets bogged down.
This takes longer, but completes in ~5hrs.

``` bash
# Load bash variables into memory
source .bashvars

############ RUN FASTQC ############

# Create array of trimmed FastQs
trimmed_fastqs_array=(${output_dir_trim}/${trimmed_fastq_pattern})

echo "Files in trimmed_fastqs_array: ${trimmed_fastqs_array[@]}"
echo ""
echo "Beginning FastQC on trimmed reads..."
echo ""
echo "Output directory: ${output_dir_trim}/fastqc"
echo ""

# Run FastQC sequentially on each file in the array
for fastq_file in "${trimmed_fastqs_array[@]}"; do
    echo "Processing file: ${fastq_file}"
    ${fastqc} \
    --threads 40 \
    --outdir ${output_dir_trim}/fastqc \
    --quiet \
    "${fastq_file}"
    echo "Completed FastQC for ${fastq_file}"
done

echo ""
echo "FastQC on all trimmed reads complete!"
echo ""

############ END FASTQC ############
```

# Compile report of trimmed reads with `MultiQC`

``` bash

# Load bash variables into memory
source .bashvars

############ RUN MULTIQC ############
echo "Beginning MultiQC on trimmed FastQC..."
echo ""

${multiqc} ${output_dir_fastqc} -o ${output_dir_fastqc}

echo ""
echo "MultiQC on trimmed FastQs complete."
echo ""

############ END MULTIQC ############

echo "Removing FastQC zip files."
echo ""
rm ${output_dir_fastqc}/*.zip
echo "FastQC zip files removed."
echo ""
```

# Summary & Next Steps

Trimmed reads are now ready for alignment in step 4, `04_align`

> [!IMPORTANT]
>
> ### Don’t forget to always rsync backup!
>
>     rsync -avz /media/4TB_JPG_ext/stanja/gitprojects \
>     stanja@gannet.fish.washington.edu:/volume2/web/stanja/ravenbackup

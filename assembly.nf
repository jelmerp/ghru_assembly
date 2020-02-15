#!/usr/bin/env nextflow
// Pipeline version
version = '1.5.5'
/*

========================================================================================
                          Assembly Pipeline
========================================================================================
 Assembly pipeline based on Shovill (Thanks to Torsten Seemann @torstenseemann): https://github.com/tseemann/shovill
 #### Authors
 Anthony Underwood @bioinformant <au3@sanger.ac.uk>
----------------------------------------------------------------------------------------
*/

def versionMessage(){
  log.info"""
    =========================================
     Assembly Pipeline version ${version}
    =========================================
  """.stripIndent()
}

def helpMessage() {
    log.info"""
    Based on Shovill (Thanks to Torsten Seemann @torstenseemann): https://github.com/tseemann/shovill
    Usage:
    The typical command for running the pipeline is as follows:
    Mandatory arguments:
      --output_dir     Path to output directory
      --adapter_file  The path to a fasta file containing adapter sequences to trim from reads
    Optional arguments:
      Either input_dir and fastq_pattern must be specified if using local short reads or accession_number_file if specifying samples in the SRA for which the fastqs will be downloaded
      --input_dir      Path to input directory containing the fastq files to be assembled
      --fastq_pattern  The regular expression that will match fastq files e.g '*{R,_}{1,2}*.fastq.gz'
      --accession_number_file Path to a text file containing a list of accession numbers (1 per line)
      --depth_cutoff The estimated depth to downsample each sample to. If not specified no downsampling will occur
      --careful Turn on the SPAdes careful option which improves assembly by mapping the reads back to the contigs
      --minimum_scaffold_length The minimum length of a scaffold to keep. Others will be filtered out. Default 500
      --minimum_scaffold_depth The minimum depth of coverage a scaffold must have to be kept. Others will be filtered out. Default 3
      --confindr_db_path The path to the confindr database. If not set assumes using Docker image where the path is '/home/bio/software_data/confindr_database'
      --qc_conditions Path to a YAML file containing pass/warning/fail conditions used by QualiFyr (https://gitlab.com/cgps/qualifyr)
      --prescreen_genome_size_check Size in bp of the maximum estimated genome to assemble. Without this any size genome assembly will be attempted
      --prescreen_file_size_check Minumum size in Mb for the input fastq files. Without this any size of file will be attempted (this and prescreen_genome_size_check are mutually exclusive)
      --full_output Output pre_trimming fastqc reports, merged_fastqs and corrected_fastqs. These take up signficant disk space
   """.stripIndent()
}

//  print help if required
params.help = false
if (params.help){
  versionMessage()
  helpMessage()
  exit 0
}

// Show version number
params.version = false
if (params.version){
  versionMessage()
  exit 0
}

/***************** Setup inputs and channels ************************/
// Defaults for configurable variables
params.input_dir = false
params.output_dir = false
params.fastq_pattern = false
params.accession_number_file = false
params.adapter_file = false
params.depth_cutoff = false
params.careful = false
params.minimum_scaffold_length = false
params.minimum_scaffold_depth = false
params.confindr_db_path = false
params.qc_conditions = false
params.prescreen_genome_size_check = false
params.prescreen_file_size_check = false
params.full_output = false

// check if getting data either locally or from SRA
Helper.check_optional_parameters(params, ['input_dir', 'accession_number_file'])

// set up output directory
output_dir = Helper.check_mandatory_parameter(params, 'output_dir') - ~/\/$/
// check that the user has specified an output_dir that differs from the default './.default_output' required for --help to work
if (output_dir =~ /.default_output/) {
  println "You must specifiy an output directory with --output_dir"
  System.exit(1)
}

//  check a pattern has been specified
if (params.input_dir){
  fastq_pattern = Helper.check_mandatory_parameter(params, 'fastq_pattern')
}

//  check an adapter_file has been specified
adapter_file = Helper.check_mandatory_parameter(params, 'adapter_file')
// assign depth_cutoff from params
depth_cutoff = params.depth_cutoff

// set careful 
if (params.careful) {
  careful = true
} else {
  careful = false
}
// assign minimum scaffold length
if ( params.minimum_scaffold_length ) {
  minimum_scaffold_length = params.minimum_scaffold_length
} else {
  minimum_scaffold_length = 500
}
// assign minimum scaffold depth
if ( params.minimum_scaffold_depth ) {
  minimum_scaffold_depth = params.minimum_scaffold_depth
} else {
  minimum_scaffold_depth = 3
}

// set up confindr database path
if ( params.confindr_db_path ) {
  confindr_db_path = params.confindr_db_path
} else {
  // path in Docker image
  confindr_db_path = "/home/bio/software_data/confindr_database"
}

if (params.prescreen_genome_size_check){
  prescreen_genome_size_check = params.prescreen_genome_size_check
} else {
  prescreen_genome_size_check = false
}

if (params.prescreen_file_size_check){
  prescreen_file_size_check = params.prescreen_file_size_check
} else {
  prescreen_file_size_check = 15
}

// set full_output 
if (params.full_output) {
  full_output = true
} else {
  full_output = false
}

// set up read_pair channel
/*
  Creates the `read_pairs` channel that emits for each read-pair a tuple containing
  three elements: the pair ID, the first read-pair file and the second read-pair file
*/


log.info "======================================================================"
log.info "                  GHRU assembly pipeline"
log.info "======================================================================"
log.info "Running version   : ${version}"
if (params.accession_number_file){
  log.info "Accession File    : ${params.accession_number_file}"
} else if (params.input_dir){
  log.info "Fastq inputs      : ${params.input_dir}/${fastq_pattern}"
}
log.info "Adapter file      : ${adapter_file}"
if (params.depth_cutoff){
  log.info "Depth cutoff      : ${depth_cutoff}"
} else{
  log.info "Depth cutoff      : None"
}

log.info "======================================================================"
log.info "Outputs written to path '${output_dir}'"
log.info "======================================================================"
log.info ""

if (params.accession_number_file){
  accession_number_file = params.accession_number_file - ~/\/$/
  // Fetch samples from ENA
  Channel
      .fromPath(accession_number_file)
      .splitText()
      .map{ x -> x.trim()}
      .set { accession_numbers }

  process fetch_from_ena {
    tag { accession_number }
    
    publishDir "${output_dir}/fastqs",
      mode: 'copy',
      saveAs: { file -> file.split('\\/')[-1] }

    input:
    val accession_number from accession_numbers

    output:
    set accession_number, file("${accession_number}/*.fastq.gz") into raw_fastqs

    """
    enaDataGet -f fastq -as /home/bio/.aspera/aspera.ini ${accession_number}
    """
  }
 
} else if (params.input_dir) {
  input_dir = params.input_dir - ~/\/$/
  fastqs = input_dir + '/' + fastq_pattern
  Channel
    .fromFilePairs( fastqs )
    .ifEmpty { error "Cannot find any reads matching: ${fastqs}" }
    .set { raw_fastqs }
}


def find_genome_size(pair_id, mash_output) {
  m = mash_output =~ /Estimated genome size: (.+)/
  genome_size = Float.parseFloat(m[0][1]).toInteger()
  return [pair_id, genome_size]
}

if ( prescreen_genome_size_check ) {
  // Prescreen Genome Size Estimation

  raw_fastqs.into {raw_fastqs_for_genome_size_prescreening; raw_fastqs_for_filtering}

  process pre_screen_genome_size_estimation {
    tag { pair_id }
    
    input:
    set pair_id, file(file_pair)  from raw_fastqs_for_genome_size_prescreening

    output:
    set pair_id, file('mash_stats.out') into pre_screen_mash_output

    """
    kat hist --mer_len 21  --thread 1 --output_prefix ${pair_id} ${file_pair[0]} > /dev/null 2>&1 \
    && minima=`cat  ${pair_id}.dist_analysis.json | jq '.global_minima .freq' | tr -d '\\n'`
    mash sketch -o sketch_${pair_id}  -k 32 -m \$minima -r ${file_pair[0]}  2> mash_stats.out
    """
  }

  pre_screen_mash_output.map { pair_id, file -> find_genome_size(pair_id, file.text) }.into{ genome_size_estimation_for_inclusion; genome_size_estimation_for_exclusion }

  // determine genomes to take for further processing
  included_genomes_based_on_size = genome_size_estimation_for_inclusion.filter { it[1] <= prescreen_genome_size_check }
  raw_fastqs_for_processing = raw_fastqs_for_filtering.join(included_genomes_based_on_size).map { items -> [items[0], items[1]] }

  // write out excluded genomes to file
  excluded_genomes_based_on_size = genome_size_estimation_for_exclusion.filter { it[1] > prescreen_genome_size_check }
  process write_out_excluded_genomes {
    tag { pair_id }
    
    publishDir "${output_dir}/estimated_size_of_excluded_genomes"
    input:
    set pair_id, genome_size from excluded_genomes_based_on_size

    output:
    file("${pair_id}.estimated_genome_size.txt") 

    """
    echo ${genome_size} > ${pair_id}.estimated_genome_size.txt
    """
  }
} else if (prescreen_file_size_check){
  raw_fastqs.into {raw_fastqs_for_genome_size_prescreening; raw_fastqs_for_filtering}

  process determine_fastq_filesize {
    tag { pair_id }
    
    input:
    set pair_id, file(file_pair)  from raw_fastqs_for_genome_size_prescreening

    output:
    set pair_id, stdout into pre_screen_file_sizes_for_output_file_creation, pre_screen_file_sizes_for_inclusion, pre_screen_file_sizes_for_exclusion

    script:
    """
    stat -Lc %s ${file_pair[0]} |  awk '{printf "%f", \$1/1000000}'
    """
  }
  process write_out_filesize_check {
    tag { pair_id }

    input:
    set pair_id, file_size from pre_screen_file_sizes_for_output_file_creation

    output:
    set pair_id, file("${pair_id}.file_size_check.tsv") into file_size_checks

    script:
    """
    echo "file\tsize\n${pair_id}\t${file_size}" > ${pair_id}.file_size_check.tsv
    """
  }
  
  // filter files based on size
  // included genomes
  included_genomes_based_on_file_size = pre_screen_file_sizes_for_inclusion.filter { it[1].toFloat() >= prescreen_file_size_check }

  // excluded genomes
  excluded_genomes_based_on_file_size = pre_screen_file_sizes_for_exclusion.filter { it[1].toFloat() < prescreen_file_size_check }

  raw_fastqs_for_processing = raw_fastqs_for_filtering.join(included_genomes_based_on_file_size).map { items -> [items[0], items[1]] } 
} else {
  raw_fastqs_for_processing = raw_fastqs
}


// duplicate raw fastq channel for trimming and qc
raw_fastqs_for_processing.into {raw_fastqs_for_qc; raw_fastqs_for_trimming; raw_fastqs_for_length_assessment}

// Assess read length and make MIN LEN for trimmomatic 1/3 of this value
process determine_min_read_length {
  tag { pair_id }

  input:
  set pair_id, file(file_pair) from raw_fastqs_for_length_assessment

  output:
  set pair_id, stdout into min_read_length_for_trimming, min_read_length_for_assembly


  """
  gzip -cd ${file_pair[0]} | head -n 400000 | printf "%.0f" \$(awk 'NR%4==2{sum+=length(\$0)}END{print sum/(NR/4)/3}')
  """
}

// Pre-Trimming QC
process qc_pre_trimming {
  tag { pair_id }
  
  if (full_output){
    publishDir "${output_dir}/fastqc/pre_trimming",
      mode: 'copy',
      pattern: "*.html"

  }

  input:
  set pair_id, file(file_pair) from raw_fastqs_for_qc

  output:
  file('*.html')

  """
  fastqc ${file_pair[0]} ${file_pair[1]}
  """
}

min_read_length_and_raw_fastqs = min_read_length_for_trimming.join(raw_fastqs_for_trimming)

// Trimming
process trimming {
  memory '4 GB'
  
  tag { pair_id }
  
  input:
  set pair_id, min_read_length, file(file_pair) from min_read_length_and_raw_fastqs
  file('adapter_file.fas') from adapter_file

  output:
  set pair_id, file('trimmed_fastqs/*.f*q.gz') into trimmed_fastqs_for_qc, trimmed_fastqs_for_correction, trimmed_fastqs_for_genome_size_estimation, trimmed_fastqs_for_base_counting, trimmed_fastqs_for_species_id

  """
  mkdir trimmed_fastqs
  trimmomatic PE -threads 1 -phred33 ${file_pair[0]} ${file_pair[1]} trimmed_fastqs/${file_pair[0]} /dev/null trimmed_fastqs/${file_pair[1]} /dev/null ILLUMINACLIP:adapter_file.fas:2:30:10 SLIDINGWINDOW:4:20 LEADING:25 TRAILING:25 MINLEN:${min_read_length}  
  """
}

// Post-Trimming QC
process qc_post_trimming {
  tag { pair_id }

  publishDir "${output_dir}/fastqc/post_trimming",
    mode: 'copy',
    pattern: "*.html"
  
  input:
  set pair_id, file(file_pair)  from trimmed_fastqs_for_qc

  output:
  file('*.html')
  set pair_id, file("${pair_id}_R1_fastqc.txt"), file("${pair_id}_R2_fastqc.txt") into qc_post_trimming_files
  set file("${r1_prefix}_fastqc_data"), file("${r2_prefix}_fastqc_data") into fastqc_directories

  script:
  r1_prefix = file_pair[0].baseName.replaceFirst(/\\.gz$/, '').split('\\.')[0..-2].join('.')
  r2_prefix = file_pair[1].baseName.replaceFirst(/\\.gz$/, '').split('\\.')[0..-2].join('.')
  """
  fastqc ${file_pair[0]} ${file_pair[1]} --extract
  # rename files
  mv ${r1_prefix}_fastqc/summary.txt ${pair_id}_R1_fastqc.txt
  mv ${r2_prefix}_fastqc/summary.txt ${pair_id}_R2_fastqc.txt

  # move files for fastqc
  mkdir ${r1_prefix}_fastqc_data
  mkdir ${r2_prefix}_fastqc_data
  mv ${r1_prefix}_fastqc/fastqc_data.txt ${r1_prefix}_fastqc_data
  mv ${r2_prefix}_fastqc/fastqc_data.txt ${r2_prefix}_fastqc_data
  """
}

//FastQC MultiQC
process fastqc_multiqc {
  tag { 'multiqc for fastqc' }
  memory { 4.GB * task.attempt }

  publishDir "${output_dir}/quality_reports",
    mode: 'copy',
    pattern: "multiqc_report.html",
    saveAs: { "fastqc_multiqc_report.html" }

  input:
  file(fastqc_directories) from fastqc_directories.collect { it }

  output:
  file("multiqc_report.html")

  script:
  """
  multiqc --interactive .
  """

}
// Species ID
process species_identification {
  tag{pair_id}

  input:
  set pair_id, file(file_pair)  from trimmed_fastqs_for_species_id

  output:
  set pair_id, file("species_investigation*.tsv") into bactinspector_output

  script:
  """
  bactinspector check_species -fq ${file_pair[0]} 
  """

}

// Genome Size Estimation
process genome_size_estimation {
  tag { pair_id }
  
  input:
  set pair_id, file(file_pair)  from trimmed_fastqs_for_genome_size_estimation

  output:
  set pair_id, file('mash_stats.out') into mash_output

  """
  kat hist --mer_len 21  --thread 1 --output_prefix ${pair_id} ${file_pair[0]} > /dev/null 2>&1 \
  && minima=`cat  ${pair_id}.dist_analysis.json | jq '.global_minima .freq' | tr -d '\\n'`
  mash sketch -o sketch_${pair_id}  -k 32 -m \$minima -r ${file_pair[0]}  2> mash_stats.out
  """
}

// channel to output genome size from mash output
mash_output.map { pair_id, file -> find_genome_size(pair_id, file.text) }.into{genome_size_estimation_for_read_correction; genome_size_estimation_for_downsampling}

trimmed_fastqs_and_genome_size = trimmed_fastqs_for_correction.join(genome_size_estimation_for_read_correction).map{ tuple -> [tuple[0], tuple[1], tuple[2]]}

// Read Corection
process read_correction {
  tag { pair_id }
  
  if (full_output){
    publishDir "${output_dir}",
      mode: 'copy',
      pattern: "corrected_fastqs/*.fastq.gz"
  }

  input:
  set pair_id, file(file_pair), genome_size from trimmed_fastqs_and_genome_size

  output:
  set pair_id, file('lighter.out') into read_correction_output
  set pair_id, file('corrected_fastqs/*.fastq.gz') into corrected_fastqs_for_merging, corrected_fastqs_for_contamination_check

  """
  lighter -od corrected_fastqs -r  ${file_pair[0]} -r  ${file_pair[1]} -K 32 ${genome_size}  -maxcor 1 2> lighter.out
  for file in corrected_fastqs/*.cor.fq.gz
  do
      new_file=\${file%.cor.fq.gz}.fastq.gz
      mv \${file} \${new_file}
  done
  """
}

def find_average_depth(pair_id, lighter_output){
  m = lighter_output =~  /.+Average coverage is ([0-9]+\.[0-9]+)\s.+/
  average_depth = Float.parseFloat(m[0][1])
  return [pair_id, average_depth]
}

// Check for contamination
process check_for_contamination {
  tag {pair_id}
  cpus 2

  publishDir "${output_dir}/confindr",
    mode: 'copy',
    saveAs: { file -> "${pair_id}_${file}"}

  input:
  set pair_id, file(file_pair) from corrected_fastqs_for_contamination_check

  output:
  set pair_id, file('confindr_report.csv') into confindr_files

  script:
  if (file_pair[0] =~ /_R1/){ // files with _R1 and _R2
    """
    confindr.py -i . -o . -d ${confindr_db_path} -t 2 -bf 0.025 -b 2 --cross_detail -Xmx 1500m
    """
  } else { // files with _1 and _2
    """
    confindr.py -i . -o . -d ${confindr_db_path} -t 2 -bf 0.025 -b 2 --cross_detail -Xmx 1500m -fid _1 -rid _2
    """  
  }

}


// Estimate total number of bases
process count_number_of_bases {
  tag { pair_id }
  
  input:
  set pair_id, file(file_pair) from trimmed_fastqs_for_base_counting

  output:
  set pair_id, file('seqtk_fqchk.out') into seqtk_fqchk_output

  """
  seqtk fqchk -q 25 ${file_pair[0]} > seqtk_fqchk.out
  """
}

def find_total_number_of_bases(pair_id, seqtk_fqchk_ouput){
  m = seqtk_fqchk_ouput =~ /ALL\s+(\d+)\s/
  total_bases = new BigDecimal(m[0][1]).toLong() * 2 // the *2 is an estimate since number of reads >q25 in R2 may not be the same
  return [pair_id, total_bases]
}
base_counts = seqtk_fqchk_output.map { pair_id, file -> find_total_number_of_bases(pair_id, file.text) }
corrected_fastqs_and_genome_size_and_base_count = corrected_fastqs_for_merging.join(genome_size_estimation_for_downsampling).join(base_counts).map{ tuple -> [tuple[0], tuple[1], tuple[2], tuple[3]]}

// merge reads and potentially downsample
process merge_reads{
  tag { pair_id }

  if (full_output){
    publishDir "${output_dir}",
      mode: 'copy',
      pattern: "merged_fastqs/*.fastq.gz"
  }
  
  input:
  set pair_id, file(file_pair), genome_size, base_count from corrected_fastqs_and_genome_size_and_base_count

  output:
  set pair_id, file('merged_fastqs/*.fastq.gz') into merged_fastqs

  script:
  
  if (depth_cutoff  && base_count/genome_size > depth_cutoff.toInteger()){
    downsampling_factor = depth_cutoff.toInteger()/(base_count/genome_size)
    """
    mkdir downsampled_fastqs
    seqtk sample  -s 12345 ${file_pair[0]} ${downsampling_factor} | gzip > downsampled_fastqs/${file_pair[0]}
    seqtk sample  -s 12345 ${file_pair[1]} ${downsampling_factor} | gzip > downsampled_fastqs/${file_pair[1]}
    flash -m 20 -M 100 -t 1 -d merged_fastqs -o ${pair_id} -z downsampled_fastqs/${file_pair[0]} downsampled_fastqs/${file_pair[1]} 
    """
  } else {
    """
    flash -m 20 -M 100 -t 1 -d merged_fastqs -o ${pair_id} -z ${file_pair[0]} ${file_pair[1]} 
    """
  }

}

 min_read_length_and_raw_fastqs = min_read_length_for_assembly.join(merged_fastqs)

// assemble reads
process spades_assembly {
  memory { 4.GB * task.attempt }
  
  tag { pair_id }

  input:
  set pair_id, min_read_length, file(file_triplet) from min_read_length_and_raw_fastqs


  output:
  set pair_id, file("scaffolds.fasta") into scaffolds

  script:
  
  spades_memory = 4 * task.attempt
  
  if (min_read_length.toInteger() < 25 ) { // this is the read length divided by 3 see trimming step
    kmers = '21,33,43,53'
  }
  else if (min_read_length.toInteger() < 50) {
    kmers = '21,33,43,53,63,75'
  } else {
    kmers = '21,33,43,55,77,99'
  }

  if (careful) {
    """
    spades.py --pe1-1 ${file_triplet[1]} --pe1-2 ${file_triplet[2]} --pe1-m ${file_triplet[0]} --only-assembler --careful -o . --tmp-dir /tmp/${pair_id}_assembly -k ${kmers} --threads 1 --memory ${spades_memory}
    """
  } else {
    """
    spades.py --pe1-1 ${file_triplet[1]} --pe1-2 ${file_triplet[2]} --pe1-m ${file_triplet[0]} --only-assembler -o . --tmp-dir /tmp/${pair_id}_assembly -k ${kmers} --threads 1 --memory ${spades_memory}
    """
  }

}

// filter scaffolds to remove small and low coverage contigs
process filter_scaffolds {
  tag { pair_id }

  input:
  set pair_id, file(scaffold_file) from scaffolds

  output:
  set pair_id, file("${pair_id}.fasta") into scaffolds_for_single_analysis, scaffolds_for_qc
  file("${pair_id}.fasta") into scaffolds_for_combined_analysis
  
  """
  contig-tools filter -l ${minimum_scaffold_length} -c ${minimum_scaffold_depth} -f ${scaffold_file}
  ln -s scaffolds.filter_gt_${minimum_scaffold_length}bp_gt_${minimum_scaffold_depth}.0cov.fasta ${pair_id}.fasta
  """

}

// assess assembly with Quast
process quast {
  tag { pair_id }
    
  publishDir "${output_dir}/quast",
    mode: 'copy',
    pattern: "report.tsv",
    saveAs: { file -> "${pair_id}_quast_" + file.split('\\/')[-1] }

  input:
  set pair_id, file(contig_file) from scaffolds_for_single_analysis

  output:
  set pair_id, file("${pair_id}") into quast_files_for_multiqc
  set pair_id, file("${pair_id}/report.tsv") into quast_files_for_qualifyr

  """
  quast.py ${contig_file} -o .
  mkdir ${pair_id}
  ln -s \$PWD/report.tsv ${pair_id}/report.tsv
  """
}


// assess assembly with Quast but in a single file
process quast_summary {
  tag { 'quast summary' }
  memory { 4.GB * task.attempt }
  
  publishDir "${output_dir}/quast",
    mode: 'copy',
    pattern: "*report.tsv",
    saveAs: { file -> "combined_${file}"}

  input:
  file(contig_files) from scaffolds_for_combined_analysis.collect( sort: {a, b -> a.getBaseName() <=> b.getBaseName()} )

  output:
  file("*report.tsv") optional true

  """
  quast.py ${contig_files} -o .
  """
}



//QUAST MultiQC
process quast_multiqc {
  tag { 'multiqc for quast' }

  publishDir "${output_dir}/quality_reports",
    mode: 'copy',
    pattern: "multiqc_report.html",
    saveAs: { "quast_multiqc_report.html" }

  input:
  file(quast_reports) from quast_files_for_multiqc.collect { it }

  output:
  file("multiqc_report.html")

  script:
  """
  multiqc --interactive .
  """

}

// determine overall quality of sample
if (params.qc_conditions) {


  qc_conditions_yml = file(params.qc_conditions)
  quality_files = qc_post_trimming_files.join(confindr_files).join(quast_files_for_qualifyr).join(scaffolds_for_qc).join(bactinspector_output).join(file_size_checks)
  process overall_quality {
    tag { pair_id }


    publishDir "${output_dir}/assemblies/pass",
      mode: 'copy',
      pattern: 'assemblies/pass/*',
      saveAs: { file -> file.split('\\/')[-1] }

    publishDir "${output_dir}/assemblies/warning",
      mode: 'copy',
      pattern: 'assemblies/warning/*',
      saveAs: { file -> file.split('\\/')[-1] }
    
    publishDir "${output_dir}/assemblies/failure",
      mode: 'copy',
      pattern: 'assemblies/failure/*',
      saveAs: { file -> file.split('\\/')[-1] }

    input:
    file(qc_conditions_yml)
    set pair_id, file(fastqc_report_r1), file(fastqc_report_r2), file(confindr_report), file(quast_report), file(scaffold_file), file(bactinspector_report), file(file_size_check_output) from quality_files

    output:
    file('assemblies/**/*')
    file("${pair_id}.qualifyr.json") into qualifyr_json_files


    """
    # extract min and max genome sizes from bactinspector output, min file size from
    # file_size_check output and replace place holder in conditions file
    MAX_GENOME_LENGTH=\$(cat ${bactinspector_report} | awk -F'\t' 'NR == 2 {print \$8}')
    # if no species match set to 0
    if [ -z \$MAX_GENOME_LENGTH ]; then MAX_GENOME_LENGTH=0; fi
    # add wobble
    MAX_GENOME_LENGTH=\$(echo \$MAX_GENOME_LENGTH  | awk '{printf("%d",  \$1 * 1.1)}')

    MIN_GENOME_LENGTH=\$(cat ${bactinspector_report} | awk -F'\t' 'NR == 2 {print \$9}')
    # if no species match set to 0
    if [ -z \$MIN_GENOME_LENGTH ]; then MIN_GENOME_LENGTH=0; fi
    # add wobble
    MIN_GENOME_LENGTH=\$(echo \$MIN_GENOME_LENGTH  | awk '{printf("%d",  \$1 * 0.9)}')    

    MIN_FILE_SIZE=\$(cat ${file_size_check_output} | awk -F'\t' 'NR == 2 {print \$2}')
  
    sed -i "s/MAX_GENOME_LENGTH/\${MAX_GENOME_LENGTH}/" ${qc_conditions_yml} 
    sed -i "s/MIN_GENOME_LENGTH/\${MIN_GENOME_LENGTH}/" ${qc_conditions_yml}
    sed -i "s/MIN_FILE_SIZE/\${MIN_FILE_SIZE}/" ${qc_conditions_yml}


    result=`qualifyr check -y ${qc_conditions_yml} -f ${fastqc_report_r1} ${fastqc_report_r2} -c ${confindr_report}  -q ${quast_report} -b ${bactinspector_report} -z ${file_size_check_output} -s ${pair_id} 2> ERR`
    return_code=\$?
    if [[ \$return_code -ne 0 ]]; then
      exit 1;
    else
      if [[ \$result == "PASS" ]]; then
        qc_level="pass"
      elif [[ \$result == "WARNING" ]]; then
        qc_level="warning"
      elif [[ \$result == "FAILURE" ]]; then
        qc_level="failure"
      fi
      mkdir -p assemblies/\${qc_level}
      mv ${scaffold_file} assemblies/\${qc_level}/

      if [[ \$result != "PASS" ]]; then
        mv ERR assemblies/\${qc_level}/${pair_id}_qc_result.tsv
      fi
    fi


    # make json file
    qualifyr check -y ${qc_conditions_yml} -f ${fastqc_report_r1} ${fastqc_report_r2} -c ${confindr_report}  -q ${quast_report} -b ${bactinspector_report} -z ${file_size_check_output} -s ${pair_id} -j -o .
    """
  }


  failed_sample_conditions_file = file("${baseDir}/templates/failed_sample_qualifyr_files/conditions_for_failed_samples.yml")
  bactinspector_file = file("${baseDir}/templates/failed_sample_qualifyr_files/bactinspector.tsv")
  confindr_file = file("${baseDir}/templates/failed_sample_qualifyr_files/confindr.csv")
  fastqc_file = file("${baseDir}/templates/failed_sample_qualifyr_files/fastqc.txt")
  file_size_check_file = file("${baseDir}/templates/failed_sample_qualifyr_files/file_size_check.tsv")
  quast_file = file("${baseDir}/templates/failed_sample_qualifyr_files/quast.txt")
  process make_qualifyr_output_for_failed_samples {
    input:
    set pair_id, file_size from excluded_genomes_based_on_file_size

    file(failed_sample_conditions_yml) from failed_sample_conditions_file
    file(bactinspector_report) from bactinspector_file
    file(confindr_report) from confindr_file
    file(fastqc_report) from fastqc_file
    file(file_size_check_report) from file_size_check_file
    file(quast_report) from quast_file

    output:
    file("${pair_id}.qualifyr.json") into file_size_failure_qualifyr_json_files

    script:
    """
    sed -i "s/FILE_SIZE/${file_size}/" ${file_size_check_report}
    sed -i "s/MIN_FILE_SIZE/${prescreen_file_size_check}/" ${failed_sample_conditions_yml}

    # make json file
    qualifyr check -y ${failed_sample_conditions_yml} -f ${fastqc_report} ${fastqc_report} -c ${confindr_report}  -q ${quast_report} -b ${bactinspector_report} -z ${file_size_check_report} -s ${pair_id} -j -o .
    """
  }
  //QualiFyr report

  combined_qualifyr_json_files = qualifyr_json_files.concat(file_size_failure_qualifyr_json_files)

  process qualifyr_report {
    tag { 'qualifyr report' }

    publishDir "${output_dir}/quality_reports",
      mode: 'copy',
      pattern: "qualifyr_report.*"

    input:
    
    file(json_files) from combined_qualifyr_json_files.collect { it }

    output:
    file("qualifyr_report.*")

    script:
    workflow_command = workflow.commandLine.replaceAll('"', '\\\\"')
    """
    qualifyr report -i . -c 'quast.N50,quast.# contigs (>= 1000 bp),quast.Total length (>= 1000 bp),confindr.contam_status,bactinspector.species' -s "Analysis with GHRU Assembly Pipeline version ${version}<br><br>Command line:<br>${workflow_command}"
    """

  }
} else {
  process write_assembly_to_dir {
    tag { pair_id }

    publishDir "${output_dir}",
      mode: 'copy'

    input:
    set pair_id, file(scaffold_file) from scaffolds_for_qc

    output:
    file("assemblies/${scaffold_file}")

    """
    mkdir assemblies
    mv ${scaffold_file} assemblies/
    """

  }
}

workflow.onComplete {
  Helper.complete_message(params, workflow, version)
}

workflow.onError {
  Helper.error_message(workflow)
}

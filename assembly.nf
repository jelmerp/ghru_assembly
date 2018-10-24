#!/usr/bin/env nextflow
// Pipeline version
version = '1.1'
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
      --adapter_sequences  The path to a fasta file containing adapter sequences to trim from reads
    Optional arguments:
      Either input_dir and fastq_pattern must be specified if using local short reads or accession_number_file if specifying samples in the SRA for which the fastqs will be downloaded
      --input_dir      Path to input directory containing the fastq files to be assembled
      --fastq_pattern  The regular expression that will match fastq files e.g '*_{1,2}.fastq.gz'
      --accession_number_file Path to a text file containing a list of accession numbers (1 per line)
      --depth_cutoff The estimated depth to downsample each sample to. If not specified no downsampling will occur
      --minimum_scaffold_length The minimum length of a scaffold to keep. Others will be filtered out. Default 500
      --minimum_scaffold_depth The minimum depth of coverage a scaffold must have to be kept. Others will be filtered out. Default 3
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
params.minimum_scaffold_length = false
params.minimum_scaffold_depth = false

// check if getting data either locally or from SRA
Helper.check_optional_parameters(params, ['input_dir', 'accession_number_file'])

// set up output directory
output_dir = Helper.check_mandatory_parameter(params, 'output_dir') - ~/\/$/
//  check a pattern has been specified
if (params.input_dir){
  fastq_pattern = Helper.check_mandatory_parameter(params, 'fastq_pattern')
}

//  check an adapter_file has been specified
adapter_file = Helper.check_mandatory_parameter(params, 'adapter_file')
// assign depth_cutoff from params
depth_cutoff = params.depth_cutoff
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



// duplicate raw fastq channel for trimming and qc
raw_fastqs.into {raw_fastqs_for_qc; raw_fastqs_for_trimming; raw_fastqs_for_length_assessment}

// Assess read length and make MIN LEN for trimmomatic 1/3 of this value
process determine_min_read_length {
  tag { pair_id }

  input:
  set pair_id, file(file_pair) from raw_fastqs_for_length_assessment

  output:
  set pair_id, stdout into min_read_length 

  """
  gzip -cd ${file_pair[0]} | head -n 400000 | printf "%.0f" \$(awk 'NR%4==2{sum+=length(\$0)}END{print sum/(NR/4)/3}')
  """
}

// Pre-Trimming QC
process qc_pre_trimming {
  tag { pair_id }
  
  publishDir "${output_dir}/qc_pre_trimming",
    mode: 'copy',
    pattern: "*.html"

  input:
  set pair_id, file(file_pair) from raw_fastqs_for_qc

  output:
  file('*.html')

  """
  fastqc ${file_pair[0]} ${file_pair[1]}
  """
}

min_read_length_and_raw_fastqs = min_read_length.join(raw_fastqs_for_trimming)

// Trimming
process trimming {
  tag { pair_id }
  
  input:
  set pair_id, min_read_length, file(file_pair) from min_read_length_and_raw_fastqs
  file('adapter_file.fas') from adapter_file

  output:
  set pair_id, file('trimmed_fastqs/*.fastq.gz') into trimmed_fastqs_for_qc, trimmed_fastqs_for_correction, trimmed_fastqs_for_genome_size_estimation, trimmed_fastqs_for_base_counting

  """
  mkdir trimmed_fastqs
  trimmomatic PE -phred33 ${file_pair[0]} ${file_pair[1]} trimmed_fastqs/${file_pair[0]} /dev/null trimmed_fastqs/${file_pair[1]} /dev/null ILLUMINACLIP:adapter_file.fas:2:30:10 SLIDINGWINDOW:4:20 LEADING:25 TRAILING:25 MINLEN:${min_read_length}  
  """
}

// Post-Trimming QC
process qc_post_trimming {
  tag { pair_id }

  publishDir "${output_dir}/qc_post_trimming",
    mode: 'copy',
    pattern: "*.html"
  
  input:
  set pair_id, file(file_pair)  from trimmed_fastqs_for_qc

  output:
  file('*.html')

  """
  fastqc ${file_pair[0]} ${file_pair[1]}
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
  mash sketch -o /tmp/sketch_${pair_id}  -k 32 -m 3 -r ${file_pair[0]}  2> mash_stats.out
  """
}

def find_genome_size(pair_id, mash_output) {
  m = mash_output =~ /Estimated genome size: (.+)/
  genome_size = Float.parseFloat(m[0][1]).toInteger()
  return [pair_id, genome_size]
}

// channel to output genome size from mash output
mash_output.map { pair_id, file -> find_genome_size(pair_id, file.text) }.into{genome_size_estimation_for_read_correction; genome_size_estimation_for_downsampling}

trimmed_fastqs_and_genome_size = trimmed_fastqs_for_correction.join(genome_size_estimation_for_read_correction).map{ tuple -> [tuple[0], tuple[1], tuple[2]]}
// trimmed_fastqs_and_genome_size.subscribe{println it}

// Read Corection
process read_correction {
  tag { pair_id }
  
  publishDir "${output_dir}",
    mode: 'copy',
    pattern: "corrected_fastqs/*.fastq.gz"

  input:
  set pair_id, file(file_pair), genome_size from trimmed_fastqs_and_genome_size

  output:
  set pair_id, file('lighter.out') into read_correction_output
  set pair_id, file('corrected_fastqs/*.fastq.gz') into corrected_fastqs

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
  total_bases = m[0][1].toInteger() * 2 // the *2 is an estimate since number of reads >q25 in R2 may not be the same
  return [pair_id, total_bases]
}
base_counts = seqtk_fqchk_output.map { pair_id, file -> find_total_number_of_bases(pair_id, file.text) }
corrected_fastqs_and_genome_size_and_base_count = corrected_fastqs.join(genome_size_estimation_for_downsampling).join(base_counts).map{ tuple -> [tuple[0], tuple[1], tuple[2], tuple[3]]}

// merge reads and potentially downsample
process merge_reads{
  tag { pair_id }

  publishDir "${output_dir}",
    mode: 'copy',
    pattern: "merged_fastqs/*.fastq.gz"
  
  input:
  set pair_id, file(file_pair), genome_size, base_count from corrected_fastqs_and_genome_size_and_base_count

  output:
  set pair_id, file('merged_fastqs/*.fastq.gz') into merged_fastqs

  script:
  
  if (depth_cutoff  && base_count/genome_size > depth_cutoff.toInteger()){
    downsampling_factor = depth_cutoff.toInteger()/(base_count/genome_size)
    """
    mkdir downsampled_fastqs
    seqtk sample  ${file_pair[0]} ${downsampling_factor} | gzip > downsampled_fastqs/${file_pair[0]}
    seqtk sample  ${file_pair[1]} ${downsampling_factor} | gzip > downsampled_fastqs/${file_pair[1]}
    flash -m 20 -M 100 -d merged_fastqs -o ${pair_id} -z downsampled_fastqs/${file_pair[0]} downsampled_fastqs/${file_pair[1]} 
    """
  } else {
    """
    flash -m 20 -M 100 -d merged_fastqs -o ${pair_id} -z ${file_pair[0]} ${file_pair[1]} 
    """
  }

}

// assemble reads
process spades_assembly {
  memory '4 GB'
  
  tag { pair_id }

  input:
  set pair_id, file(file_triplet) from merged_fastqs

  output:
  set pair_id, file("scaffolds.fasta") into scaffolds

  """
  spades.py --pe1-1 ${file_triplet[1]} --pe1-2 ${file_triplet[2]} --pe1-m ${file_triplet[0]} --only-assembler  -o . --tmp-dir /tmp/${pair_id}_assembly -k 21,33,43,53,63,75 --threads 1 --memory 4
  """

}

process filter_scaffolds {
  tag { pair_id }

  publishDir "${output_dir}/assembly",
  mode: 'copy'

  input:
  set pair_id, file(scaffold_file) from scaffolds

  output:
  set pair_id, file("${pair_id}_scaffolds.fasta") into scaffolds_for_single_analysis
  file("${pair_id}_scaffolds.fasta") into scaffolds_for_combined_analysis
  
  """
  contig-tools filter -l ${minimum_scaffold_length} -c ${minimum_scaffold_depth} -f ${scaffold_file}
  ln -s scaffolds.filter_gt_${minimum_scaffold_length}bp_gt_${minimum_scaffold_depth}.0cov.fasta ${pair_id}_scaffolds.fasta
  """

}

// assess assembly with Quast
process quast {
  tag { pair_id }
    
  publishDir "${output_dir}/quast",
    mode: 'copy',
    pattern: "report.tsv",
    saveAs: { file -> "${pair_id}_quast_${file}"}

  input:
  set pair_id, file(contig_file) from scaffolds_for_single_analysis

  output:
  set pair_id, file('report.tsv')

  """
  quast.py ${contig_file} -o .
  """
}

// assess assembly with Quast but in a single file
process quast_summary {
  tag { 'quast_summary' }
  
  publishDir "${output_dir}/quast",
    mode: 'copy',
    pattern: "report.tsv",
    saveAs: { file -> "combined_quast_report.tsv"}

  input:
  file(contig_files) from scaffolds_for_combined_analysis.toSortedList{ it.getBaseName() }

  output:
  file('report.tsv')

  """
  quast.py ${contig_files} -o .
  """
}



workflow.onComplete {
  Helper.complete_message(params, workflow, version)
}

workflow.onError {
  Helper.error_message(workflow)
}

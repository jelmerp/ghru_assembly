#!/usr/bin/env nextflow
// Pipeline version
version = '1.3'
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
      --minimum_scaffold_length The minimum length of a scaffold to keep. Others will be filtered out. Default 500
      --minimum_scaffold_depth The minimum depth of coverage a scaffold must have to be kept. Others will be filtered out. Default 3
      --confindr_db_path The path to the confindr database. If not set assumes using Docker image where the path is '/home/bio/software_data/confindr_database'
      --qc_conditions Path to a YAML file containing pass/warning/fail conditions used by QualiFyr (https://gitlab.com/cgps/qualifyr)
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
params.confindr_db_path = false
params.qc_conditions = false

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

// set up confindr database path
if ( params.confindr_db_path ) {
  confindr_db_path = params.confindr_db_path
} else {
  // path in Docker image
  confindr_db_path = "/home/bio/software_data/confindr_database"
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
  set pair_id, stdout into min_read_length_for_trimming, min_read_length_for_assembly


  """
  gzip -cd ${file_pair[0]} | head -n 400000 | printf "%.0f" \$(awk 'NR%4==2{sum+=length(\$0)}END{print sum/(NR/4)/3}')
  """
}

// Pre-Trimming QC
process qc_pre_trimming {
  tag { pair_id }
  
  publishDir "${output_dir}/fastqc/pre_trimming",
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

min_read_length_and_raw_fastqs = min_read_length_for_trimming.join(raw_fastqs_for_trimming)

// Trimming
process trimming {
  memory '4 GB'
  
  tag { pair_id }
  
  input:
  set pair_id, min_read_length, file(file_pair) from min_read_length_and_raw_fastqs
  file('adapter_file.fas') from adapter_file

  output:
  set pair_id, file('trimmed_fastqs/*.fastq.gz') into trimmed_fastqs_for_qc, trimmed_fastqs_for_correction, trimmed_fastqs_for_genome_size_estimation, trimmed_fastqs_for_base_counting

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
  set file("${r1_prefix}_fastqc"), file("${r2_prefix}_fastqc") into fastqc_directories

  script:
  r1_prefix = file_pair[0].baseName.split('\\.')[0]
  r2_prefix = file_pair[1].baseName.split('\\.')[0]
  """
  fastqc ${file_pair[0]} ${file_pair[1]} --extract
  mv ${r1_prefix}_fastqc/summary.txt ${pair_id}_R1_fastqc.txt
  mv ${r2_prefix}_fastqc/summary.txt ${pair_id}_R2_fastqc.txt
  """
}

//FastQC MultiQC
process fastqc_multiqc {
  tag { 'multiqc for fastqc' }

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
  multiqc .
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
    confindr.py -i . -o . -d ${confindr_db_path} -t 1 -bf 0.025 -b 2 -Xmx 1500m
    """
  } else { // files with _1 and _2
    """
    confindr.py -i . -o . -d ${confindr_db_path} -t 1 -bf 0.025 -b 2 -Xmx 1500m -fid _1 -rid _2
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
  total_bases = m[0][1].toInteger() * 2 // the *2 is an estimate since number of reads >q25 in R2 may not be the same
  return [pair_id, total_bases]
}
base_counts = seqtk_fqchk_output.map { pair_id, file -> find_total_number_of_bases(pair_id, file.text) }
corrected_fastqs_and_genome_size_and_base_count = corrected_fastqs_for_merging.join(genome_size_estimation_for_downsampling).join(base_counts).map{ tuple -> [tuple[0], tuple[1], tuple[2], tuple[3]]}

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

 min_read_length_and_raw_fastqs = min_read_length_for_assembly.join(merged_fastqs)

// assemble reads
process spades_assembly {
  memory '4 GB'
  
  tag { pair_id }

  input:
  set pair_id, min_read_length, file(file_triplet) from min_read_length_and_raw_fastqs


  output:
  set pair_id, file("scaffolds.fasta") into scaffolds

  script:
  if (min_read_length.toInteger() < 27 ) {
    kmers = '21,33,43,53'
  } else {
    kmers = '21,33,43,53,63,75'
  }

  """
  spades.py --pe1-1 ${file_triplet[1]} --pe1-2 ${file_triplet[2]} --pe1-m ${file_triplet[0]} --only-assembler  -o . --tmp-dir /tmp/${pair_id}_assembly -k ${kmers} --threads 1 --memory 4

  """

}

// filter scaffolds to remove small and low coverage contigs
process filter_scaffolds {
  tag { pair_id }

  input:
  set pair_id, file(scaffold_file) from scaffolds

  output:
  set pair_id, file("${pair_id}_scaffolds.fasta") into scaffolds_for_single_analysis, scaffolds_for_qc
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
  memory '4 GB'
  
  publishDir "${output_dir}/quast",
    mode: 'copy',
    pattern: "*report.tsv",
    saveAs: { file -> "combined_${file}"}

  input:
  file(contig_files) from scaffolds_for_combined_analysis.toSortedList{ it.getBaseName() }

  output:
  file("*report.tsv")

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
  multiqc .
  """

}

// determine overall quality of sample
if (params.qc_conditions) {


  qc_conditions_yml = file(params.qc_conditions)
  quality_files = qc_post_trimming_files.join(confindr_files).join(quast_files_for_qualifyr).join(scaffolds_for_qc)
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
    set pair_id, file(fastqc_report_r1), file(fastqc_report_r2), file(confindr_report), file(quast_report), file(scaffold_file) from quality_files

    output:
    file('assemblies/**/*')
    file("${pair_id}.qualifyr.json") into qualifyr_json_files


    """
    result=`qualifyr check -y ${qc_conditions_yml} -f ${fastqc_report_r1} ${fastqc_report_r2} -c ${confindr_report}  -q ${quast_report} -s ${pair_id} 2> ERR`
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
    qualifyr check -y ${qc_conditions_yml} -f ${fastqc_report_r1} ${fastqc_report_r2} -c ${confindr_report}  -q ${quast_report} -s ${pair_id} -j -o .
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

//QualiFyr report
process qualifyr_report {
  tag { 'qualifyr report' }

  publishDir "${output_dir}/quality_reports",
    mode: 'copy',
    pattern: "qualifyr_report.html"

  input:
  
  file(json_files) from qualifyr_json_files.collect { it }

  output:
  file("qualifyr_report.html")

  script:
  """
  qualifyr report -i . -c 'quast.N50,quast.# contigs (>= 1000 bp),confindr.contam_status'
  """

}

workflow.onComplete {
  Helper.complete_message(params, workflow, version)
}

workflow.onError {
  Helper.error_message(workflow)
}

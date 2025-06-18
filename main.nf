#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ========= PARAMETERS =========
params.ref_fasta = null
params.reads     = '*.fastq.gz'

// ========= CHECK PARAMETERS =========
if ( !params.ref_fasta ) {
    exit 1, "You must specify the reference FASTA: --ref_fasta ref.fa"
}

workflow {
    // Create channel from reads glob pattern
    reads_ch = Channel.fromPath(params.reads, checkIfExists: true)
    
    // Create channel for reference FASTA
    ref_ch = Channel.fromPath(params.ref_fasta, checkIfExists: true)
    
    // Step 1: Prepare all reads
    all_reads = PREPARE_READS(reads_ch.collect())
    
    // Step 2: Build reference index and align
    (aligned_sam, alignment_stats, target_reads) = ALIGN_AND_EXTRACT(all_reads, ref_ch)
    
    // Step 3: Miniasm assembly
    (assembly_fasta, assembly_gfa, overlaps_paf) = ASSEMBLE(target_reads)
    
    // Step 4: Polish - round 1
    polished_round1 = POLISH1(assembly_fasta, target_reads)
    
    // Step 5: Polish - round 2
    polished_round2 = POLISH2(polished_round1, target_reads)
    
    // Step 6: Validate final assembly and publish results
    FINAL_VALIDATION(polished_round2, ref_ch)
}

// ========= PROCESSES =========

// Step 1: Concatenate reads
process PREPARE_READS {
    input:
    path read_files
    
    output:
    path "all_reads.fastq"
    
    script:
    """
    zcat ${read_files} > all_reads.fastq
    """
}

// Step 2: Alignment and extract aligned reads
process ALIGN_AND_EXTRACT {
    cpus 8
    
    input:
    path all_reads
    path ref_fasta
    
    output:
    path "aligned.sam"
    path "alignment_stats.txt"
    path "target_reads.fastq"
    
    script:
    """
    echo "Files in work directory:"
    ls -la
    
    echo "Using reference file: ${ref_fasta}"
    cp "${ref_fasta}" reference.fa
    
    if command -v minimap2 &> /dev/null; then
        echo "Using minimap2 for long read alignment..."
        minimap2 -ax map-ont -t ${task.cpus} reference.fa all_reads.fastq > aligned.sam 2> alignment_stats.txt
    else
        echo "minimap2 not found, using bowtie2..."
        bowtie2-build reference.fa reference
        bowtie2 --local --very-sensitive-local -p ${task.cpus} -x reference -U all_reads.fastq -S aligned.sam 2> alignment_stats.txt
    fi
    
    echo "Extracting aligned reads..."
    samtools view -h -F 4 aligned.sam | samtools fastq - > target_reads.fastq
    
    echo "Alignment complete!"
    """
}

// Step 3: Miniasm assembly
process ASSEMBLE {
    cpus 4
    
    input:
    path target_reads
    
    output:
    path "assembly.fasta"
    path "assembly.gfa"
    path "overlaps.paf"
    
    script:
    """
    echo "Starting overlap detection..."
    minimap2 -x ava-ont -t ${task.cpus} target_reads.fastq target_reads.fastq > overlaps.paf
    
    echo "Starting assembly..."
    miniasm -f target_reads.fastq overlaps.paf > assembly.gfa
    
    echo "Converting GFA to FASTA..."
    awk '/^S/{print ">"\$2; print \$3}' assembly.gfa > assembly.fasta
    
    echo "Assembly complete!"
    """
}

// Step 4: First round of polishing with Racon
process POLISH1 {
    cpus 4
    
    input:
    path assembly
    path target_reads
    
    output:
    path "polished_round1.fasta"
    
    script:
    """
    echo "Starting first polishing round..."
    minimap2 -ax map-ont -t ${task.cpus} ${assembly} target_reads.fastq > polish.sam
    racon -m 8 -x -6 -g -8 -w 500 -t ${task.cpus} target_reads.fastq polish.sam ${assembly} > polished_round1.fasta
    echo "First polishing round complete!"
    """
}

// Step 5: Second round of polishing with Racon
process POLISH2 {
    cpus 4
    
    input:
    path polished
    path target_reads
    
    output:
    path "polished_round2.fasta"
    
    script:
    """
    echo "Starting second polishing round..."
    minimap2 -ax map-ont -t ${task.cpus} ${polished} target_reads.fastq > polish2.sam
    racon -m 8 -x -6 -g -8 -w 500 -t ${task.cpus} target_reads.fastq polish2.sam ${polished} > polished_round2.fasta
    echo "Second polishing round complete!"
    """
}

// Step 6: Final validation
process FINAL_VALIDATION {
    publishDir "./results", mode: 'copy'
    cpus 4
    
    input:
    path polished
    path ref_fasta
    
    output:
    path "final_vs_ref.sam"
    path "final_coverage.txt"
    path "final_assembly.fasta"
    
    script:
    """
    echo "Using reference file: ${ref_fasta}"
    cp "${ref_fasta}" reference.fa
    
    echo "Aligning final assembly to reference..."
    minimap2 -ax asm5 -t ${task.cpus} reference.fa ${polished} > final_vs_ref.sam
    
    echo "Calculating coverage..."
    samtools view -bS final_vs_ref.sam | samtools sort -@ ${task.cpus} - | samtools depth - > final_coverage.txt
    
    cp ${polished} final_assembly.fasta
    
    echo "Final validation complete!"
    """
}

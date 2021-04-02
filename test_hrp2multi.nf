
params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads = params.input.fastq_path+'/*_{0,1,2}.fastq'
fastq_path = Channel.fromFilePairs(params.reads,size:3)

params.genome = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt {
    publishDir "$params.output.folder/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
        set val(id),  path("${id}_trimmed_0.fq") into trim_out1
        set val(id),  path("${id}_trimmed_1.fq") into trim_out2
        set val(id),  path("${id}_trimmed_2.fq") into trim_out3

    script:
        """
        
        NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
        NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
        NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq

        """



}

process maptoreference {
    publishDir "$params.output.folder/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
        set val(id), path(trim_read0) from trim_out1
        set val(id), path(trim_read1) from trim_out2
        set val(id), path(trim_read2) from trim_out3
        path(ref_dir) from ref_dir

       
    output:
        set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1
        set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2
        set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3

    script:
        """
        bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
        bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta

        """

}


process assemble {

    publishDir "$params.output.folder/contigs/${sample}", mode : "copy"
    input:
        set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1
        set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2
        set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
        """
        $baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
       
        
        """
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate {

    publishDir "$params.output.folder/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern {

    publishDir "$params.output.folder/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}



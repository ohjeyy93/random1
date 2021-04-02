
params.reads1 = params.input.fastq_path1+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads1,size:1)
params.genome1 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt1 {
    publishDir "$params.output.folder1/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out1

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		"""
        




}

process maptoreference1 {
    publishDir "$params.output.folder1/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder1/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out1

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble1 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder1/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate1 {

    publishDir "$params.output.folder1/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern1 {

    publishDir "$params.output.folder1/resultfiles/${sample}", mode : "copy"
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


params.reads2 = params.input.fastq_path2+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads2,size:1)
params.genome2 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt2 {
    publishDir "$params.output.folder2/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out1024

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		"""
        




}

process maptoreference2 {
    publishDir "$params.output.folder2/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder2/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out1024

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1024

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble2 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder2/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1024


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1024
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate2 {

    publishDir "$params.output.folder2/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1024
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1024

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern2 {

    publishDir "$params.output.folder2/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1024
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads3 = params.input.fastq_path3+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads3,size:1)
params.genome3 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt3 {
    publishDir "$params.output.folder3/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out59049

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		"""
        




}

process maptoreference3 {
    publishDir "$params.output.folder3/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder3/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out59049

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out59049

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble3 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder3/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out59049


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout59049
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate3 {

    publishDir "$params.output.folder3/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout59049
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq59049

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern3 {

    publishDir "$params.output.folder3/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq59049
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads4 = params.input.fastq_path4+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads4,size:17)
params.genome4 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt4 {
    publishDir "$params.output.folder4/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out1048576
		set val(id),  path("${id}_trimmed_1.fq") into trim_out2097152
		set val(id),  path("${id}_trimmed_2.fq") into trim_out3145728
		set val(id),  path("${id}_trimmed_3.fq") into trim_out4194304
		set val(id),  path("${id}_trimmed_4.fq") into trim_out5242880
		set val(id),  path("${id}_trimmed_5.fq") into trim_out6291456
		set val(id),  path("${id}_trimmed_6.fq") into trim_out7340032
		set val(id),  path("${id}_trimmed_7.fq") into trim_out8388608
		set val(id),  path("${id}_trimmed_8.fq") into trim_out9437184
		set val(id),  path("${id}_trimmed_9.fq") into trim_out10485760
		set val(id),  path("${id}_trimmed_10.fq") into trim_out11534336
		set val(id),  path("${id}_trimmed_11.fq") into trim_out12582912
		set val(id),  path("${id}_trimmed_12.fq") into trim_out13631488
		set val(id),  path("${id}_trimmed_13.fq") into trim_out14680064
		set val(id),  path("${id}_trimmed_14.fq") into trim_out15728640
		set val(id),  path("${id}_trimmed_15.fq") into trim_out16777216
		set val(id),  path("${id}_trimmed_16.fq") into trim_out17825792

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
		NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq
		NanoFilt ${id}_3.fastq -l 500 -q 10 > ${id}_trimmed_3.fq
		NanoFilt ${id}_4.fastq -l 500 -q 10 > ${id}_trimmed_4.fq
		NanoFilt ${id}_5.fastq -l 500 -q 10 > ${id}_trimmed_5.fq
		NanoFilt ${id}_6.fastq -l 500 -q 10 > ${id}_trimmed_6.fq
		NanoFilt ${id}_7.fastq -l 500 -q 10 > ${id}_trimmed_7.fq
		NanoFilt ${id}_8.fastq -l 500 -q 10 > ${id}_trimmed_8.fq
		NanoFilt ${id}_9.fastq -l 500 -q 10 > ${id}_trimmed_9.fq
		NanoFilt ${id}_10.fastq -l 500 -q 10 > ${id}_trimmed_10.fq
		NanoFilt ${id}_11.fastq -l 500 -q 10 > ${id}_trimmed_11.fq
		NanoFilt ${id}_12.fastq -l 500 -q 10 > ${id}_trimmed_12.fq
		NanoFilt ${id}_13.fastq -l 500 -q 10 > ${id}_trimmed_13.fq
		NanoFilt ${id}_14.fastq -l 500 -q 10 > ${id}_trimmed_14.fq
		NanoFilt ${id}_15.fastq -l 500 -q 10 > ${id}_trimmed_15.fq
		NanoFilt ${id}_16.fastq -l 500 -q 10 > ${id}_trimmed_16.fq
		"""
        




}

process maptoreference4 {
    publishDir "$params.output.folder4/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder4/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out1048576
		set val(id), path(trim_read1) from trim_out2097152
		set val(id), path(trim_read2) from trim_out3145728
		set val(id), path(trim_read3) from trim_out4194304
		set val(id), path(trim_read4) from trim_out5242880
		set val(id), path(trim_read5) from trim_out6291456
		set val(id), path(trim_read6) from trim_out7340032
		set val(id), path(trim_read7) from trim_out8388608
		set val(id), path(trim_read8) from trim_out9437184
		set val(id), path(trim_read9) from trim_out10485760
		set val(id), path(trim_read10) from trim_out11534336
		set val(id), path(trim_read11) from trim_out12582912
		set val(id), path(trim_read12) from trim_out13631488
		set val(id), path(trim_read13) from trim_out14680064
		set val(id), path(trim_read14) from trim_out15728640
		set val(id), path(trim_read15) from trim_out16777216
		set val(id), path(trim_read16) from trim_out17825792

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1048576
		set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2097152
		set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3145728
		set val(id), path("${id}_mapped_3.fq"), path("${id}_unmapped_3.fq") into mapped_out4194304
		set val(id), path("${id}_mapped_4.fq"), path("${id}_unmapped_4.fq") into mapped_out5242880
		set val(id), path("${id}_mapped_5.fq"), path("${id}_unmapped_5.fq") into mapped_out6291456
		set val(id), path("${id}_mapped_6.fq"), path("${id}_unmapped_6.fq") into mapped_out7340032
		set val(id), path("${id}_mapped_7.fq"), path("${id}_unmapped_7.fq") into mapped_out8388608
		set val(id), path("${id}_mapped_8.fq"), path("${id}_unmapped_8.fq") into mapped_out9437184
		set val(id), path("${id}_mapped_9.fq"), path("${id}_unmapped_9.fq") into mapped_out10485760
		set val(id), path("${id}_mapped_10.fq"), path("${id}_unmapped_10.fq") into mapped_out11534336
		set val(id), path("${id}_mapped_11.fq"), path("${id}_unmapped_11.fq") into mapped_out12582912
		set val(id), path("${id}_mapped_12.fq"), path("${id}_unmapped_12.fq") into mapped_out13631488
		set val(id), path("${id}_mapped_13.fq"), path("${id}_unmapped_13.fq") into mapped_out14680064
		set val(id), path("${id}_mapped_14.fq"), path("${id}_unmapped_14.fq") into mapped_out15728640
		set val(id), path("${id}_mapped_15.fq"), path("${id}_unmapped_15.fq") into mapped_out16777216
		set val(id), path("${id}_mapped_16.fq"), path("${id}_unmapped_16.fq") into mapped_out17825792

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_3.fq outm=${id}_mapped_3.fq outu=${id}_unmapped_3.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_4.fq outm=${id}_mapped_4.fq outu=${id}_unmapped_4.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_5.fq outm=${id}_mapped_5.fq outu=${id}_unmapped_5.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_6.fq outm=${id}_mapped_6.fq outu=${id}_unmapped_6.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_7.fq outm=${id}_mapped_7.fq outu=${id}_unmapped_7.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_8.fq outm=${id}_mapped_8.fq outu=${id}_unmapped_8.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_9.fq outm=${id}_mapped_9.fq outu=${id}_unmapped_9.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_10.fq outm=${id}_mapped_10.fq outu=${id}_unmapped_10.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_11.fq outm=${id}_mapped_11.fq outu=${id}_unmapped_11.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_12.fq outm=${id}_mapped_12.fq outu=${id}_unmapped_12.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_13.fq outm=${id}_mapped_13.fq outu=${id}_unmapped_13.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_14.fq outm=${id}_mapped_14.fq outu=${id}_unmapped_14.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_15.fq outm=${id}_mapped_15.fq outu=${id}_unmapped_15.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_16.fq outm=${id}_mapped_16.fq outu=${id}_unmapped_16.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble4 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder4/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1048576
		set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2097152
		set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3145728
		set val(sample), path(mapped_read3), path(unmapped_read3) from mapped_out4194304
		set val(sample), path(mapped_read4), path(unmapped_read4) from mapped_out5242880
		set val(sample), path(mapped_read5), path(unmapped_read5) from mapped_out6291456
		set val(sample), path(mapped_read6), path(unmapped_read6) from mapped_out7340032
		set val(sample), path(mapped_read7), path(unmapped_read7) from mapped_out8388608
		set val(sample), path(mapped_read8), path(unmapped_read8) from mapped_out9437184
		set val(sample), path(mapped_read9), path(unmapped_read9) from mapped_out10485760
		set val(sample), path(mapped_read10), path(unmapped_read10) from mapped_out11534336
		set val(sample), path(mapped_read11), path(unmapped_read11) from mapped_out12582912
		set val(sample), path(mapped_read12), path(unmapped_read12) from mapped_out13631488
		set val(sample), path(mapped_read13), path(unmapped_read13) from mapped_out14680064
		set val(sample), path(mapped_read14), path(unmapped_read14) from mapped_out15728640
		set val(sample), path(mapped_read15), path(unmapped_read15) from mapped_out16777216
		set val(sample), path(mapped_read16), path(unmapped_read16) from mapped_out17825792


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1048576
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq  -s ${sample}_mapped_1.fq  -s ${sample}_mapped_2.fq  -s ${sample}_mapped_3.fq  -s ${sample}_mapped_4.fq  -s ${sample}_mapped_5.fq  -s ${sample}_mapped_6.fq  -s ${sample}_mapped_7.fq  -s ${sample}_mapped_8.fq  -s ${sample}_mapped_9.fq  -s ${sample}_mapped_10.fq  -s ${sample}_mapped_11.fq  -s ${sample}_mapped_12.fq  -s ${sample}_mapped_13.fq  -s ${sample}_mapped_14.fq  -s ${sample}_mapped_15.fq  -s ${sample}_mapped_16.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate4 {

    publishDir "$params.output.folder4/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1048576
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1048576

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern4 {

    publishDir "$params.output.folder4/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1048576
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads5 = params.input.fastq_path5+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads5,size:13)
params.genome5 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt5 {
    publishDir "$params.output.folder5/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out9765625
		set val(id),  path("${id}_trimmed_1.fq") into trim_out19531250
		set val(id),  path("${id}_trimmed_2.fq") into trim_out29296875
		set val(id),  path("${id}_trimmed_3.fq") into trim_out39062500
		set val(id),  path("${id}_trimmed_4.fq") into trim_out48828125
		set val(id),  path("${id}_trimmed_5.fq") into trim_out58593750
		set val(id),  path("${id}_trimmed_6.fq") into trim_out68359375
		set val(id),  path("${id}_trimmed_7.fq") into trim_out78125000
		set val(id),  path("${id}_trimmed_8.fq") into trim_out87890625
		set val(id),  path("${id}_trimmed_9.fq") into trim_out97656250
		set val(id),  path("${id}_trimmed_10.fq") into trim_out107421875
		set val(id),  path("${id}_trimmed_11.fq") into trim_out117187500
		set val(id),  path("${id}_trimmed_12.fq") into trim_out126953125

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
		NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq
		NanoFilt ${id}_3.fastq -l 500 -q 10 > ${id}_trimmed_3.fq
		NanoFilt ${id}_4.fastq -l 500 -q 10 > ${id}_trimmed_4.fq
		NanoFilt ${id}_5.fastq -l 500 -q 10 > ${id}_trimmed_5.fq
		NanoFilt ${id}_6.fastq -l 500 -q 10 > ${id}_trimmed_6.fq
		NanoFilt ${id}_7.fastq -l 500 -q 10 > ${id}_trimmed_7.fq
		NanoFilt ${id}_8.fastq -l 500 -q 10 > ${id}_trimmed_8.fq
		NanoFilt ${id}_9.fastq -l 500 -q 10 > ${id}_trimmed_9.fq
		NanoFilt ${id}_10.fastq -l 500 -q 10 > ${id}_trimmed_10.fq
		NanoFilt ${id}_11.fastq -l 500 -q 10 > ${id}_trimmed_11.fq
		NanoFilt ${id}_12.fastq -l 500 -q 10 > ${id}_trimmed_12.fq
		"""
        




}

process maptoreference5 {
    publishDir "$params.output.folder5/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder5/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out9765625
		set val(id), path(trim_read1) from trim_out19531250
		set val(id), path(trim_read2) from trim_out29296875
		set val(id), path(trim_read3) from trim_out39062500
		set val(id), path(trim_read4) from trim_out48828125
		set val(id), path(trim_read5) from trim_out58593750
		set val(id), path(trim_read6) from trim_out68359375
		set val(id), path(trim_read7) from trim_out78125000
		set val(id), path(trim_read8) from trim_out87890625
		set val(id), path(trim_read9) from trim_out97656250
		set val(id), path(trim_read10) from trim_out107421875
		set val(id), path(trim_read11) from trim_out117187500
		set val(id), path(trim_read12) from trim_out126953125

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out9765625
		set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out19531250
		set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out29296875
		set val(id), path("${id}_mapped_3.fq"), path("${id}_unmapped_3.fq") into mapped_out39062500
		set val(id), path("${id}_mapped_4.fq"), path("${id}_unmapped_4.fq") into mapped_out48828125
		set val(id), path("${id}_mapped_5.fq"), path("${id}_unmapped_5.fq") into mapped_out58593750
		set val(id), path("${id}_mapped_6.fq"), path("${id}_unmapped_6.fq") into mapped_out68359375
		set val(id), path("${id}_mapped_7.fq"), path("${id}_unmapped_7.fq") into mapped_out78125000
		set val(id), path("${id}_mapped_8.fq"), path("${id}_unmapped_8.fq") into mapped_out87890625
		set val(id), path("${id}_mapped_9.fq"), path("${id}_unmapped_9.fq") into mapped_out97656250
		set val(id), path("${id}_mapped_10.fq"), path("${id}_unmapped_10.fq") into mapped_out107421875
		set val(id), path("${id}_mapped_11.fq"), path("${id}_unmapped_11.fq") into mapped_out117187500
		set val(id), path("${id}_mapped_12.fq"), path("${id}_unmapped_12.fq") into mapped_out126953125

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_3.fq outm=${id}_mapped_3.fq outu=${id}_unmapped_3.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_4.fq outm=${id}_mapped_4.fq outu=${id}_unmapped_4.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_5.fq outm=${id}_mapped_5.fq outu=${id}_unmapped_5.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_6.fq outm=${id}_mapped_6.fq outu=${id}_unmapped_6.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_7.fq outm=${id}_mapped_7.fq outu=${id}_unmapped_7.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_8.fq outm=${id}_mapped_8.fq outu=${id}_unmapped_8.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_9.fq outm=${id}_mapped_9.fq outu=${id}_unmapped_9.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_10.fq outm=${id}_mapped_10.fq outu=${id}_unmapped_10.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_11.fq outm=${id}_mapped_11.fq outu=${id}_unmapped_11.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_12.fq outm=${id}_mapped_12.fq outu=${id}_unmapped_12.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble5 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder5/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out9765625
		set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out19531250
		set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out29296875
		set val(sample), path(mapped_read3), path(unmapped_read3) from mapped_out39062500
		set val(sample), path(mapped_read4), path(unmapped_read4) from mapped_out48828125
		set val(sample), path(mapped_read5), path(unmapped_read5) from mapped_out58593750
		set val(sample), path(mapped_read6), path(unmapped_read6) from mapped_out68359375
		set val(sample), path(mapped_read7), path(unmapped_read7) from mapped_out78125000
		set val(sample), path(mapped_read8), path(unmapped_read8) from mapped_out87890625
		set val(sample), path(mapped_read9), path(unmapped_read9) from mapped_out97656250
		set val(sample), path(mapped_read10), path(unmapped_read10) from mapped_out107421875
		set val(sample), path(mapped_read11), path(unmapped_read11) from mapped_out117187500
		set val(sample), path(mapped_read12), path(unmapped_read12) from mapped_out126953125


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout9765625
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq  -s ${sample}_mapped_1.fq  -s ${sample}_mapped_2.fq  -s ${sample}_mapped_3.fq  -s ${sample}_mapped_4.fq  -s ${sample}_mapped_5.fq  -s ${sample}_mapped_6.fq  -s ${sample}_mapped_7.fq  -s ${sample}_mapped_8.fq  -s ${sample}_mapped_9.fq  -s ${sample}_mapped_10.fq  -s ${sample}_mapped_11.fq  -s ${sample}_mapped_12.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate5 {

    publishDir "$params.output.folder5/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout9765625
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq9765625

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern5 {

    publishDir "$params.output.folder5/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq9765625
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads6 = params.input.fastq_path6+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads6,size:14)
params.genome6 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt6 {
    publishDir "$params.output.folder6/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out60466176
		set val(id),  path("${id}_trimmed_1.fq") into trim_out120932352
		set val(id),  path("${id}_trimmed_2.fq") into trim_out181398528
		set val(id),  path("${id}_trimmed_3.fq") into trim_out241864704
		set val(id),  path("${id}_trimmed_4.fq") into trim_out302330880
		set val(id),  path("${id}_trimmed_5.fq") into trim_out362797056
		set val(id),  path("${id}_trimmed_6.fq") into trim_out423263232
		set val(id),  path("${id}_trimmed_7.fq") into trim_out483729408
		set val(id),  path("${id}_trimmed_8.fq") into trim_out544195584
		set val(id),  path("${id}_trimmed_9.fq") into trim_out604661760
		set val(id),  path("${id}_trimmed_10.fq") into trim_out665127936
		set val(id),  path("${id}_trimmed_11.fq") into trim_out725594112
		set val(id),  path("${id}_trimmed_12.fq") into trim_out786060288
		set val(id),  path("${id}_trimmed_13.fq") into trim_out846526464

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
		NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq
		NanoFilt ${id}_3.fastq -l 500 -q 10 > ${id}_trimmed_3.fq
		NanoFilt ${id}_4.fastq -l 500 -q 10 > ${id}_trimmed_4.fq
		NanoFilt ${id}_5.fastq -l 500 -q 10 > ${id}_trimmed_5.fq
		NanoFilt ${id}_6.fastq -l 500 -q 10 > ${id}_trimmed_6.fq
		NanoFilt ${id}_7.fastq -l 500 -q 10 > ${id}_trimmed_7.fq
		NanoFilt ${id}_8.fastq -l 500 -q 10 > ${id}_trimmed_8.fq
		NanoFilt ${id}_9.fastq -l 500 -q 10 > ${id}_trimmed_9.fq
		NanoFilt ${id}_10.fastq -l 500 -q 10 > ${id}_trimmed_10.fq
		NanoFilt ${id}_11.fastq -l 500 -q 10 > ${id}_trimmed_11.fq
		NanoFilt ${id}_12.fastq -l 500 -q 10 > ${id}_trimmed_12.fq
		NanoFilt ${id}_13.fastq -l 500 -q 10 > ${id}_trimmed_13.fq
		"""
        




}

process maptoreference6 {
    publishDir "$params.output.folder6/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder6/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out60466176
		set val(id), path(trim_read1) from trim_out120932352
		set val(id), path(trim_read2) from trim_out181398528
		set val(id), path(trim_read3) from trim_out241864704
		set val(id), path(trim_read4) from trim_out302330880
		set val(id), path(trim_read5) from trim_out362797056
		set val(id), path(trim_read6) from trim_out423263232
		set val(id), path(trim_read7) from trim_out483729408
		set val(id), path(trim_read8) from trim_out544195584
		set val(id), path(trim_read9) from trim_out604661760
		set val(id), path(trim_read10) from trim_out665127936
		set val(id), path(trim_read11) from trim_out725594112
		set val(id), path(trim_read12) from trim_out786060288
		set val(id), path(trim_read13) from trim_out846526464

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out60466176
		set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out120932352
		set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out181398528
		set val(id), path("${id}_mapped_3.fq"), path("${id}_unmapped_3.fq") into mapped_out241864704
		set val(id), path("${id}_mapped_4.fq"), path("${id}_unmapped_4.fq") into mapped_out302330880
		set val(id), path("${id}_mapped_5.fq"), path("${id}_unmapped_5.fq") into mapped_out362797056
		set val(id), path("${id}_mapped_6.fq"), path("${id}_unmapped_6.fq") into mapped_out423263232
		set val(id), path("${id}_mapped_7.fq"), path("${id}_unmapped_7.fq") into mapped_out483729408
		set val(id), path("${id}_mapped_8.fq"), path("${id}_unmapped_8.fq") into mapped_out544195584
		set val(id), path("${id}_mapped_9.fq"), path("${id}_unmapped_9.fq") into mapped_out604661760
		set val(id), path("${id}_mapped_10.fq"), path("${id}_unmapped_10.fq") into mapped_out665127936
		set val(id), path("${id}_mapped_11.fq"), path("${id}_unmapped_11.fq") into mapped_out725594112
		set val(id), path("${id}_mapped_12.fq"), path("${id}_unmapped_12.fq") into mapped_out786060288
		set val(id), path("${id}_mapped_13.fq"), path("${id}_unmapped_13.fq") into mapped_out846526464

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_3.fq outm=${id}_mapped_3.fq outu=${id}_unmapped_3.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_4.fq outm=${id}_mapped_4.fq outu=${id}_unmapped_4.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_5.fq outm=${id}_mapped_5.fq outu=${id}_unmapped_5.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_6.fq outm=${id}_mapped_6.fq outu=${id}_unmapped_6.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_7.fq outm=${id}_mapped_7.fq outu=${id}_unmapped_7.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_8.fq outm=${id}_mapped_8.fq outu=${id}_unmapped_8.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_9.fq outm=${id}_mapped_9.fq outu=${id}_unmapped_9.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_10.fq outm=${id}_mapped_10.fq outu=${id}_unmapped_10.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_11.fq outm=${id}_mapped_11.fq outu=${id}_unmapped_11.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_12.fq outm=${id}_mapped_12.fq outu=${id}_unmapped_12.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_13.fq outm=${id}_mapped_13.fq outu=${id}_unmapped_13.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble6 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder6/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out60466176
		set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out120932352
		set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out181398528
		set val(sample), path(mapped_read3), path(unmapped_read3) from mapped_out241864704
		set val(sample), path(mapped_read4), path(unmapped_read4) from mapped_out302330880
		set val(sample), path(mapped_read5), path(unmapped_read5) from mapped_out362797056
		set val(sample), path(mapped_read6), path(unmapped_read6) from mapped_out423263232
		set val(sample), path(mapped_read7), path(unmapped_read7) from mapped_out483729408
		set val(sample), path(mapped_read8), path(unmapped_read8) from mapped_out544195584
		set val(sample), path(mapped_read9), path(unmapped_read9) from mapped_out604661760
		set val(sample), path(mapped_read10), path(unmapped_read10) from mapped_out665127936
		set val(sample), path(mapped_read11), path(unmapped_read11) from mapped_out725594112
		set val(sample), path(mapped_read12), path(unmapped_read12) from mapped_out786060288
		set val(sample), path(mapped_read13), path(unmapped_read13) from mapped_out846526464


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout60466176
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq  -s ${sample}_mapped_1.fq  -s ${sample}_mapped_2.fq  -s ${sample}_mapped_3.fq  -s ${sample}_mapped_4.fq  -s ${sample}_mapped_5.fq  -s ${sample}_mapped_6.fq  -s ${sample}_mapped_7.fq  -s ${sample}_mapped_8.fq  -s ${sample}_mapped_9.fq  -s ${sample}_mapped_10.fq  -s ${sample}_mapped_11.fq  -s ${sample}_mapped_12.fq  -s ${sample}_mapped_13.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate6 {

    publishDir "$params.output.folder6/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout60466176
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq60466176

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern6 {

    publishDir "$params.output.folder6/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq60466176
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads7 = params.input.fastq_path7+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads7,size:9)
params.genome7 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt7 {
    publishDir "$params.output.folder7/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out282475249
		set val(id),  path("${id}_trimmed_1.fq") into trim_out564950498
		set val(id),  path("${id}_trimmed_2.fq") into trim_out847425747
		set val(id),  path("${id}_trimmed_3.fq") into trim_out1129900996
		set val(id),  path("${id}_trimmed_4.fq") into trim_out1412376245
		set val(id),  path("${id}_trimmed_5.fq") into trim_out1694851494
		set val(id),  path("${id}_trimmed_6.fq") into trim_out1977326743
		set val(id),  path("${id}_trimmed_7.fq") into trim_out2259801992
		set val(id),  path("${id}_trimmed_8.fq") into trim_out2542277241

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
		NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq
		NanoFilt ${id}_3.fastq -l 500 -q 10 > ${id}_trimmed_3.fq
		NanoFilt ${id}_4.fastq -l 500 -q 10 > ${id}_trimmed_4.fq
		NanoFilt ${id}_5.fastq -l 500 -q 10 > ${id}_trimmed_5.fq
		NanoFilt ${id}_6.fastq -l 500 -q 10 > ${id}_trimmed_6.fq
		NanoFilt ${id}_7.fastq -l 500 -q 10 > ${id}_trimmed_7.fq
		NanoFilt ${id}_8.fastq -l 500 -q 10 > ${id}_trimmed_8.fq
		"""
        




}

process maptoreference7 {
    publishDir "$params.output.folder7/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder7/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out282475249
		set val(id), path(trim_read1) from trim_out564950498
		set val(id), path(trim_read2) from trim_out847425747
		set val(id), path(trim_read3) from trim_out1129900996
		set val(id), path(trim_read4) from trim_out1412376245
		set val(id), path(trim_read5) from trim_out1694851494
		set val(id), path(trim_read6) from trim_out1977326743
		set val(id), path(trim_read7) from trim_out2259801992
		set val(id), path(trim_read8) from trim_out2542277241

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out282475249
		set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out564950498
		set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out847425747
		set val(id), path("${id}_mapped_3.fq"), path("${id}_unmapped_3.fq") into mapped_out1129900996
		set val(id), path("${id}_mapped_4.fq"), path("${id}_unmapped_4.fq") into mapped_out1412376245
		set val(id), path("${id}_mapped_5.fq"), path("${id}_unmapped_5.fq") into mapped_out1694851494
		set val(id), path("${id}_mapped_6.fq"), path("${id}_unmapped_6.fq") into mapped_out1977326743
		set val(id), path("${id}_mapped_7.fq"), path("${id}_unmapped_7.fq") into mapped_out2259801992
		set val(id), path("${id}_mapped_8.fq"), path("${id}_unmapped_8.fq") into mapped_out2542277241

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_3.fq outm=${id}_mapped_3.fq outu=${id}_unmapped_3.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_4.fq outm=${id}_mapped_4.fq outu=${id}_unmapped_4.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_5.fq outm=${id}_mapped_5.fq outu=${id}_unmapped_5.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_6.fq outm=${id}_mapped_6.fq outu=${id}_unmapped_6.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_7.fq outm=${id}_mapped_7.fq outu=${id}_unmapped_7.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_8.fq outm=${id}_mapped_8.fq outu=${id}_unmapped_8.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble7 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder7/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out282475249
		set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out564950498
		set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out847425747
		set val(sample), path(mapped_read3), path(unmapped_read3) from mapped_out1129900996
		set val(sample), path(mapped_read4), path(unmapped_read4) from mapped_out1412376245
		set val(sample), path(mapped_read5), path(unmapped_read5) from mapped_out1694851494
		set val(sample), path(mapped_read6), path(unmapped_read6) from mapped_out1977326743
		set val(sample), path(mapped_read7), path(unmapped_read7) from mapped_out2259801992
		set val(sample), path(mapped_read8), path(unmapped_read8) from mapped_out2542277241


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout282475249
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq  -s ${sample}_mapped_1.fq  -s ${sample}_mapped_2.fq  -s ${sample}_mapped_3.fq  -s ${sample}_mapped_4.fq  -s ${sample}_mapped_5.fq  -s ${sample}_mapped_6.fq  -s ${sample}_mapped_7.fq  -s ${sample}_mapped_8.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate7 {

    publishDir "$params.output.folder7/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout282475249
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq282475249

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern7 {

    publishDir "$params.output.folder7/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq282475249
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads8 = params.input.fastq_path8+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads8,size:78)
params.genome8 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt8 {
    publishDir "$params.output.folder8/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out1073741824
		set val(id),  path("${id}_trimmed_1.fq") into trim_out2147483648
		set val(id),  path("${id}_trimmed_2.fq") into trim_out3221225472
		set val(id),  path("${id}_trimmed_3.fq") into trim_out4294967296
		set val(id),  path("${id}_trimmed_4.fq") into trim_out5368709120
		set val(id),  path("${id}_trimmed_5.fq") into trim_out6442450944
		set val(id),  path("${id}_trimmed_6.fq") into trim_out7516192768
		set val(id),  path("${id}_trimmed_7.fq") into trim_out8589934592
		set val(id),  path("${id}_trimmed_8.fq") into trim_out9663676416
		set val(id),  path("${id}_trimmed_9.fq") into trim_out10737418240
		set val(id),  path("${id}_trimmed_10.fq") into trim_out11811160064
		set val(id),  path("${id}_trimmed_11.fq") into trim_out12884901888
		set val(id),  path("${id}_trimmed_12.fq") into trim_out13958643712
		set val(id),  path("${id}_trimmed_13.fq") into trim_out15032385536
		set val(id),  path("${id}_trimmed_14.fq") into trim_out16106127360
		set val(id),  path("${id}_trimmed_15.fq") into trim_out17179869184
		set val(id),  path("${id}_trimmed_16.fq") into trim_out18253611008
		set val(id),  path("${id}_trimmed_17.fq") into trim_out19327352832
		set val(id),  path("${id}_trimmed_18.fq") into trim_out20401094656
		set val(id),  path("${id}_trimmed_19.fq") into trim_out21474836480
		set val(id),  path("${id}_trimmed_20.fq") into trim_out22548578304
		set val(id),  path("${id}_trimmed_21.fq") into trim_out23622320128
		set val(id),  path("${id}_trimmed_22.fq") into trim_out24696061952
		set val(id),  path("${id}_trimmed_23.fq") into trim_out25769803776
		set val(id),  path("${id}_trimmed_24.fq") into trim_out26843545600
		set val(id),  path("${id}_trimmed_25.fq") into trim_out27917287424
		set val(id),  path("${id}_trimmed_26.fq") into trim_out28991029248
		set val(id),  path("${id}_trimmed_27.fq") into trim_out30064771072
		set val(id),  path("${id}_trimmed_28.fq") into trim_out31138512896
		set val(id),  path("${id}_trimmed_29.fq") into trim_out32212254720
		set val(id),  path("${id}_trimmed_30.fq") into trim_out33285996544
		set val(id),  path("${id}_trimmed_31.fq") into trim_out34359738368
		set val(id),  path("${id}_trimmed_32.fq") into trim_out35433480192
		set val(id),  path("${id}_trimmed_33.fq") into trim_out36507222016
		set val(id),  path("${id}_trimmed_34.fq") into trim_out37580963840
		set val(id),  path("${id}_trimmed_35.fq") into trim_out38654705664
		set val(id),  path("${id}_trimmed_36.fq") into trim_out39728447488
		set val(id),  path("${id}_trimmed_37.fq") into trim_out40802189312
		set val(id),  path("${id}_trimmed_38.fq") into trim_out41875931136
		set val(id),  path("${id}_trimmed_39.fq") into trim_out42949672960
		set val(id),  path("${id}_trimmed_40.fq") into trim_out44023414784
		set val(id),  path("${id}_trimmed_41.fq") into trim_out45097156608
		set val(id),  path("${id}_trimmed_42.fq") into trim_out46170898432
		set val(id),  path("${id}_trimmed_43.fq") into trim_out47244640256
		set val(id),  path("${id}_trimmed_44.fq") into trim_out48318382080
		set val(id),  path("${id}_trimmed_45.fq") into trim_out49392123904
		set val(id),  path("${id}_trimmed_46.fq") into trim_out50465865728
		set val(id),  path("${id}_trimmed_47.fq") into trim_out51539607552
		set val(id),  path("${id}_trimmed_48.fq") into trim_out52613349376
		set val(id),  path("${id}_trimmed_49.fq") into trim_out53687091200
		set val(id),  path("${id}_trimmed_50.fq") into trim_out54760833024
		set val(id),  path("${id}_trimmed_51.fq") into trim_out55834574848
		set val(id),  path("${id}_trimmed_52.fq") into trim_out56908316672
		set val(id),  path("${id}_trimmed_53.fq") into trim_out57982058496
		set val(id),  path("${id}_trimmed_54.fq") into trim_out59055800320
		set val(id),  path("${id}_trimmed_55.fq") into trim_out60129542144
		set val(id),  path("${id}_trimmed_56.fq") into trim_out61203283968
		set val(id),  path("${id}_trimmed_57.fq") into trim_out62277025792
		set val(id),  path("${id}_trimmed_58.fq") into trim_out63350767616
		set val(id),  path("${id}_trimmed_59.fq") into trim_out64424509440
		set val(id),  path("${id}_trimmed_60.fq") into trim_out65498251264
		set val(id),  path("${id}_trimmed_61.fq") into trim_out66571993088
		set val(id),  path("${id}_trimmed_62.fq") into trim_out67645734912
		set val(id),  path("${id}_trimmed_63.fq") into trim_out68719476736
		set val(id),  path("${id}_trimmed_64.fq") into trim_out69793218560
		set val(id),  path("${id}_trimmed_65.fq") into trim_out70866960384
		set val(id),  path("${id}_trimmed_66.fq") into trim_out71940702208
		set val(id),  path("${id}_trimmed_67.fq") into trim_out73014444032
		set val(id),  path("${id}_trimmed_68.fq") into trim_out74088185856
		set val(id),  path("${id}_trimmed_69.fq") into trim_out75161927680
		set val(id),  path("${id}_trimmed_70.fq") into trim_out76235669504
		set val(id),  path("${id}_trimmed_71.fq") into trim_out77309411328
		set val(id),  path("${id}_trimmed_72.fq") into trim_out78383153152
		set val(id),  path("${id}_trimmed_73.fq") into trim_out79456894976
		set val(id),  path("${id}_trimmed_74.fq") into trim_out80530636800
		set val(id),  path("${id}_trimmed_75.fq") into trim_out81604378624
		set val(id),  path("${id}_trimmed_76.fq") into trim_out82678120448
		set val(id),  path("${id}_trimmed_77.fq") into trim_out83751862272

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
		NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq
		NanoFilt ${id}_3.fastq -l 500 -q 10 > ${id}_trimmed_3.fq
		NanoFilt ${id}_4.fastq -l 500 -q 10 > ${id}_trimmed_4.fq
		NanoFilt ${id}_5.fastq -l 500 -q 10 > ${id}_trimmed_5.fq
		NanoFilt ${id}_6.fastq -l 500 -q 10 > ${id}_trimmed_6.fq
		NanoFilt ${id}_7.fastq -l 500 -q 10 > ${id}_trimmed_7.fq
		NanoFilt ${id}_8.fastq -l 500 -q 10 > ${id}_trimmed_8.fq
		NanoFilt ${id}_9.fastq -l 500 -q 10 > ${id}_trimmed_9.fq
		NanoFilt ${id}_10.fastq -l 500 -q 10 > ${id}_trimmed_10.fq
		NanoFilt ${id}_11.fastq -l 500 -q 10 > ${id}_trimmed_11.fq
		NanoFilt ${id}_12.fastq -l 500 -q 10 > ${id}_trimmed_12.fq
		NanoFilt ${id}_13.fastq -l 500 -q 10 > ${id}_trimmed_13.fq
		NanoFilt ${id}_14.fastq -l 500 -q 10 > ${id}_trimmed_14.fq
		NanoFilt ${id}_15.fastq -l 500 -q 10 > ${id}_trimmed_15.fq
		NanoFilt ${id}_16.fastq -l 500 -q 10 > ${id}_trimmed_16.fq
		NanoFilt ${id}_17.fastq -l 500 -q 10 > ${id}_trimmed_17.fq
		NanoFilt ${id}_18.fastq -l 500 -q 10 > ${id}_trimmed_18.fq
		NanoFilt ${id}_19.fastq -l 500 -q 10 > ${id}_trimmed_19.fq
		NanoFilt ${id}_20.fastq -l 500 -q 10 > ${id}_trimmed_20.fq
		NanoFilt ${id}_21.fastq -l 500 -q 10 > ${id}_trimmed_21.fq
		NanoFilt ${id}_22.fastq -l 500 -q 10 > ${id}_trimmed_22.fq
		NanoFilt ${id}_23.fastq -l 500 -q 10 > ${id}_trimmed_23.fq
		NanoFilt ${id}_24.fastq -l 500 -q 10 > ${id}_trimmed_24.fq
		NanoFilt ${id}_25.fastq -l 500 -q 10 > ${id}_trimmed_25.fq
		NanoFilt ${id}_26.fastq -l 500 -q 10 > ${id}_trimmed_26.fq
		NanoFilt ${id}_27.fastq -l 500 -q 10 > ${id}_trimmed_27.fq
		NanoFilt ${id}_28.fastq -l 500 -q 10 > ${id}_trimmed_28.fq
		NanoFilt ${id}_29.fastq -l 500 -q 10 > ${id}_trimmed_29.fq
		NanoFilt ${id}_30.fastq -l 500 -q 10 > ${id}_trimmed_30.fq
		NanoFilt ${id}_31.fastq -l 500 -q 10 > ${id}_trimmed_31.fq
		NanoFilt ${id}_32.fastq -l 500 -q 10 > ${id}_trimmed_32.fq
		NanoFilt ${id}_33.fastq -l 500 -q 10 > ${id}_trimmed_33.fq
		NanoFilt ${id}_34.fastq -l 500 -q 10 > ${id}_trimmed_34.fq
		NanoFilt ${id}_35.fastq -l 500 -q 10 > ${id}_trimmed_35.fq
		NanoFilt ${id}_36.fastq -l 500 -q 10 > ${id}_trimmed_36.fq
		NanoFilt ${id}_37.fastq -l 500 -q 10 > ${id}_trimmed_37.fq
		NanoFilt ${id}_38.fastq -l 500 -q 10 > ${id}_trimmed_38.fq
		NanoFilt ${id}_39.fastq -l 500 -q 10 > ${id}_trimmed_39.fq
		NanoFilt ${id}_40.fastq -l 500 -q 10 > ${id}_trimmed_40.fq
		NanoFilt ${id}_41.fastq -l 500 -q 10 > ${id}_trimmed_41.fq
		NanoFilt ${id}_42.fastq -l 500 -q 10 > ${id}_trimmed_42.fq
		NanoFilt ${id}_43.fastq -l 500 -q 10 > ${id}_trimmed_43.fq
		NanoFilt ${id}_44.fastq -l 500 -q 10 > ${id}_trimmed_44.fq
		NanoFilt ${id}_45.fastq -l 500 -q 10 > ${id}_trimmed_45.fq
		NanoFilt ${id}_46.fastq -l 500 -q 10 > ${id}_trimmed_46.fq
		NanoFilt ${id}_47.fastq -l 500 -q 10 > ${id}_trimmed_47.fq
		NanoFilt ${id}_48.fastq -l 500 -q 10 > ${id}_trimmed_48.fq
		NanoFilt ${id}_49.fastq -l 500 -q 10 > ${id}_trimmed_49.fq
		NanoFilt ${id}_50.fastq -l 500 -q 10 > ${id}_trimmed_50.fq
		NanoFilt ${id}_51.fastq -l 500 -q 10 > ${id}_trimmed_51.fq
		NanoFilt ${id}_52.fastq -l 500 -q 10 > ${id}_trimmed_52.fq
		NanoFilt ${id}_53.fastq -l 500 -q 10 > ${id}_trimmed_53.fq
		NanoFilt ${id}_54.fastq -l 500 -q 10 > ${id}_trimmed_54.fq
		NanoFilt ${id}_55.fastq -l 500 -q 10 > ${id}_trimmed_55.fq
		NanoFilt ${id}_56.fastq -l 500 -q 10 > ${id}_trimmed_56.fq
		NanoFilt ${id}_57.fastq -l 500 -q 10 > ${id}_trimmed_57.fq
		NanoFilt ${id}_58.fastq -l 500 -q 10 > ${id}_trimmed_58.fq
		NanoFilt ${id}_59.fastq -l 500 -q 10 > ${id}_trimmed_59.fq
		NanoFilt ${id}_60.fastq -l 500 -q 10 > ${id}_trimmed_60.fq
		NanoFilt ${id}_61.fastq -l 500 -q 10 > ${id}_trimmed_61.fq
		NanoFilt ${id}_62.fastq -l 500 -q 10 > ${id}_trimmed_62.fq
		NanoFilt ${id}_63.fastq -l 500 -q 10 > ${id}_trimmed_63.fq
		NanoFilt ${id}_64.fastq -l 500 -q 10 > ${id}_trimmed_64.fq
		NanoFilt ${id}_65.fastq -l 500 -q 10 > ${id}_trimmed_65.fq
		NanoFilt ${id}_66.fastq -l 500 -q 10 > ${id}_trimmed_66.fq
		NanoFilt ${id}_67.fastq -l 500 -q 10 > ${id}_trimmed_67.fq
		NanoFilt ${id}_68.fastq -l 500 -q 10 > ${id}_trimmed_68.fq
		NanoFilt ${id}_69.fastq -l 500 -q 10 > ${id}_trimmed_69.fq
		NanoFilt ${id}_70.fastq -l 500 -q 10 > ${id}_trimmed_70.fq
		NanoFilt ${id}_71.fastq -l 500 -q 10 > ${id}_trimmed_71.fq
		NanoFilt ${id}_72.fastq -l 500 -q 10 > ${id}_trimmed_72.fq
		NanoFilt ${id}_73.fastq -l 500 -q 10 > ${id}_trimmed_73.fq
		NanoFilt ${id}_74.fastq -l 500 -q 10 > ${id}_trimmed_74.fq
		NanoFilt ${id}_75.fastq -l 500 -q 10 > ${id}_trimmed_75.fq
		NanoFilt ${id}_76.fastq -l 500 -q 10 > ${id}_trimmed_76.fq
		NanoFilt ${id}_77.fastq -l 500 -q 10 > ${id}_trimmed_77.fq
		"""
        




}

process maptoreference8 {
    publishDir "$params.output.folder8/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder8/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out1073741824
		set val(id), path(trim_read1) from trim_out2147483648
		set val(id), path(trim_read2) from trim_out3221225472
		set val(id), path(trim_read3) from trim_out4294967296
		set val(id), path(trim_read4) from trim_out5368709120
		set val(id), path(trim_read5) from trim_out6442450944
		set val(id), path(trim_read6) from trim_out7516192768
		set val(id), path(trim_read7) from trim_out8589934592
		set val(id), path(trim_read8) from trim_out9663676416
		set val(id), path(trim_read9) from trim_out10737418240
		set val(id), path(trim_read10) from trim_out11811160064
		set val(id), path(trim_read11) from trim_out12884901888
		set val(id), path(trim_read12) from trim_out13958643712
		set val(id), path(trim_read13) from trim_out15032385536
		set val(id), path(trim_read14) from trim_out16106127360
		set val(id), path(trim_read15) from trim_out17179869184
		set val(id), path(trim_read16) from trim_out18253611008
		set val(id), path(trim_read17) from trim_out19327352832
		set val(id), path(trim_read18) from trim_out20401094656
		set val(id), path(trim_read19) from trim_out21474836480
		set val(id), path(trim_read20) from trim_out22548578304
		set val(id), path(trim_read21) from trim_out23622320128
		set val(id), path(trim_read22) from trim_out24696061952
		set val(id), path(trim_read23) from trim_out25769803776
		set val(id), path(trim_read24) from trim_out26843545600
		set val(id), path(trim_read25) from trim_out27917287424
		set val(id), path(trim_read26) from trim_out28991029248
		set val(id), path(trim_read27) from trim_out30064771072
		set val(id), path(trim_read28) from trim_out31138512896
		set val(id), path(trim_read29) from trim_out32212254720
		set val(id), path(trim_read30) from trim_out33285996544
		set val(id), path(trim_read31) from trim_out34359738368
		set val(id), path(trim_read32) from trim_out35433480192
		set val(id), path(trim_read33) from trim_out36507222016
		set val(id), path(trim_read34) from trim_out37580963840
		set val(id), path(trim_read35) from trim_out38654705664
		set val(id), path(trim_read36) from trim_out39728447488
		set val(id), path(trim_read37) from trim_out40802189312
		set val(id), path(trim_read38) from trim_out41875931136
		set val(id), path(trim_read39) from trim_out42949672960
		set val(id), path(trim_read40) from trim_out44023414784
		set val(id), path(trim_read41) from trim_out45097156608
		set val(id), path(trim_read42) from trim_out46170898432
		set val(id), path(trim_read43) from trim_out47244640256
		set val(id), path(trim_read44) from trim_out48318382080
		set val(id), path(trim_read45) from trim_out49392123904
		set val(id), path(trim_read46) from trim_out50465865728
		set val(id), path(trim_read47) from trim_out51539607552
		set val(id), path(trim_read48) from trim_out52613349376
		set val(id), path(trim_read49) from trim_out53687091200
		set val(id), path(trim_read50) from trim_out54760833024
		set val(id), path(trim_read51) from trim_out55834574848
		set val(id), path(trim_read52) from trim_out56908316672
		set val(id), path(trim_read53) from trim_out57982058496
		set val(id), path(trim_read54) from trim_out59055800320
		set val(id), path(trim_read55) from trim_out60129542144
		set val(id), path(trim_read56) from trim_out61203283968
		set val(id), path(trim_read57) from trim_out62277025792
		set val(id), path(trim_read58) from trim_out63350767616
		set val(id), path(trim_read59) from trim_out64424509440
		set val(id), path(trim_read60) from trim_out65498251264
		set val(id), path(trim_read61) from trim_out66571993088
		set val(id), path(trim_read62) from trim_out67645734912
		set val(id), path(trim_read63) from trim_out68719476736
		set val(id), path(trim_read64) from trim_out69793218560
		set val(id), path(trim_read65) from trim_out70866960384
		set val(id), path(trim_read66) from trim_out71940702208
		set val(id), path(trim_read67) from trim_out73014444032
		set val(id), path(trim_read68) from trim_out74088185856
		set val(id), path(trim_read69) from trim_out75161927680
		set val(id), path(trim_read70) from trim_out76235669504
		set val(id), path(trim_read71) from trim_out77309411328
		set val(id), path(trim_read72) from trim_out78383153152
		set val(id), path(trim_read73) from trim_out79456894976
		set val(id), path(trim_read74) from trim_out80530636800
		set val(id), path(trim_read75) from trim_out81604378624
		set val(id), path(trim_read76) from trim_out82678120448
		set val(id), path(trim_read77) from trim_out83751862272

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1073741824
		set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2147483648
		set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3221225472
		set val(id), path("${id}_mapped_3.fq"), path("${id}_unmapped_3.fq") into mapped_out4294967296
		set val(id), path("${id}_mapped_4.fq"), path("${id}_unmapped_4.fq") into mapped_out5368709120
		set val(id), path("${id}_mapped_5.fq"), path("${id}_unmapped_5.fq") into mapped_out6442450944
		set val(id), path("${id}_mapped_6.fq"), path("${id}_unmapped_6.fq") into mapped_out7516192768
		set val(id), path("${id}_mapped_7.fq"), path("${id}_unmapped_7.fq") into mapped_out8589934592
		set val(id), path("${id}_mapped_8.fq"), path("${id}_unmapped_8.fq") into mapped_out9663676416
		set val(id), path("${id}_mapped_9.fq"), path("${id}_unmapped_9.fq") into mapped_out10737418240
		set val(id), path("${id}_mapped_10.fq"), path("${id}_unmapped_10.fq") into mapped_out11811160064
		set val(id), path("${id}_mapped_11.fq"), path("${id}_unmapped_11.fq") into mapped_out12884901888
		set val(id), path("${id}_mapped_12.fq"), path("${id}_unmapped_12.fq") into mapped_out13958643712
		set val(id), path("${id}_mapped_13.fq"), path("${id}_unmapped_13.fq") into mapped_out15032385536
		set val(id), path("${id}_mapped_14.fq"), path("${id}_unmapped_14.fq") into mapped_out16106127360
		set val(id), path("${id}_mapped_15.fq"), path("${id}_unmapped_15.fq") into mapped_out17179869184
		set val(id), path("${id}_mapped_16.fq"), path("${id}_unmapped_16.fq") into mapped_out18253611008
		set val(id), path("${id}_mapped_17.fq"), path("${id}_unmapped_17.fq") into mapped_out19327352832
		set val(id), path("${id}_mapped_18.fq"), path("${id}_unmapped_18.fq") into mapped_out20401094656
		set val(id), path("${id}_mapped_19.fq"), path("${id}_unmapped_19.fq") into mapped_out21474836480
		set val(id), path("${id}_mapped_20.fq"), path("${id}_unmapped_20.fq") into mapped_out22548578304
		set val(id), path("${id}_mapped_21.fq"), path("${id}_unmapped_21.fq") into mapped_out23622320128
		set val(id), path("${id}_mapped_22.fq"), path("${id}_unmapped_22.fq") into mapped_out24696061952
		set val(id), path("${id}_mapped_23.fq"), path("${id}_unmapped_23.fq") into mapped_out25769803776
		set val(id), path("${id}_mapped_24.fq"), path("${id}_unmapped_24.fq") into mapped_out26843545600
		set val(id), path("${id}_mapped_25.fq"), path("${id}_unmapped_25.fq") into mapped_out27917287424
		set val(id), path("${id}_mapped_26.fq"), path("${id}_unmapped_26.fq") into mapped_out28991029248
		set val(id), path("${id}_mapped_27.fq"), path("${id}_unmapped_27.fq") into mapped_out30064771072
		set val(id), path("${id}_mapped_28.fq"), path("${id}_unmapped_28.fq") into mapped_out31138512896
		set val(id), path("${id}_mapped_29.fq"), path("${id}_unmapped_29.fq") into mapped_out32212254720
		set val(id), path("${id}_mapped_30.fq"), path("${id}_unmapped_30.fq") into mapped_out33285996544
		set val(id), path("${id}_mapped_31.fq"), path("${id}_unmapped_31.fq") into mapped_out34359738368
		set val(id), path("${id}_mapped_32.fq"), path("${id}_unmapped_32.fq") into mapped_out35433480192
		set val(id), path("${id}_mapped_33.fq"), path("${id}_unmapped_33.fq") into mapped_out36507222016
		set val(id), path("${id}_mapped_34.fq"), path("${id}_unmapped_34.fq") into mapped_out37580963840
		set val(id), path("${id}_mapped_35.fq"), path("${id}_unmapped_35.fq") into mapped_out38654705664
		set val(id), path("${id}_mapped_36.fq"), path("${id}_unmapped_36.fq") into mapped_out39728447488
		set val(id), path("${id}_mapped_37.fq"), path("${id}_unmapped_37.fq") into mapped_out40802189312
		set val(id), path("${id}_mapped_38.fq"), path("${id}_unmapped_38.fq") into mapped_out41875931136
		set val(id), path("${id}_mapped_39.fq"), path("${id}_unmapped_39.fq") into mapped_out42949672960
		set val(id), path("${id}_mapped_40.fq"), path("${id}_unmapped_40.fq") into mapped_out44023414784
		set val(id), path("${id}_mapped_41.fq"), path("${id}_unmapped_41.fq") into mapped_out45097156608
		set val(id), path("${id}_mapped_42.fq"), path("${id}_unmapped_42.fq") into mapped_out46170898432
		set val(id), path("${id}_mapped_43.fq"), path("${id}_unmapped_43.fq") into mapped_out47244640256
		set val(id), path("${id}_mapped_44.fq"), path("${id}_unmapped_44.fq") into mapped_out48318382080
		set val(id), path("${id}_mapped_45.fq"), path("${id}_unmapped_45.fq") into mapped_out49392123904
		set val(id), path("${id}_mapped_46.fq"), path("${id}_unmapped_46.fq") into mapped_out50465865728
		set val(id), path("${id}_mapped_47.fq"), path("${id}_unmapped_47.fq") into mapped_out51539607552
		set val(id), path("${id}_mapped_48.fq"), path("${id}_unmapped_48.fq") into mapped_out52613349376
		set val(id), path("${id}_mapped_49.fq"), path("${id}_unmapped_49.fq") into mapped_out53687091200
		set val(id), path("${id}_mapped_50.fq"), path("${id}_unmapped_50.fq") into mapped_out54760833024
		set val(id), path("${id}_mapped_51.fq"), path("${id}_unmapped_51.fq") into mapped_out55834574848
		set val(id), path("${id}_mapped_52.fq"), path("${id}_unmapped_52.fq") into mapped_out56908316672
		set val(id), path("${id}_mapped_53.fq"), path("${id}_unmapped_53.fq") into mapped_out57982058496
		set val(id), path("${id}_mapped_54.fq"), path("${id}_unmapped_54.fq") into mapped_out59055800320
		set val(id), path("${id}_mapped_55.fq"), path("${id}_unmapped_55.fq") into mapped_out60129542144
		set val(id), path("${id}_mapped_56.fq"), path("${id}_unmapped_56.fq") into mapped_out61203283968
		set val(id), path("${id}_mapped_57.fq"), path("${id}_unmapped_57.fq") into mapped_out62277025792
		set val(id), path("${id}_mapped_58.fq"), path("${id}_unmapped_58.fq") into mapped_out63350767616
		set val(id), path("${id}_mapped_59.fq"), path("${id}_unmapped_59.fq") into mapped_out64424509440
		set val(id), path("${id}_mapped_60.fq"), path("${id}_unmapped_60.fq") into mapped_out65498251264
		set val(id), path("${id}_mapped_61.fq"), path("${id}_unmapped_61.fq") into mapped_out66571993088
		set val(id), path("${id}_mapped_62.fq"), path("${id}_unmapped_62.fq") into mapped_out67645734912
		set val(id), path("${id}_mapped_63.fq"), path("${id}_unmapped_63.fq") into mapped_out68719476736
		set val(id), path("${id}_mapped_64.fq"), path("${id}_unmapped_64.fq") into mapped_out69793218560
		set val(id), path("${id}_mapped_65.fq"), path("${id}_unmapped_65.fq") into mapped_out70866960384
		set val(id), path("${id}_mapped_66.fq"), path("${id}_unmapped_66.fq") into mapped_out71940702208
		set val(id), path("${id}_mapped_67.fq"), path("${id}_unmapped_67.fq") into mapped_out73014444032
		set val(id), path("${id}_mapped_68.fq"), path("${id}_unmapped_68.fq") into mapped_out74088185856
		set val(id), path("${id}_mapped_69.fq"), path("${id}_unmapped_69.fq") into mapped_out75161927680
		set val(id), path("${id}_mapped_70.fq"), path("${id}_unmapped_70.fq") into mapped_out76235669504
		set val(id), path("${id}_mapped_71.fq"), path("${id}_unmapped_71.fq") into mapped_out77309411328
		set val(id), path("${id}_mapped_72.fq"), path("${id}_unmapped_72.fq") into mapped_out78383153152
		set val(id), path("${id}_mapped_73.fq"), path("${id}_unmapped_73.fq") into mapped_out79456894976
		set val(id), path("${id}_mapped_74.fq"), path("${id}_unmapped_74.fq") into mapped_out80530636800
		set val(id), path("${id}_mapped_75.fq"), path("${id}_unmapped_75.fq") into mapped_out81604378624
		set val(id), path("${id}_mapped_76.fq"), path("${id}_unmapped_76.fq") into mapped_out82678120448
		set val(id), path("${id}_mapped_77.fq"), path("${id}_unmapped_77.fq") into mapped_out83751862272

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_3.fq outm=${id}_mapped_3.fq outu=${id}_unmapped_3.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_4.fq outm=${id}_mapped_4.fq outu=${id}_unmapped_4.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_5.fq outm=${id}_mapped_5.fq outu=${id}_unmapped_5.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_6.fq outm=${id}_mapped_6.fq outu=${id}_unmapped_6.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_7.fq outm=${id}_mapped_7.fq outu=${id}_unmapped_7.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_8.fq outm=${id}_mapped_8.fq outu=${id}_unmapped_8.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_9.fq outm=${id}_mapped_9.fq outu=${id}_unmapped_9.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_10.fq outm=${id}_mapped_10.fq outu=${id}_unmapped_10.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_11.fq outm=${id}_mapped_11.fq outu=${id}_unmapped_11.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_12.fq outm=${id}_mapped_12.fq outu=${id}_unmapped_12.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_13.fq outm=${id}_mapped_13.fq outu=${id}_unmapped_13.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_14.fq outm=${id}_mapped_14.fq outu=${id}_unmapped_14.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_15.fq outm=${id}_mapped_15.fq outu=${id}_unmapped_15.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_16.fq outm=${id}_mapped_16.fq outu=${id}_unmapped_16.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_17.fq outm=${id}_mapped_17.fq outu=${id}_unmapped_17.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_18.fq outm=${id}_mapped_18.fq outu=${id}_unmapped_18.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_19.fq outm=${id}_mapped_19.fq outu=${id}_unmapped_19.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_20.fq outm=${id}_mapped_20.fq outu=${id}_unmapped_20.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_21.fq outm=${id}_mapped_21.fq outu=${id}_unmapped_21.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_22.fq outm=${id}_mapped_22.fq outu=${id}_unmapped_22.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_23.fq outm=${id}_mapped_23.fq outu=${id}_unmapped_23.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_24.fq outm=${id}_mapped_24.fq outu=${id}_unmapped_24.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_25.fq outm=${id}_mapped_25.fq outu=${id}_unmapped_25.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_26.fq outm=${id}_mapped_26.fq outu=${id}_unmapped_26.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_27.fq outm=${id}_mapped_27.fq outu=${id}_unmapped_27.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_28.fq outm=${id}_mapped_28.fq outu=${id}_unmapped_28.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_29.fq outm=${id}_mapped_29.fq outu=${id}_unmapped_29.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_30.fq outm=${id}_mapped_30.fq outu=${id}_unmapped_30.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_31.fq outm=${id}_mapped_31.fq outu=${id}_unmapped_31.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_32.fq outm=${id}_mapped_32.fq outu=${id}_unmapped_32.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_33.fq outm=${id}_mapped_33.fq outu=${id}_unmapped_33.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_34.fq outm=${id}_mapped_34.fq outu=${id}_unmapped_34.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_35.fq outm=${id}_mapped_35.fq outu=${id}_unmapped_35.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_36.fq outm=${id}_mapped_36.fq outu=${id}_unmapped_36.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_37.fq outm=${id}_mapped_37.fq outu=${id}_unmapped_37.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_38.fq outm=${id}_mapped_38.fq outu=${id}_unmapped_38.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_39.fq outm=${id}_mapped_39.fq outu=${id}_unmapped_39.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_40.fq outm=${id}_mapped_40.fq outu=${id}_unmapped_40.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_41.fq outm=${id}_mapped_41.fq outu=${id}_unmapped_41.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_42.fq outm=${id}_mapped_42.fq outu=${id}_unmapped_42.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_43.fq outm=${id}_mapped_43.fq outu=${id}_unmapped_43.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_44.fq outm=${id}_mapped_44.fq outu=${id}_unmapped_44.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_45.fq outm=${id}_mapped_45.fq outu=${id}_unmapped_45.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_46.fq outm=${id}_mapped_46.fq outu=${id}_unmapped_46.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_47.fq outm=${id}_mapped_47.fq outu=${id}_unmapped_47.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_48.fq outm=${id}_mapped_48.fq outu=${id}_unmapped_48.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_49.fq outm=${id}_mapped_49.fq outu=${id}_unmapped_49.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_50.fq outm=${id}_mapped_50.fq outu=${id}_unmapped_50.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_51.fq outm=${id}_mapped_51.fq outu=${id}_unmapped_51.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_52.fq outm=${id}_mapped_52.fq outu=${id}_unmapped_52.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_53.fq outm=${id}_mapped_53.fq outu=${id}_unmapped_53.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_54.fq outm=${id}_mapped_54.fq outu=${id}_unmapped_54.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_55.fq outm=${id}_mapped_55.fq outu=${id}_unmapped_55.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_56.fq outm=${id}_mapped_56.fq outu=${id}_unmapped_56.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_57.fq outm=${id}_mapped_57.fq outu=${id}_unmapped_57.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_58.fq outm=${id}_mapped_58.fq outu=${id}_unmapped_58.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_59.fq outm=${id}_mapped_59.fq outu=${id}_unmapped_59.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_60.fq outm=${id}_mapped_60.fq outu=${id}_unmapped_60.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_61.fq outm=${id}_mapped_61.fq outu=${id}_unmapped_61.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_62.fq outm=${id}_mapped_62.fq outu=${id}_unmapped_62.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_63.fq outm=${id}_mapped_63.fq outu=${id}_unmapped_63.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_64.fq outm=${id}_mapped_64.fq outu=${id}_unmapped_64.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_65.fq outm=${id}_mapped_65.fq outu=${id}_unmapped_65.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_66.fq outm=${id}_mapped_66.fq outu=${id}_unmapped_66.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_67.fq outm=${id}_mapped_67.fq outu=${id}_unmapped_67.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_68.fq outm=${id}_mapped_68.fq outu=${id}_unmapped_68.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_69.fq outm=${id}_mapped_69.fq outu=${id}_unmapped_69.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_70.fq outm=${id}_mapped_70.fq outu=${id}_unmapped_70.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_71.fq outm=${id}_mapped_71.fq outu=${id}_unmapped_71.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_72.fq outm=${id}_mapped_72.fq outu=${id}_unmapped_72.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_73.fq outm=${id}_mapped_73.fq outu=${id}_unmapped_73.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_74.fq outm=${id}_mapped_74.fq outu=${id}_unmapped_74.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_75.fq outm=${id}_mapped_75.fq outu=${id}_unmapped_75.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_76.fq outm=${id}_mapped_76.fq outu=${id}_unmapped_76.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_77.fq outm=${id}_mapped_77.fq outu=${id}_unmapped_77.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble8 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder8/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1073741824
		set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2147483648
		set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3221225472
		set val(sample), path(mapped_read3), path(unmapped_read3) from mapped_out4294967296
		set val(sample), path(mapped_read4), path(unmapped_read4) from mapped_out5368709120
		set val(sample), path(mapped_read5), path(unmapped_read5) from mapped_out6442450944
		set val(sample), path(mapped_read6), path(unmapped_read6) from mapped_out7516192768
		set val(sample), path(mapped_read7), path(unmapped_read7) from mapped_out8589934592
		set val(sample), path(mapped_read8), path(unmapped_read8) from mapped_out9663676416
		set val(sample), path(mapped_read9), path(unmapped_read9) from mapped_out10737418240
		set val(sample), path(mapped_read10), path(unmapped_read10) from mapped_out11811160064
		set val(sample), path(mapped_read11), path(unmapped_read11) from mapped_out12884901888
		set val(sample), path(mapped_read12), path(unmapped_read12) from mapped_out13958643712
		set val(sample), path(mapped_read13), path(unmapped_read13) from mapped_out15032385536
		set val(sample), path(mapped_read14), path(unmapped_read14) from mapped_out16106127360
		set val(sample), path(mapped_read15), path(unmapped_read15) from mapped_out17179869184
		set val(sample), path(mapped_read16), path(unmapped_read16) from mapped_out18253611008
		set val(sample), path(mapped_read17), path(unmapped_read17) from mapped_out19327352832
		set val(sample), path(mapped_read18), path(unmapped_read18) from mapped_out20401094656
		set val(sample), path(mapped_read19), path(unmapped_read19) from mapped_out21474836480
		set val(sample), path(mapped_read20), path(unmapped_read20) from mapped_out22548578304
		set val(sample), path(mapped_read21), path(unmapped_read21) from mapped_out23622320128
		set val(sample), path(mapped_read22), path(unmapped_read22) from mapped_out24696061952
		set val(sample), path(mapped_read23), path(unmapped_read23) from mapped_out25769803776
		set val(sample), path(mapped_read24), path(unmapped_read24) from mapped_out26843545600
		set val(sample), path(mapped_read25), path(unmapped_read25) from mapped_out27917287424
		set val(sample), path(mapped_read26), path(unmapped_read26) from mapped_out28991029248
		set val(sample), path(mapped_read27), path(unmapped_read27) from mapped_out30064771072
		set val(sample), path(mapped_read28), path(unmapped_read28) from mapped_out31138512896
		set val(sample), path(mapped_read29), path(unmapped_read29) from mapped_out32212254720
		set val(sample), path(mapped_read30), path(unmapped_read30) from mapped_out33285996544
		set val(sample), path(mapped_read31), path(unmapped_read31) from mapped_out34359738368
		set val(sample), path(mapped_read32), path(unmapped_read32) from mapped_out35433480192
		set val(sample), path(mapped_read33), path(unmapped_read33) from mapped_out36507222016
		set val(sample), path(mapped_read34), path(unmapped_read34) from mapped_out37580963840
		set val(sample), path(mapped_read35), path(unmapped_read35) from mapped_out38654705664
		set val(sample), path(mapped_read36), path(unmapped_read36) from mapped_out39728447488
		set val(sample), path(mapped_read37), path(unmapped_read37) from mapped_out40802189312
		set val(sample), path(mapped_read38), path(unmapped_read38) from mapped_out41875931136
		set val(sample), path(mapped_read39), path(unmapped_read39) from mapped_out42949672960
		set val(sample), path(mapped_read40), path(unmapped_read40) from mapped_out44023414784
		set val(sample), path(mapped_read41), path(unmapped_read41) from mapped_out45097156608
		set val(sample), path(mapped_read42), path(unmapped_read42) from mapped_out46170898432
		set val(sample), path(mapped_read43), path(unmapped_read43) from mapped_out47244640256
		set val(sample), path(mapped_read44), path(unmapped_read44) from mapped_out48318382080
		set val(sample), path(mapped_read45), path(unmapped_read45) from mapped_out49392123904
		set val(sample), path(mapped_read46), path(unmapped_read46) from mapped_out50465865728
		set val(sample), path(mapped_read47), path(unmapped_read47) from mapped_out51539607552
		set val(sample), path(mapped_read48), path(unmapped_read48) from mapped_out52613349376
		set val(sample), path(mapped_read49), path(unmapped_read49) from mapped_out53687091200
		set val(sample), path(mapped_read50), path(unmapped_read50) from mapped_out54760833024
		set val(sample), path(mapped_read51), path(unmapped_read51) from mapped_out55834574848
		set val(sample), path(mapped_read52), path(unmapped_read52) from mapped_out56908316672
		set val(sample), path(mapped_read53), path(unmapped_read53) from mapped_out57982058496
		set val(sample), path(mapped_read54), path(unmapped_read54) from mapped_out59055800320
		set val(sample), path(mapped_read55), path(unmapped_read55) from mapped_out60129542144
		set val(sample), path(mapped_read56), path(unmapped_read56) from mapped_out61203283968
		set val(sample), path(mapped_read57), path(unmapped_read57) from mapped_out62277025792
		set val(sample), path(mapped_read58), path(unmapped_read58) from mapped_out63350767616
		set val(sample), path(mapped_read59), path(unmapped_read59) from mapped_out64424509440
		set val(sample), path(mapped_read60), path(unmapped_read60) from mapped_out65498251264
		set val(sample), path(mapped_read61), path(unmapped_read61) from mapped_out66571993088
		set val(sample), path(mapped_read62), path(unmapped_read62) from mapped_out67645734912
		set val(sample), path(mapped_read63), path(unmapped_read63) from mapped_out68719476736
		set val(sample), path(mapped_read64), path(unmapped_read64) from mapped_out69793218560
		set val(sample), path(mapped_read65), path(unmapped_read65) from mapped_out70866960384
		set val(sample), path(mapped_read66), path(unmapped_read66) from mapped_out71940702208
		set val(sample), path(mapped_read67), path(unmapped_read67) from mapped_out73014444032
		set val(sample), path(mapped_read68), path(unmapped_read68) from mapped_out74088185856
		set val(sample), path(mapped_read69), path(unmapped_read69) from mapped_out75161927680
		set val(sample), path(mapped_read70), path(unmapped_read70) from mapped_out76235669504
		set val(sample), path(mapped_read71), path(unmapped_read71) from mapped_out77309411328
		set val(sample), path(mapped_read72), path(unmapped_read72) from mapped_out78383153152
		set val(sample), path(mapped_read73), path(unmapped_read73) from mapped_out79456894976
		set val(sample), path(mapped_read74), path(unmapped_read74) from mapped_out80530636800
		set val(sample), path(mapped_read75), path(unmapped_read75) from mapped_out81604378624
		set val(sample), path(mapped_read76), path(unmapped_read76) from mapped_out82678120448
		set val(sample), path(mapped_read77), path(unmapped_read77) from mapped_out83751862272


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1073741824
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq  -s ${sample}_mapped_1.fq  -s ${sample}_mapped_2.fq  -s ${sample}_mapped_3.fq  -s ${sample}_mapped_4.fq  -s ${sample}_mapped_5.fq  -s ${sample}_mapped_6.fq  -s ${sample}_mapped_7.fq  -s ${sample}_mapped_8.fq  -s ${sample}_mapped_9.fq  -s ${sample}_mapped_10.fq  -s ${sample}_mapped_11.fq  -s ${sample}_mapped_12.fq  -s ${sample}_mapped_13.fq  -s ${sample}_mapped_14.fq  -s ${sample}_mapped_15.fq  -s ${sample}_mapped_16.fq  -s ${sample}_mapped_17.fq  -s ${sample}_mapped_18.fq  -s ${sample}_mapped_19.fq  -s ${sample}_mapped_20.fq  -s ${sample}_mapped_21.fq  -s ${sample}_mapped_22.fq  -s ${sample}_mapped_23.fq  -s ${sample}_mapped_24.fq  -s ${sample}_mapped_25.fq  -s ${sample}_mapped_26.fq  -s ${sample}_mapped_27.fq  -s ${sample}_mapped_28.fq  -s ${sample}_mapped_29.fq  -s ${sample}_mapped_30.fq  -s ${sample}_mapped_31.fq  -s ${sample}_mapped_32.fq  -s ${sample}_mapped_33.fq  -s ${sample}_mapped_34.fq  -s ${sample}_mapped_35.fq  -s ${sample}_mapped_36.fq  -s ${sample}_mapped_37.fq  -s ${sample}_mapped_38.fq  -s ${sample}_mapped_39.fq  -s ${sample}_mapped_40.fq  -s ${sample}_mapped_41.fq  -s ${sample}_mapped_42.fq  -s ${sample}_mapped_43.fq  -s ${sample}_mapped_44.fq  -s ${sample}_mapped_45.fq  -s ${sample}_mapped_46.fq  -s ${sample}_mapped_47.fq  -s ${sample}_mapped_48.fq  -s ${sample}_mapped_49.fq  -s ${sample}_mapped_50.fq  -s ${sample}_mapped_51.fq  -s ${sample}_mapped_52.fq  -s ${sample}_mapped_53.fq  -s ${sample}_mapped_54.fq  -s ${sample}_mapped_55.fq  -s ${sample}_mapped_56.fq  -s ${sample}_mapped_57.fq  -s ${sample}_mapped_58.fq  -s ${sample}_mapped_59.fq  -s ${sample}_mapped_60.fq  -s ${sample}_mapped_61.fq  -s ${sample}_mapped_62.fq  -s ${sample}_mapped_63.fq  -s ${sample}_mapped_64.fq  -s ${sample}_mapped_65.fq  -s ${sample}_mapped_66.fq  -s ${sample}_mapped_67.fq  -s ${sample}_mapped_68.fq  -s ${sample}_mapped_69.fq  -s ${sample}_mapped_70.fq  -s ${sample}_mapped_71.fq  -s ${sample}_mapped_72.fq  -s ${sample}_mapped_73.fq  -s ${sample}_mapped_74.fq  -s ${sample}_mapped_75.fq  -s ${sample}_mapped_76.fq  -s ${sample}_mapped_77.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate8 {

    publishDir "$params.output.folder8/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1073741824
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1073741824

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern8 {

    publishDir "$params.output.folder8/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1073741824
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads9 = params.input.fastq_path9+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads9,size:1)
params.genome9 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt9 {
    publishDir "$params.output.folder9/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out3486784401

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		"""
        




}

process maptoreference9 {
    publishDir "$params.output.folder9/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder9/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out3486784401

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out3486784401

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble9 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder9/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out3486784401


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout3486784401
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate9 {

    publishDir "$params.output.folder9/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout3486784401
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq3486784401

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern9 {

    publishDir "$params.output.folder9/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq3486784401
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads10 = params.input.fastq_path10+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads10,size:1)
params.genome10 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt10 {
    publishDir "$params.output.folder10/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out10000000000

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		"""
        




}

process maptoreference10 {
    publishDir "$params.output.folder10/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder10/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out10000000000

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out10000000000

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble10 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder10/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out10000000000


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout10000000000
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate10 {

    publishDir "$params.output.folder10/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout10000000000
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq10000000000

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern10 {

    publishDir "$params.output.folder10/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq10000000000
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads11 = params.input.fastq_path11+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads11,size:1)
params.genome11 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt11 {
    publishDir "$params.output.folder11/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out25937424601

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		"""
        




}

process maptoreference11 {
    publishDir "$params.output.folder11/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder11/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out25937424601

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out25937424601

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble11 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder11/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out25937424601


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout25937424601
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate11 {

    publishDir "$params.output.folder11/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout25937424601
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq25937424601

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern11 {

    publishDir "$params.output.folder11/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq25937424601
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads12 = params.input.fastq_path12+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads12,size:1)
params.genome12 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt12 {
    publishDir "$params.output.folder12/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out61917364224

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		"""
        




}

process maptoreference12 {
    publishDir "$params.output.folder12/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder12/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out61917364224

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out61917364224

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble12 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder12/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out61917364224


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout61917364224
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate12 {

    publishDir "$params.output.folder12/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout61917364224
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq61917364224

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern12 {

    publishDir "$params.output.folder12/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq61917364224
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads13 = params.input.fastq_path13+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads13,size:1)
params.genome13 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt13 {
    publishDir "$params.output.folder13/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out137858491849

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		"""
        




}

process maptoreference13 {
    publishDir "$params.output.folder13/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder13/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out137858491849

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out137858491849

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble13 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder13/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out137858491849


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout137858491849
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate13 {

    publishDir "$params.output.folder13/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout137858491849
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq137858491849

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern13 {

    publishDir "$params.output.folder13/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq137858491849
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads14 = params.input.fastq_path14+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads14,size:3)
params.genome14 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt14 {
    publishDir "$params.output.folder14/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out289254654976
		set val(id),  path("${id}_trimmed_1.fq") into trim_out578509309952
		set val(id),  path("${id}_trimmed_2.fq") into trim_out867763964928

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
		NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq
		"""
        




}

process maptoreference14 {
    publishDir "$params.output.folder14/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder14/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out289254654976
		set val(id), path(trim_read1) from trim_out578509309952
		set val(id), path(trim_read2) from trim_out867763964928

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out289254654976
		set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out578509309952
		set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out867763964928

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble14 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder14/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out289254654976
		set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out578509309952
		set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out867763964928


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout289254654976
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq  -s ${sample}_mapped_1.fq  -s ${sample}_mapped_2.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate14 {

    publishDir "$params.output.folder14/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout289254654976
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq289254654976

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern14 {

    publishDir "$params.output.folder14/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq289254654976
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads15 = params.input.fastq_path15+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads15,size:3)
params.genome15 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt15 {
    publishDir "$params.output.folder15/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out576650390625
		set val(id),  path("${id}_trimmed_1.fq") into trim_out1153300781250
		set val(id),  path("${id}_trimmed_2.fq") into trim_out1729951171875

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
		NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq
		"""
        




}

process maptoreference15 {
    publishDir "$params.output.folder15/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder15/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out576650390625
		set val(id), path(trim_read1) from trim_out1153300781250
		set val(id), path(trim_read2) from trim_out1729951171875

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out576650390625
		set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out1153300781250
		set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out1729951171875

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble15 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder15/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out576650390625
		set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out1153300781250
		set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out1729951171875


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout576650390625
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq  -s ${sample}_mapped_1.fq  -s ${sample}_mapped_2.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate15 {

    publishDir "$params.output.folder15/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout576650390625
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq576650390625

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern15 {

    publishDir "$params.output.folder15/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq576650390625
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads16 = params.input.fastq_path16+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads16,size:3)
params.genome16 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt16 {
    publishDir "$params.output.folder16/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out1099511627776
		set val(id),  path("${id}_trimmed_1.fq") into trim_out2199023255552
		set val(id),  path("${id}_trimmed_2.fq") into trim_out3298534883328

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
		NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq
		"""
        




}

process maptoreference16 {
    publishDir "$params.output.folder16/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder16/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out1099511627776
		set val(id), path(trim_read1) from trim_out2199023255552
		set val(id), path(trim_read2) from trim_out3298534883328

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1099511627776
		set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2199023255552
		set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3298534883328

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble16 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder16/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1099511627776
		set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2199023255552
		set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3298534883328


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1099511627776
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq  -s ${sample}_mapped_1.fq  -s ${sample}_mapped_2.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate16 {

    publishDir "$params.output.folder16/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1099511627776
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1099511627776

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern16 {

    publishDir "$params.output.folder16/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1099511627776
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads17 = params.input.fastq_path17+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads17,size:1)
params.genome17 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt17 {
    publishDir "$params.output.folder17/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out2015993900449

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		"""
        




}

process maptoreference17 {
    publishDir "$params.output.folder17/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder17/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out2015993900449

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out2015993900449

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble17 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder17/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out2015993900449


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout2015993900449
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate17 {

    publishDir "$params.output.folder17/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout2015993900449
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq2015993900449

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern17 {

    publishDir "$params.output.folder17/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq2015993900449
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads18 = params.input.fastq_path18+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads18,size:1)
params.genome18 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt18 {
    publishDir "$params.output.folder18/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out3570467226624

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		"""
        




}

process maptoreference18 {
    publishDir "$params.output.folder18/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder18/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out3570467226624

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out3570467226624

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble18 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder18/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out3570467226624


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout3570467226624
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate18 {

    publishDir "$params.output.folder18/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout3570467226624
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq3570467226624

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern18 {

    publishDir "$params.output.folder18/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq3570467226624
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads19 = params.input.fastq_path19+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads19,size:3)
params.genome19 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt19 {
    publishDir "$params.output.folder19/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out6131066257801
		set val(id),  path("${id}_trimmed_1.fq") into trim_out12262132515602
		set val(id),  path("${id}_trimmed_2.fq") into trim_out18393198773403

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
		NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq
		"""
        




}

process maptoreference19 {
    publishDir "$params.output.folder19/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder19/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out6131066257801
		set val(id), path(trim_read1) from trim_out12262132515602
		set val(id), path(trim_read2) from trim_out18393198773403

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out6131066257801
		set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out12262132515602
		set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out18393198773403

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble19 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder19/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out6131066257801
		set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out12262132515602
		set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out18393198773403


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout6131066257801
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq  -s ${sample}_mapped_1.fq  -s ${sample}_mapped_2.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate19 {

    publishDir "$params.output.folder19/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout6131066257801
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq6131066257801

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern19 {

    publishDir "$params.output.folder19/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq6131066257801
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads20 = params.input.fastq_path20+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads20,size:40)
params.genome20 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt20 {
    publishDir "$params.output.folder20/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out10240000000000
		set val(id),  path("${id}_trimmed_1.fq") into trim_out20480000000000
		set val(id),  path("${id}_trimmed_2.fq") into trim_out30720000000000
		set val(id),  path("${id}_trimmed_3.fq") into trim_out40960000000000
		set val(id),  path("${id}_trimmed_4.fq") into trim_out51200000000000
		set val(id),  path("${id}_trimmed_5.fq") into trim_out61440000000000
		set val(id),  path("${id}_trimmed_6.fq") into trim_out71680000000000
		set val(id),  path("${id}_trimmed_7.fq") into trim_out81920000000000
		set val(id),  path("${id}_trimmed_8.fq") into trim_out92160000000000
		set val(id),  path("${id}_trimmed_9.fq") into trim_out102400000000000
		set val(id),  path("${id}_trimmed_10.fq") into trim_out112640000000000
		set val(id),  path("${id}_trimmed_11.fq") into trim_out122880000000000
		set val(id),  path("${id}_trimmed_12.fq") into trim_out133120000000000
		set val(id),  path("${id}_trimmed_13.fq") into trim_out143360000000000
		set val(id),  path("${id}_trimmed_14.fq") into trim_out153600000000000
		set val(id),  path("${id}_trimmed_15.fq") into trim_out163840000000000
		set val(id),  path("${id}_trimmed_16.fq") into trim_out174080000000000
		set val(id),  path("${id}_trimmed_17.fq") into trim_out184320000000000
		set val(id),  path("${id}_trimmed_18.fq") into trim_out194560000000000
		set val(id),  path("${id}_trimmed_19.fq") into trim_out204800000000000
		set val(id),  path("${id}_trimmed_20.fq") into trim_out215040000000000
		set val(id),  path("${id}_trimmed_21.fq") into trim_out225280000000000
		set val(id),  path("${id}_trimmed_22.fq") into trim_out235520000000000
		set val(id),  path("${id}_trimmed_23.fq") into trim_out245760000000000
		set val(id),  path("${id}_trimmed_24.fq") into trim_out256000000000000
		set val(id),  path("${id}_trimmed_25.fq") into trim_out266240000000000
		set val(id),  path("${id}_trimmed_26.fq") into trim_out276480000000000
		set val(id),  path("${id}_trimmed_27.fq") into trim_out286720000000000
		set val(id),  path("${id}_trimmed_28.fq") into trim_out296960000000000
		set val(id),  path("${id}_trimmed_29.fq") into trim_out307200000000000
		set val(id),  path("${id}_trimmed_30.fq") into trim_out317440000000000
		set val(id),  path("${id}_trimmed_31.fq") into trim_out327680000000000
		set val(id),  path("${id}_trimmed_32.fq") into trim_out337920000000000
		set val(id),  path("${id}_trimmed_33.fq") into trim_out348160000000000
		set val(id),  path("${id}_trimmed_34.fq") into trim_out358400000000000
		set val(id),  path("${id}_trimmed_35.fq") into trim_out368640000000000
		set val(id),  path("${id}_trimmed_36.fq") into trim_out378880000000000
		set val(id),  path("${id}_trimmed_37.fq") into trim_out389120000000000
		set val(id),  path("${id}_trimmed_38.fq") into trim_out399360000000000
		set val(id),  path("${id}_trimmed_39.fq") into trim_out409600000000000

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
		NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq
		NanoFilt ${id}_3.fastq -l 500 -q 10 > ${id}_trimmed_3.fq
		NanoFilt ${id}_4.fastq -l 500 -q 10 > ${id}_trimmed_4.fq
		NanoFilt ${id}_5.fastq -l 500 -q 10 > ${id}_trimmed_5.fq
		NanoFilt ${id}_6.fastq -l 500 -q 10 > ${id}_trimmed_6.fq
		NanoFilt ${id}_7.fastq -l 500 -q 10 > ${id}_trimmed_7.fq
		NanoFilt ${id}_8.fastq -l 500 -q 10 > ${id}_trimmed_8.fq
		NanoFilt ${id}_9.fastq -l 500 -q 10 > ${id}_trimmed_9.fq
		NanoFilt ${id}_10.fastq -l 500 -q 10 > ${id}_trimmed_10.fq
		NanoFilt ${id}_11.fastq -l 500 -q 10 > ${id}_trimmed_11.fq
		NanoFilt ${id}_12.fastq -l 500 -q 10 > ${id}_trimmed_12.fq
		NanoFilt ${id}_13.fastq -l 500 -q 10 > ${id}_trimmed_13.fq
		NanoFilt ${id}_14.fastq -l 500 -q 10 > ${id}_trimmed_14.fq
		NanoFilt ${id}_15.fastq -l 500 -q 10 > ${id}_trimmed_15.fq
		NanoFilt ${id}_16.fastq -l 500 -q 10 > ${id}_trimmed_16.fq
		NanoFilt ${id}_17.fastq -l 500 -q 10 > ${id}_trimmed_17.fq
		NanoFilt ${id}_18.fastq -l 500 -q 10 > ${id}_trimmed_18.fq
		NanoFilt ${id}_19.fastq -l 500 -q 10 > ${id}_trimmed_19.fq
		NanoFilt ${id}_20.fastq -l 500 -q 10 > ${id}_trimmed_20.fq
		NanoFilt ${id}_21.fastq -l 500 -q 10 > ${id}_trimmed_21.fq
		NanoFilt ${id}_22.fastq -l 500 -q 10 > ${id}_trimmed_22.fq
		NanoFilt ${id}_23.fastq -l 500 -q 10 > ${id}_trimmed_23.fq
		NanoFilt ${id}_24.fastq -l 500 -q 10 > ${id}_trimmed_24.fq
		NanoFilt ${id}_25.fastq -l 500 -q 10 > ${id}_trimmed_25.fq
		NanoFilt ${id}_26.fastq -l 500 -q 10 > ${id}_trimmed_26.fq
		NanoFilt ${id}_27.fastq -l 500 -q 10 > ${id}_trimmed_27.fq
		NanoFilt ${id}_28.fastq -l 500 -q 10 > ${id}_trimmed_28.fq
		NanoFilt ${id}_29.fastq -l 500 -q 10 > ${id}_trimmed_29.fq
		NanoFilt ${id}_30.fastq -l 500 -q 10 > ${id}_trimmed_30.fq
		NanoFilt ${id}_31.fastq -l 500 -q 10 > ${id}_trimmed_31.fq
		NanoFilt ${id}_32.fastq -l 500 -q 10 > ${id}_trimmed_32.fq
		NanoFilt ${id}_33.fastq -l 500 -q 10 > ${id}_trimmed_33.fq
		NanoFilt ${id}_34.fastq -l 500 -q 10 > ${id}_trimmed_34.fq
		NanoFilt ${id}_35.fastq -l 500 -q 10 > ${id}_trimmed_35.fq
		NanoFilt ${id}_36.fastq -l 500 -q 10 > ${id}_trimmed_36.fq
		NanoFilt ${id}_37.fastq -l 500 -q 10 > ${id}_trimmed_37.fq
		NanoFilt ${id}_38.fastq -l 500 -q 10 > ${id}_trimmed_38.fq
		NanoFilt ${id}_39.fastq -l 500 -q 10 > ${id}_trimmed_39.fq
		"""
        




}

process maptoreference20 {
    publishDir "$params.output.folder20/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder20/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out10240000000000
		set val(id), path(trim_read1) from trim_out20480000000000
		set val(id), path(trim_read2) from trim_out30720000000000
		set val(id), path(trim_read3) from trim_out40960000000000
		set val(id), path(trim_read4) from trim_out51200000000000
		set val(id), path(trim_read5) from trim_out61440000000000
		set val(id), path(trim_read6) from trim_out71680000000000
		set val(id), path(trim_read7) from trim_out81920000000000
		set val(id), path(trim_read8) from trim_out92160000000000
		set val(id), path(trim_read9) from trim_out102400000000000
		set val(id), path(trim_read10) from trim_out112640000000000
		set val(id), path(trim_read11) from trim_out122880000000000
		set val(id), path(trim_read12) from trim_out133120000000000
		set val(id), path(trim_read13) from trim_out143360000000000
		set val(id), path(trim_read14) from trim_out153600000000000
		set val(id), path(trim_read15) from trim_out163840000000000
		set val(id), path(trim_read16) from trim_out174080000000000
		set val(id), path(trim_read17) from trim_out184320000000000
		set val(id), path(trim_read18) from trim_out194560000000000
		set val(id), path(trim_read19) from trim_out204800000000000
		set val(id), path(trim_read20) from trim_out215040000000000
		set val(id), path(trim_read21) from trim_out225280000000000
		set val(id), path(trim_read22) from trim_out235520000000000
		set val(id), path(trim_read23) from trim_out245760000000000
		set val(id), path(trim_read24) from trim_out256000000000000
		set val(id), path(trim_read25) from trim_out266240000000000
		set val(id), path(trim_read26) from trim_out276480000000000
		set val(id), path(trim_read27) from trim_out286720000000000
		set val(id), path(trim_read28) from trim_out296960000000000
		set val(id), path(trim_read29) from trim_out307200000000000
		set val(id), path(trim_read30) from trim_out317440000000000
		set val(id), path(trim_read31) from trim_out327680000000000
		set val(id), path(trim_read32) from trim_out337920000000000
		set val(id), path(trim_read33) from trim_out348160000000000
		set val(id), path(trim_read34) from trim_out358400000000000
		set val(id), path(trim_read35) from trim_out368640000000000
		set val(id), path(trim_read36) from trim_out378880000000000
		set val(id), path(trim_read37) from trim_out389120000000000
		set val(id), path(trim_read38) from trim_out399360000000000
		set val(id), path(trim_read39) from trim_out409600000000000

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out10240000000000
		set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out20480000000000
		set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out30720000000000
		set val(id), path("${id}_mapped_3.fq"), path("${id}_unmapped_3.fq") into mapped_out40960000000000
		set val(id), path("${id}_mapped_4.fq"), path("${id}_unmapped_4.fq") into mapped_out51200000000000
		set val(id), path("${id}_mapped_5.fq"), path("${id}_unmapped_5.fq") into mapped_out61440000000000
		set val(id), path("${id}_mapped_6.fq"), path("${id}_unmapped_6.fq") into mapped_out71680000000000
		set val(id), path("${id}_mapped_7.fq"), path("${id}_unmapped_7.fq") into mapped_out81920000000000
		set val(id), path("${id}_mapped_8.fq"), path("${id}_unmapped_8.fq") into mapped_out92160000000000
		set val(id), path("${id}_mapped_9.fq"), path("${id}_unmapped_9.fq") into mapped_out102400000000000
		set val(id), path("${id}_mapped_10.fq"), path("${id}_unmapped_10.fq") into mapped_out112640000000000
		set val(id), path("${id}_mapped_11.fq"), path("${id}_unmapped_11.fq") into mapped_out122880000000000
		set val(id), path("${id}_mapped_12.fq"), path("${id}_unmapped_12.fq") into mapped_out133120000000000
		set val(id), path("${id}_mapped_13.fq"), path("${id}_unmapped_13.fq") into mapped_out143360000000000
		set val(id), path("${id}_mapped_14.fq"), path("${id}_unmapped_14.fq") into mapped_out153600000000000
		set val(id), path("${id}_mapped_15.fq"), path("${id}_unmapped_15.fq") into mapped_out163840000000000
		set val(id), path("${id}_mapped_16.fq"), path("${id}_unmapped_16.fq") into mapped_out174080000000000
		set val(id), path("${id}_mapped_17.fq"), path("${id}_unmapped_17.fq") into mapped_out184320000000000
		set val(id), path("${id}_mapped_18.fq"), path("${id}_unmapped_18.fq") into mapped_out194560000000000
		set val(id), path("${id}_mapped_19.fq"), path("${id}_unmapped_19.fq") into mapped_out204800000000000
		set val(id), path("${id}_mapped_20.fq"), path("${id}_unmapped_20.fq") into mapped_out215040000000000
		set val(id), path("${id}_mapped_21.fq"), path("${id}_unmapped_21.fq") into mapped_out225280000000000
		set val(id), path("${id}_mapped_22.fq"), path("${id}_unmapped_22.fq") into mapped_out235520000000000
		set val(id), path("${id}_mapped_23.fq"), path("${id}_unmapped_23.fq") into mapped_out245760000000000
		set val(id), path("${id}_mapped_24.fq"), path("${id}_unmapped_24.fq") into mapped_out256000000000000
		set val(id), path("${id}_mapped_25.fq"), path("${id}_unmapped_25.fq") into mapped_out266240000000000
		set val(id), path("${id}_mapped_26.fq"), path("${id}_unmapped_26.fq") into mapped_out276480000000000
		set val(id), path("${id}_mapped_27.fq"), path("${id}_unmapped_27.fq") into mapped_out286720000000000
		set val(id), path("${id}_mapped_28.fq"), path("${id}_unmapped_28.fq") into mapped_out296960000000000
		set val(id), path("${id}_mapped_29.fq"), path("${id}_unmapped_29.fq") into mapped_out307200000000000
		set val(id), path("${id}_mapped_30.fq"), path("${id}_unmapped_30.fq") into mapped_out317440000000000
		set val(id), path("${id}_mapped_31.fq"), path("${id}_unmapped_31.fq") into mapped_out327680000000000
		set val(id), path("${id}_mapped_32.fq"), path("${id}_unmapped_32.fq") into mapped_out337920000000000
		set val(id), path("${id}_mapped_33.fq"), path("${id}_unmapped_33.fq") into mapped_out348160000000000
		set val(id), path("${id}_mapped_34.fq"), path("${id}_unmapped_34.fq") into mapped_out358400000000000
		set val(id), path("${id}_mapped_35.fq"), path("${id}_unmapped_35.fq") into mapped_out368640000000000
		set val(id), path("${id}_mapped_36.fq"), path("${id}_unmapped_36.fq") into mapped_out378880000000000
		set val(id), path("${id}_mapped_37.fq"), path("${id}_unmapped_37.fq") into mapped_out389120000000000
		set val(id), path("${id}_mapped_38.fq"), path("${id}_unmapped_38.fq") into mapped_out399360000000000
		set val(id), path("${id}_mapped_39.fq"), path("${id}_unmapped_39.fq") into mapped_out409600000000000

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_3.fq outm=${id}_mapped_3.fq outu=${id}_unmapped_3.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_4.fq outm=${id}_mapped_4.fq outu=${id}_unmapped_4.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_5.fq outm=${id}_mapped_5.fq outu=${id}_unmapped_5.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_6.fq outm=${id}_mapped_6.fq outu=${id}_unmapped_6.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_7.fq outm=${id}_mapped_7.fq outu=${id}_unmapped_7.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_8.fq outm=${id}_mapped_8.fq outu=${id}_unmapped_8.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_9.fq outm=${id}_mapped_9.fq outu=${id}_unmapped_9.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_10.fq outm=${id}_mapped_10.fq outu=${id}_unmapped_10.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_11.fq outm=${id}_mapped_11.fq outu=${id}_unmapped_11.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_12.fq outm=${id}_mapped_12.fq outu=${id}_unmapped_12.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_13.fq outm=${id}_mapped_13.fq outu=${id}_unmapped_13.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_14.fq outm=${id}_mapped_14.fq outu=${id}_unmapped_14.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_15.fq outm=${id}_mapped_15.fq outu=${id}_unmapped_15.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_16.fq outm=${id}_mapped_16.fq outu=${id}_unmapped_16.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_17.fq outm=${id}_mapped_17.fq outu=${id}_unmapped_17.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_18.fq outm=${id}_mapped_18.fq outu=${id}_unmapped_18.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_19.fq outm=${id}_mapped_19.fq outu=${id}_unmapped_19.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_20.fq outm=${id}_mapped_20.fq outu=${id}_unmapped_20.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_21.fq outm=${id}_mapped_21.fq outu=${id}_unmapped_21.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_22.fq outm=${id}_mapped_22.fq outu=${id}_unmapped_22.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_23.fq outm=${id}_mapped_23.fq outu=${id}_unmapped_23.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_24.fq outm=${id}_mapped_24.fq outu=${id}_unmapped_24.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_25.fq outm=${id}_mapped_25.fq outu=${id}_unmapped_25.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_26.fq outm=${id}_mapped_26.fq outu=${id}_unmapped_26.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_27.fq outm=${id}_mapped_27.fq outu=${id}_unmapped_27.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_28.fq outm=${id}_mapped_28.fq outu=${id}_unmapped_28.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_29.fq outm=${id}_mapped_29.fq outu=${id}_unmapped_29.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_30.fq outm=${id}_mapped_30.fq outu=${id}_unmapped_30.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_31.fq outm=${id}_mapped_31.fq outu=${id}_unmapped_31.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_32.fq outm=${id}_mapped_32.fq outu=${id}_unmapped_32.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_33.fq outm=${id}_mapped_33.fq outu=${id}_unmapped_33.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_34.fq outm=${id}_mapped_34.fq outu=${id}_unmapped_34.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_35.fq outm=${id}_mapped_35.fq outu=${id}_unmapped_35.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_36.fq outm=${id}_mapped_36.fq outu=${id}_unmapped_36.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_37.fq outm=${id}_mapped_37.fq outu=${id}_unmapped_37.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_38.fq outm=${id}_mapped_38.fq outu=${id}_unmapped_38.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_39.fq outm=${id}_mapped_39.fq outu=${id}_unmapped_39.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble20 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder20/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out10240000000000
		set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out20480000000000
		set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out30720000000000
		set val(sample), path(mapped_read3), path(unmapped_read3) from mapped_out40960000000000
		set val(sample), path(mapped_read4), path(unmapped_read4) from mapped_out51200000000000
		set val(sample), path(mapped_read5), path(unmapped_read5) from mapped_out61440000000000
		set val(sample), path(mapped_read6), path(unmapped_read6) from mapped_out71680000000000
		set val(sample), path(mapped_read7), path(unmapped_read7) from mapped_out81920000000000
		set val(sample), path(mapped_read8), path(unmapped_read8) from mapped_out92160000000000
		set val(sample), path(mapped_read9), path(unmapped_read9) from mapped_out102400000000000
		set val(sample), path(mapped_read10), path(unmapped_read10) from mapped_out112640000000000
		set val(sample), path(mapped_read11), path(unmapped_read11) from mapped_out122880000000000
		set val(sample), path(mapped_read12), path(unmapped_read12) from mapped_out133120000000000
		set val(sample), path(mapped_read13), path(unmapped_read13) from mapped_out143360000000000
		set val(sample), path(mapped_read14), path(unmapped_read14) from mapped_out153600000000000
		set val(sample), path(mapped_read15), path(unmapped_read15) from mapped_out163840000000000
		set val(sample), path(mapped_read16), path(unmapped_read16) from mapped_out174080000000000
		set val(sample), path(mapped_read17), path(unmapped_read17) from mapped_out184320000000000
		set val(sample), path(mapped_read18), path(unmapped_read18) from mapped_out194560000000000
		set val(sample), path(mapped_read19), path(unmapped_read19) from mapped_out204800000000000
		set val(sample), path(mapped_read20), path(unmapped_read20) from mapped_out215040000000000
		set val(sample), path(mapped_read21), path(unmapped_read21) from mapped_out225280000000000
		set val(sample), path(mapped_read22), path(unmapped_read22) from mapped_out235520000000000
		set val(sample), path(mapped_read23), path(unmapped_read23) from mapped_out245760000000000
		set val(sample), path(mapped_read24), path(unmapped_read24) from mapped_out256000000000000
		set val(sample), path(mapped_read25), path(unmapped_read25) from mapped_out266240000000000
		set val(sample), path(mapped_read26), path(unmapped_read26) from mapped_out276480000000000
		set val(sample), path(mapped_read27), path(unmapped_read27) from mapped_out286720000000000
		set val(sample), path(mapped_read28), path(unmapped_read28) from mapped_out296960000000000
		set val(sample), path(mapped_read29), path(unmapped_read29) from mapped_out307200000000000
		set val(sample), path(mapped_read30), path(unmapped_read30) from mapped_out317440000000000
		set val(sample), path(mapped_read31), path(unmapped_read31) from mapped_out327680000000000
		set val(sample), path(mapped_read32), path(unmapped_read32) from mapped_out337920000000000
		set val(sample), path(mapped_read33), path(unmapped_read33) from mapped_out348160000000000
		set val(sample), path(mapped_read34), path(unmapped_read34) from mapped_out358400000000000
		set val(sample), path(mapped_read35), path(unmapped_read35) from mapped_out368640000000000
		set val(sample), path(mapped_read36), path(unmapped_read36) from mapped_out378880000000000
		set val(sample), path(mapped_read37), path(unmapped_read37) from mapped_out389120000000000
		set val(sample), path(mapped_read38), path(unmapped_read38) from mapped_out399360000000000
		set val(sample), path(mapped_read39), path(unmapped_read39) from mapped_out409600000000000


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout10240000000000
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq  -s ${sample}_mapped_1.fq  -s ${sample}_mapped_2.fq  -s ${sample}_mapped_3.fq  -s ${sample}_mapped_4.fq  -s ${sample}_mapped_5.fq  -s ${sample}_mapped_6.fq  -s ${sample}_mapped_7.fq  -s ${sample}_mapped_8.fq  -s ${sample}_mapped_9.fq  -s ${sample}_mapped_10.fq  -s ${sample}_mapped_11.fq  -s ${sample}_mapped_12.fq  -s ${sample}_mapped_13.fq  -s ${sample}_mapped_14.fq  -s ${sample}_mapped_15.fq  -s ${sample}_mapped_16.fq  -s ${sample}_mapped_17.fq  -s ${sample}_mapped_18.fq  -s ${sample}_mapped_19.fq  -s ${sample}_mapped_20.fq  -s ${sample}_mapped_21.fq  -s ${sample}_mapped_22.fq  -s ${sample}_mapped_23.fq  -s ${sample}_mapped_24.fq  -s ${sample}_mapped_25.fq  -s ${sample}_mapped_26.fq  -s ${sample}_mapped_27.fq  -s ${sample}_mapped_28.fq  -s ${sample}_mapped_29.fq  -s ${sample}_mapped_30.fq  -s ${sample}_mapped_31.fq  -s ${sample}_mapped_32.fq  -s ${sample}_mapped_33.fq  -s ${sample}_mapped_34.fq  -s ${sample}_mapped_35.fq  -s ${sample}_mapped_36.fq  -s ${sample}_mapped_37.fq  -s ${sample}_mapped_38.fq  -s ${sample}_mapped_39.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate20 {

    publishDir "$params.output.folder20/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout10240000000000
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq10240000000000

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern20 {

    publishDir "$params.output.folder20/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq10240000000000
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads21 = params.input.fastq_path21+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads21,size:1)
params.genome21 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt21 {
    publishDir "$params.output.folder21/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out16679880978201

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		"""
        




}

process maptoreference21 {
    publishDir "$params.output.folder21/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder21/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out16679880978201

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out16679880978201

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble21 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder21/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out16679880978201


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout16679880978201
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate21 {

    publishDir "$params.output.folder21/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout16679880978201
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq16679880978201

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern21 {

    publishDir "$params.output.folder21/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq16679880978201
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads22 = params.input.fastq_path22+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads22,size:1)
params.genome22 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt22 {
    publishDir "$params.output.folder22/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out26559922791424

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		"""
        




}

process maptoreference22 {
    publishDir "$params.output.folder22/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder22/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out26559922791424

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out26559922791424

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble22 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder22/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out26559922791424


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout26559922791424
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate22 {

    publishDir "$params.output.folder22/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout26559922791424
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq26559922791424

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern22 {

    publishDir "$params.output.folder22/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq26559922791424
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads23 = params.input.fastq_path23+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads23,size:1)
params.genome23 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt23 {
    publishDir "$params.output.folder23/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out41426511213649

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		"""
        




}

process maptoreference23 {
    publishDir "$params.output.folder23/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder23/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out41426511213649

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out41426511213649

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble23 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder23/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out41426511213649


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout41426511213649
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate23 {

    publishDir "$params.output.folder23/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout41426511213649
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq41426511213649

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern23 {

    publishDir "$params.output.folder23/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq41426511213649
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads24 = params.input.fastq_path24+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads24,size:1)
params.genome24 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt24 {
    publishDir "$params.output.folder24/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out63403380965376

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		"""
        




}

process maptoreference24 {
    publishDir "$params.output.folder24/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder24/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out63403380965376

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out63403380965376

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble24 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder24/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out63403380965376


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout63403380965376
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate24 {

    publishDir "$params.output.folder24/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout63403380965376
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq63403380965376

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern24 {

    publishDir "$params.output.folder24/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq63403380965376
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads25 = params.input.fastq_path25+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads25,size:1)
params.genome25 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt25 {
    publishDir "$params.output.folder25/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out95367431640625

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		"""
        




}

process maptoreference25 {
    publishDir "$params.output.folder25/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder25/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out95367431640625

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out95367431640625

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble25 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder25/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out95367431640625


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout95367431640625
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate25 {

    publishDir "$params.output.folder25/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout95367431640625
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq95367431640625

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern25 {

    publishDir "$params.output.folder25/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq95367431640625
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads26 = params.input.fastq_path26+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads26,size:31)
params.genome26 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt26 {
    publishDir "$params.output.folder26/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out141167095653376
		set val(id),  path("${id}_trimmed_1.fq") into trim_out282334191306752
		set val(id),  path("${id}_trimmed_2.fq") into trim_out423501286960128
		set val(id),  path("${id}_trimmed_3.fq") into trim_out564668382613504
		set val(id),  path("${id}_trimmed_4.fq") into trim_out705835478266880
		set val(id),  path("${id}_trimmed_5.fq") into trim_out847002573920256
		set val(id),  path("${id}_trimmed_6.fq") into trim_out988169669573632
		set val(id),  path("${id}_trimmed_7.fq") into trim_out1129336765227008
		set val(id),  path("${id}_trimmed_8.fq") into trim_out1270503860880384
		set val(id),  path("${id}_trimmed_9.fq") into trim_out1411670956533760
		set val(id),  path("${id}_trimmed_10.fq") into trim_out1552838052187136
		set val(id),  path("${id}_trimmed_11.fq") into trim_out1694005147840512
		set val(id),  path("${id}_trimmed_12.fq") into trim_out1835172243493888
		set val(id),  path("${id}_trimmed_13.fq") into trim_out1976339339147264
		set val(id),  path("${id}_trimmed_14.fq") into trim_out2117506434800640
		set val(id),  path("${id}_trimmed_15.fq") into trim_out2258673530454016
		set val(id),  path("${id}_trimmed_16.fq") into trim_out2399840626107392
		set val(id),  path("${id}_trimmed_17.fq") into trim_out2541007721760768
		set val(id),  path("${id}_trimmed_18.fq") into trim_out2682174817414144
		set val(id),  path("${id}_trimmed_19.fq") into trim_out2823341913067520
		set val(id),  path("${id}_trimmed_20.fq") into trim_out2964509008720896
		set val(id),  path("${id}_trimmed_21.fq") into trim_out3105676104374272
		set val(id),  path("${id}_trimmed_22.fq") into trim_out3246843200027648
		set val(id),  path("${id}_trimmed_23.fq") into trim_out3388010295681024
		set val(id),  path("${id}_trimmed_24.fq") into trim_out3529177391334400
		set val(id),  path("${id}_trimmed_25.fq") into trim_out3670344486987776
		set val(id),  path("${id}_trimmed_26.fq") into trim_out3811511582641152
		set val(id),  path("${id}_trimmed_27.fq") into trim_out3952678678294528
		set val(id),  path("${id}_trimmed_28.fq") into trim_out4093845773947904
		set val(id),  path("${id}_trimmed_29.fq") into trim_out4235012869601280
		set val(id),  path("${id}_trimmed_30.fq") into trim_out4376179965254656

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
		NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq
		NanoFilt ${id}_3.fastq -l 500 -q 10 > ${id}_trimmed_3.fq
		NanoFilt ${id}_4.fastq -l 500 -q 10 > ${id}_trimmed_4.fq
		NanoFilt ${id}_5.fastq -l 500 -q 10 > ${id}_trimmed_5.fq
		NanoFilt ${id}_6.fastq -l 500 -q 10 > ${id}_trimmed_6.fq
		NanoFilt ${id}_7.fastq -l 500 -q 10 > ${id}_trimmed_7.fq
		NanoFilt ${id}_8.fastq -l 500 -q 10 > ${id}_trimmed_8.fq
		NanoFilt ${id}_9.fastq -l 500 -q 10 > ${id}_trimmed_9.fq
		NanoFilt ${id}_10.fastq -l 500 -q 10 > ${id}_trimmed_10.fq
		NanoFilt ${id}_11.fastq -l 500 -q 10 > ${id}_trimmed_11.fq
		NanoFilt ${id}_12.fastq -l 500 -q 10 > ${id}_trimmed_12.fq
		NanoFilt ${id}_13.fastq -l 500 -q 10 > ${id}_trimmed_13.fq
		NanoFilt ${id}_14.fastq -l 500 -q 10 > ${id}_trimmed_14.fq
		NanoFilt ${id}_15.fastq -l 500 -q 10 > ${id}_trimmed_15.fq
		NanoFilt ${id}_16.fastq -l 500 -q 10 > ${id}_trimmed_16.fq
		NanoFilt ${id}_17.fastq -l 500 -q 10 > ${id}_trimmed_17.fq
		NanoFilt ${id}_18.fastq -l 500 -q 10 > ${id}_trimmed_18.fq
		NanoFilt ${id}_19.fastq -l 500 -q 10 > ${id}_trimmed_19.fq
		NanoFilt ${id}_20.fastq -l 500 -q 10 > ${id}_trimmed_20.fq
		NanoFilt ${id}_21.fastq -l 500 -q 10 > ${id}_trimmed_21.fq
		NanoFilt ${id}_22.fastq -l 500 -q 10 > ${id}_trimmed_22.fq
		NanoFilt ${id}_23.fastq -l 500 -q 10 > ${id}_trimmed_23.fq
		NanoFilt ${id}_24.fastq -l 500 -q 10 > ${id}_trimmed_24.fq
		NanoFilt ${id}_25.fastq -l 500 -q 10 > ${id}_trimmed_25.fq
		NanoFilt ${id}_26.fastq -l 500 -q 10 > ${id}_trimmed_26.fq
		NanoFilt ${id}_27.fastq -l 500 -q 10 > ${id}_trimmed_27.fq
		NanoFilt ${id}_28.fastq -l 500 -q 10 > ${id}_trimmed_28.fq
		NanoFilt ${id}_29.fastq -l 500 -q 10 > ${id}_trimmed_29.fq
		NanoFilt ${id}_30.fastq -l 500 -q 10 > ${id}_trimmed_30.fq
		"""
        




}

process maptoreference26 {
    publishDir "$params.output.folder26/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder26/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out141167095653376
		set val(id), path(trim_read1) from trim_out282334191306752
		set val(id), path(trim_read2) from trim_out423501286960128
		set val(id), path(trim_read3) from trim_out564668382613504
		set val(id), path(trim_read4) from trim_out705835478266880
		set val(id), path(trim_read5) from trim_out847002573920256
		set val(id), path(trim_read6) from trim_out988169669573632
		set val(id), path(trim_read7) from trim_out1129336765227008
		set val(id), path(trim_read8) from trim_out1270503860880384
		set val(id), path(trim_read9) from trim_out1411670956533760
		set val(id), path(trim_read10) from trim_out1552838052187136
		set val(id), path(trim_read11) from trim_out1694005147840512
		set val(id), path(trim_read12) from trim_out1835172243493888
		set val(id), path(trim_read13) from trim_out1976339339147264
		set val(id), path(trim_read14) from trim_out2117506434800640
		set val(id), path(trim_read15) from trim_out2258673530454016
		set val(id), path(trim_read16) from trim_out2399840626107392
		set val(id), path(trim_read17) from trim_out2541007721760768
		set val(id), path(trim_read18) from trim_out2682174817414144
		set val(id), path(trim_read19) from trim_out2823341913067520
		set val(id), path(trim_read20) from trim_out2964509008720896
		set val(id), path(trim_read21) from trim_out3105676104374272
		set val(id), path(trim_read22) from trim_out3246843200027648
		set val(id), path(trim_read23) from trim_out3388010295681024
		set val(id), path(trim_read24) from trim_out3529177391334400
		set val(id), path(trim_read25) from trim_out3670344486987776
		set val(id), path(trim_read26) from trim_out3811511582641152
		set val(id), path(trim_read27) from trim_out3952678678294528
		set val(id), path(trim_read28) from trim_out4093845773947904
		set val(id), path(trim_read29) from trim_out4235012869601280
		set val(id), path(trim_read30) from trim_out4376179965254656

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out141167095653376
		set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out282334191306752
		set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out423501286960128
		set val(id), path("${id}_mapped_3.fq"), path("${id}_unmapped_3.fq") into mapped_out564668382613504
		set val(id), path("${id}_mapped_4.fq"), path("${id}_unmapped_4.fq") into mapped_out705835478266880
		set val(id), path("${id}_mapped_5.fq"), path("${id}_unmapped_5.fq") into mapped_out847002573920256
		set val(id), path("${id}_mapped_6.fq"), path("${id}_unmapped_6.fq") into mapped_out988169669573632
		set val(id), path("${id}_mapped_7.fq"), path("${id}_unmapped_7.fq") into mapped_out1129336765227008
		set val(id), path("${id}_mapped_8.fq"), path("${id}_unmapped_8.fq") into mapped_out1270503860880384
		set val(id), path("${id}_mapped_9.fq"), path("${id}_unmapped_9.fq") into mapped_out1411670956533760
		set val(id), path("${id}_mapped_10.fq"), path("${id}_unmapped_10.fq") into mapped_out1552838052187136
		set val(id), path("${id}_mapped_11.fq"), path("${id}_unmapped_11.fq") into mapped_out1694005147840512
		set val(id), path("${id}_mapped_12.fq"), path("${id}_unmapped_12.fq") into mapped_out1835172243493888
		set val(id), path("${id}_mapped_13.fq"), path("${id}_unmapped_13.fq") into mapped_out1976339339147264
		set val(id), path("${id}_mapped_14.fq"), path("${id}_unmapped_14.fq") into mapped_out2117506434800640
		set val(id), path("${id}_mapped_15.fq"), path("${id}_unmapped_15.fq") into mapped_out2258673530454016
		set val(id), path("${id}_mapped_16.fq"), path("${id}_unmapped_16.fq") into mapped_out2399840626107392
		set val(id), path("${id}_mapped_17.fq"), path("${id}_unmapped_17.fq") into mapped_out2541007721760768
		set val(id), path("${id}_mapped_18.fq"), path("${id}_unmapped_18.fq") into mapped_out2682174817414144
		set val(id), path("${id}_mapped_19.fq"), path("${id}_unmapped_19.fq") into mapped_out2823341913067520
		set val(id), path("${id}_mapped_20.fq"), path("${id}_unmapped_20.fq") into mapped_out2964509008720896
		set val(id), path("${id}_mapped_21.fq"), path("${id}_unmapped_21.fq") into mapped_out3105676104374272
		set val(id), path("${id}_mapped_22.fq"), path("${id}_unmapped_22.fq") into mapped_out3246843200027648
		set val(id), path("${id}_mapped_23.fq"), path("${id}_unmapped_23.fq") into mapped_out3388010295681024
		set val(id), path("${id}_mapped_24.fq"), path("${id}_unmapped_24.fq") into mapped_out3529177391334400
		set val(id), path("${id}_mapped_25.fq"), path("${id}_unmapped_25.fq") into mapped_out3670344486987776
		set val(id), path("${id}_mapped_26.fq"), path("${id}_unmapped_26.fq") into mapped_out3811511582641152
		set val(id), path("${id}_mapped_27.fq"), path("${id}_unmapped_27.fq") into mapped_out3952678678294528
		set val(id), path("${id}_mapped_28.fq"), path("${id}_unmapped_28.fq") into mapped_out4093845773947904
		set val(id), path("${id}_mapped_29.fq"), path("${id}_unmapped_29.fq") into mapped_out4235012869601280
		set val(id), path("${id}_mapped_30.fq"), path("${id}_unmapped_30.fq") into mapped_out4376179965254656

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_3.fq outm=${id}_mapped_3.fq outu=${id}_unmapped_3.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_4.fq outm=${id}_mapped_4.fq outu=${id}_unmapped_4.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_5.fq outm=${id}_mapped_5.fq outu=${id}_unmapped_5.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_6.fq outm=${id}_mapped_6.fq outu=${id}_unmapped_6.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_7.fq outm=${id}_mapped_7.fq outu=${id}_unmapped_7.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_8.fq outm=${id}_mapped_8.fq outu=${id}_unmapped_8.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_9.fq outm=${id}_mapped_9.fq outu=${id}_unmapped_9.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_10.fq outm=${id}_mapped_10.fq outu=${id}_unmapped_10.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_11.fq outm=${id}_mapped_11.fq outu=${id}_unmapped_11.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_12.fq outm=${id}_mapped_12.fq outu=${id}_unmapped_12.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_13.fq outm=${id}_mapped_13.fq outu=${id}_unmapped_13.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_14.fq outm=${id}_mapped_14.fq outu=${id}_unmapped_14.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_15.fq outm=${id}_mapped_15.fq outu=${id}_unmapped_15.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_16.fq outm=${id}_mapped_16.fq outu=${id}_unmapped_16.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_17.fq outm=${id}_mapped_17.fq outu=${id}_unmapped_17.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_18.fq outm=${id}_mapped_18.fq outu=${id}_unmapped_18.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_19.fq outm=${id}_mapped_19.fq outu=${id}_unmapped_19.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_20.fq outm=${id}_mapped_20.fq outu=${id}_unmapped_20.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_21.fq outm=${id}_mapped_21.fq outu=${id}_unmapped_21.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_22.fq outm=${id}_mapped_22.fq outu=${id}_unmapped_22.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_23.fq outm=${id}_mapped_23.fq outu=${id}_unmapped_23.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_24.fq outm=${id}_mapped_24.fq outu=${id}_unmapped_24.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_25.fq outm=${id}_mapped_25.fq outu=${id}_unmapped_25.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_26.fq outm=${id}_mapped_26.fq outu=${id}_unmapped_26.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_27.fq outm=${id}_mapped_27.fq outu=${id}_unmapped_27.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_28.fq outm=${id}_mapped_28.fq outu=${id}_unmapped_28.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_29.fq outm=${id}_mapped_29.fq outu=${id}_unmapped_29.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_30.fq outm=${id}_mapped_30.fq outu=${id}_unmapped_30.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble26 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder26/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out141167095653376
		set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out282334191306752
		set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out423501286960128
		set val(sample), path(mapped_read3), path(unmapped_read3) from mapped_out564668382613504
		set val(sample), path(mapped_read4), path(unmapped_read4) from mapped_out705835478266880
		set val(sample), path(mapped_read5), path(unmapped_read5) from mapped_out847002573920256
		set val(sample), path(mapped_read6), path(unmapped_read6) from mapped_out988169669573632
		set val(sample), path(mapped_read7), path(unmapped_read7) from mapped_out1129336765227008
		set val(sample), path(mapped_read8), path(unmapped_read8) from mapped_out1270503860880384
		set val(sample), path(mapped_read9), path(unmapped_read9) from mapped_out1411670956533760
		set val(sample), path(mapped_read10), path(unmapped_read10) from mapped_out1552838052187136
		set val(sample), path(mapped_read11), path(unmapped_read11) from mapped_out1694005147840512
		set val(sample), path(mapped_read12), path(unmapped_read12) from mapped_out1835172243493888
		set val(sample), path(mapped_read13), path(unmapped_read13) from mapped_out1976339339147264
		set val(sample), path(mapped_read14), path(unmapped_read14) from mapped_out2117506434800640
		set val(sample), path(mapped_read15), path(unmapped_read15) from mapped_out2258673530454016
		set val(sample), path(mapped_read16), path(unmapped_read16) from mapped_out2399840626107392
		set val(sample), path(mapped_read17), path(unmapped_read17) from mapped_out2541007721760768
		set val(sample), path(mapped_read18), path(unmapped_read18) from mapped_out2682174817414144
		set val(sample), path(mapped_read19), path(unmapped_read19) from mapped_out2823341913067520
		set val(sample), path(mapped_read20), path(unmapped_read20) from mapped_out2964509008720896
		set val(sample), path(mapped_read21), path(unmapped_read21) from mapped_out3105676104374272
		set val(sample), path(mapped_read22), path(unmapped_read22) from mapped_out3246843200027648
		set val(sample), path(mapped_read23), path(unmapped_read23) from mapped_out3388010295681024
		set val(sample), path(mapped_read24), path(unmapped_read24) from mapped_out3529177391334400
		set val(sample), path(mapped_read25), path(unmapped_read25) from mapped_out3670344486987776
		set val(sample), path(mapped_read26), path(unmapped_read26) from mapped_out3811511582641152
		set val(sample), path(mapped_read27), path(unmapped_read27) from mapped_out3952678678294528
		set val(sample), path(mapped_read28), path(unmapped_read28) from mapped_out4093845773947904
		set val(sample), path(mapped_read29), path(unmapped_read29) from mapped_out4235012869601280
		set val(sample), path(mapped_read30), path(unmapped_read30) from mapped_out4376179965254656


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout141167095653376
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq  -s ${sample}_mapped_1.fq  -s ${sample}_mapped_2.fq  -s ${sample}_mapped_3.fq  -s ${sample}_mapped_4.fq  -s ${sample}_mapped_5.fq  -s ${sample}_mapped_6.fq  -s ${sample}_mapped_7.fq  -s ${sample}_mapped_8.fq  -s ${sample}_mapped_9.fq  -s ${sample}_mapped_10.fq  -s ${sample}_mapped_11.fq  -s ${sample}_mapped_12.fq  -s ${sample}_mapped_13.fq  -s ${sample}_mapped_14.fq  -s ${sample}_mapped_15.fq  -s ${sample}_mapped_16.fq  -s ${sample}_mapped_17.fq  -s ${sample}_mapped_18.fq  -s ${sample}_mapped_19.fq  -s ${sample}_mapped_20.fq  -s ${sample}_mapped_21.fq  -s ${sample}_mapped_22.fq  -s ${sample}_mapped_23.fq  -s ${sample}_mapped_24.fq  -s ${sample}_mapped_25.fq  -s ${sample}_mapped_26.fq  -s ${sample}_mapped_27.fq  -s ${sample}_mapped_28.fq  -s ${sample}_mapped_29.fq  -s ${sample}_mapped_30.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate26 {

    publishDir "$params.output.folder26/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout141167095653376
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq141167095653376

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern26 {

    publishDir "$params.output.folder26/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq141167095653376
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads27 = params.input.fastq_path27+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads27,size:1)
params.genome27 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt27 {
    publishDir "$params.output.folder27/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out205891132094649

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		"""
        




}

process maptoreference27 {
    publishDir "$params.output.folder27/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder27/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out205891132094649

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out205891132094649

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble27 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder27/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out205891132094649


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout205891132094649
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate27 {

    publishDir "$params.output.folder27/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout205891132094649
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq205891132094649

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern27 {

    publishDir "$params.output.folder27/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq205891132094649
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads28 = params.input.fastq_path28+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads28,size:1)
params.genome28 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt28 {
    publishDir "$params.output.folder28/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out296196766695424

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		"""
        




}

process maptoreference28 {
    publishDir "$params.output.folder28/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder28/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out296196766695424

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out296196766695424

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble28 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder28/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out296196766695424


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout296196766695424
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate28 {

    publishDir "$params.output.folder28/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout296196766695424
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq296196766695424

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern28 {

    publishDir "$params.output.folder28/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq296196766695424
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads29 = params.input.fastq_path29+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads29,size:8)
params.genome29 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt29 {
    publishDir "$params.output.folder29/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out420707233300201
		set val(id),  path("${id}_trimmed_1.fq") into trim_out841414466600402
		set val(id),  path("${id}_trimmed_2.fq") into trim_out1262121699900603
		set val(id),  path("${id}_trimmed_3.fq") into trim_out1682828933200804
		set val(id),  path("${id}_trimmed_4.fq") into trim_out2103536166501005
		set val(id),  path("${id}_trimmed_5.fq") into trim_out2524243399801206
		set val(id),  path("${id}_trimmed_6.fq") into trim_out2944950633101407
		set val(id),  path("${id}_trimmed_7.fq") into trim_out3365657866401608

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
		NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq
		NanoFilt ${id}_3.fastq -l 500 -q 10 > ${id}_trimmed_3.fq
		NanoFilt ${id}_4.fastq -l 500 -q 10 > ${id}_trimmed_4.fq
		NanoFilt ${id}_5.fastq -l 500 -q 10 > ${id}_trimmed_5.fq
		NanoFilt ${id}_6.fastq -l 500 -q 10 > ${id}_trimmed_6.fq
		NanoFilt ${id}_7.fastq -l 500 -q 10 > ${id}_trimmed_7.fq
		"""
        




}

process maptoreference29 {
    publishDir "$params.output.folder29/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder29/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out420707233300201
		set val(id), path(trim_read1) from trim_out841414466600402
		set val(id), path(trim_read2) from trim_out1262121699900603
		set val(id), path(trim_read3) from trim_out1682828933200804
		set val(id), path(trim_read4) from trim_out2103536166501005
		set val(id), path(trim_read5) from trim_out2524243399801206
		set val(id), path(trim_read6) from trim_out2944950633101407
		set val(id), path(trim_read7) from trim_out3365657866401608

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out420707233300201
		set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out841414466600402
		set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out1262121699900603
		set val(id), path("${id}_mapped_3.fq"), path("${id}_unmapped_3.fq") into mapped_out1682828933200804
		set val(id), path("${id}_mapped_4.fq"), path("${id}_unmapped_4.fq") into mapped_out2103536166501005
		set val(id), path("${id}_mapped_5.fq"), path("${id}_unmapped_5.fq") into mapped_out2524243399801206
		set val(id), path("${id}_mapped_6.fq"), path("${id}_unmapped_6.fq") into mapped_out2944950633101407
		set val(id), path("${id}_mapped_7.fq"), path("${id}_unmapped_7.fq") into mapped_out3365657866401608

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_3.fq outm=${id}_mapped_3.fq outu=${id}_unmapped_3.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_4.fq outm=${id}_mapped_4.fq outu=${id}_unmapped_4.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_5.fq outm=${id}_mapped_5.fq outu=${id}_unmapped_5.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_6.fq outm=${id}_mapped_6.fq outu=${id}_unmapped_6.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_7.fq outm=${id}_mapped_7.fq outu=${id}_unmapped_7.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble29 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder29/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out420707233300201
		set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out841414466600402
		set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out1262121699900603
		set val(sample), path(mapped_read3), path(unmapped_read3) from mapped_out1682828933200804
		set val(sample), path(mapped_read4), path(unmapped_read4) from mapped_out2103536166501005
		set val(sample), path(mapped_read5), path(unmapped_read5) from mapped_out2524243399801206
		set val(sample), path(mapped_read6), path(unmapped_read6) from mapped_out2944950633101407
		set val(sample), path(mapped_read7), path(unmapped_read7) from mapped_out3365657866401608


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout420707233300201
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq  -s ${sample}_mapped_1.fq  -s ${sample}_mapped_2.fq  -s ${sample}_mapped_3.fq  -s ${sample}_mapped_4.fq  -s ${sample}_mapped_5.fq  -s ${sample}_mapped_6.fq  -s ${sample}_mapped_7.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate29 {

    publishDir "$params.output.folder29/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout420707233300201
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq420707233300201

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern29 {

    publishDir "$params.output.folder29/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq420707233300201
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads30 = params.input.fastq_path30+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads30,size:14)
params.genome30 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt30 {
    publishDir "$params.output.folder30/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out590490000000000
		set val(id),  path("${id}_trimmed_1.fq") into trim_out1180980000000000
		set val(id),  path("${id}_trimmed_2.fq") into trim_out1771470000000000
		set val(id),  path("${id}_trimmed_3.fq") into trim_out2361960000000000
		set val(id),  path("${id}_trimmed_4.fq") into trim_out2952450000000000
		set val(id),  path("${id}_trimmed_5.fq") into trim_out3542940000000000
		set val(id),  path("${id}_trimmed_6.fq") into trim_out4133430000000000
		set val(id),  path("${id}_trimmed_7.fq") into trim_out4723920000000000
		set val(id),  path("${id}_trimmed_8.fq") into trim_out5314410000000000
		set val(id),  path("${id}_trimmed_9.fq") into trim_out5904900000000000
		set val(id),  path("${id}_trimmed_10.fq") into trim_out6495390000000000
		set val(id),  path("${id}_trimmed_11.fq") into trim_out7085880000000000
		set val(id),  path("${id}_trimmed_12.fq") into trim_out7676370000000000
		set val(id),  path("${id}_trimmed_13.fq") into trim_out8266860000000000

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
		NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq
		NanoFilt ${id}_3.fastq -l 500 -q 10 > ${id}_trimmed_3.fq
		NanoFilt ${id}_4.fastq -l 500 -q 10 > ${id}_trimmed_4.fq
		NanoFilt ${id}_5.fastq -l 500 -q 10 > ${id}_trimmed_5.fq
		NanoFilt ${id}_6.fastq -l 500 -q 10 > ${id}_trimmed_6.fq
		NanoFilt ${id}_7.fastq -l 500 -q 10 > ${id}_trimmed_7.fq
		NanoFilt ${id}_8.fastq -l 500 -q 10 > ${id}_trimmed_8.fq
		NanoFilt ${id}_9.fastq -l 500 -q 10 > ${id}_trimmed_9.fq
		NanoFilt ${id}_10.fastq -l 500 -q 10 > ${id}_trimmed_10.fq
		NanoFilt ${id}_11.fastq -l 500 -q 10 > ${id}_trimmed_11.fq
		NanoFilt ${id}_12.fastq -l 500 -q 10 > ${id}_trimmed_12.fq
		NanoFilt ${id}_13.fastq -l 500 -q 10 > ${id}_trimmed_13.fq
		"""
        




}

process maptoreference30 {
    publishDir "$params.output.folder30/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder30/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out590490000000000
		set val(id), path(trim_read1) from trim_out1180980000000000
		set val(id), path(trim_read2) from trim_out1771470000000000
		set val(id), path(trim_read3) from trim_out2361960000000000
		set val(id), path(trim_read4) from trim_out2952450000000000
		set val(id), path(trim_read5) from trim_out3542940000000000
		set val(id), path(trim_read6) from trim_out4133430000000000
		set val(id), path(trim_read7) from trim_out4723920000000000
		set val(id), path(trim_read8) from trim_out5314410000000000
		set val(id), path(trim_read9) from trim_out5904900000000000
		set val(id), path(trim_read10) from trim_out6495390000000000
		set val(id), path(trim_read11) from trim_out7085880000000000
		set val(id), path(trim_read12) from trim_out7676370000000000
		set val(id), path(trim_read13) from trim_out8266860000000000

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out590490000000000
		set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out1180980000000000
		set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out1771470000000000
		set val(id), path("${id}_mapped_3.fq"), path("${id}_unmapped_3.fq") into mapped_out2361960000000000
		set val(id), path("${id}_mapped_4.fq"), path("${id}_unmapped_4.fq") into mapped_out2952450000000000
		set val(id), path("${id}_mapped_5.fq"), path("${id}_unmapped_5.fq") into mapped_out3542940000000000
		set val(id), path("${id}_mapped_6.fq"), path("${id}_unmapped_6.fq") into mapped_out4133430000000000
		set val(id), path("${id}_mapped_7.fq"), path("${id}_unmapped_7.fq") into mapped_out4723920000000000
		set val(id), path("${id}_mapped_8.fq"), path("${id}_unmapped_8.fq") into mapped_out5314410000000000
		set val(id), path("${id}_mapped_9.fq"), path("${id}_unmapped_9.fq") into mapped_out5904900000000000
		set val(id), path("${id}_mapped_10.fq"), path("${id}_unmapped_10.fq") into mapped_out6495390000000000
		set val(id), path("${id}_mapped_11.fq"), path("${id}_unmapped_11.fq") into mapped_out7085880000000000
		set val(id), path("${id}_mapped_12.fq"), path("${id}_unmapped_12.fq") into mapped_out7676370000000000
		set val(id), path("${id}_mapped_13.fq"), path("${id}_unmapped_13.fq") into mapped_out8266860000000000

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_3.fq outm=${id}_mapped_3.fq outu=${id}_unmapped_3.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_4.fq outm=${id}_mapped_4.fq outu=${id}_unmapped_4.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_5.fq outm=${id}_mapped_5.fq outu=${id}_unmapped_5.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_6.fq outm=${id}_mapped_6.fq outu=${id}_unmapped_6.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_7.fq outm=${id}_mapped_7.fq outu=${id}_unmapped_7.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_8.fq outm=${id}_mapped_8.fq outu=${id}_unmapped_8.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_9.fq outm=${id}_mapped_9.fq outu=${id}_unmapped_9.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_10.fq outm=${id}_mapped_10.fq outu=${id}_unmapped_10.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_11.fq outm=${id}_mapped_11.fq outu=${id}_unmapped_11.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_12.fq outm=${id}_mapped_12.fq outu=${id}_unmapped_12.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_13.fq outm=${id}_mapped_13.fq outu=${id}_unmapped_13.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble30 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder30/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out590490000000000
		set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out1180980000000000
		set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out1771470000000000
		set val(sample), path(mapped_read3), path(unmapped_read3) from mapped_out2361960000000000
		set val(sample), path(mapped_read4), path(unmapped_read4) from mapped_out2952450000000000
		set val(sample), path(mapped_read5), path(unmapped_read5) from mapped_out3542940000000000
		set val(sample), path(mapped_read6), path(unmapped_read6) from mapped_out4133430000000000
		set val(sample), path(mapped_read7), path(unmapped_read7) from mapped_out4723920000000000
		set val(sample), path(mapped_read8), path(unmapped_read8) from mapped_out5314410000000000
		set val(sample), path(mapped_read9), path(unmapped_read9) from mapped_out5904900000000000
		set val(sample), path(mapped_read10), path(unmapped_read10) from mapped_out6495390000000000
		set val(sample), path(mapped_read11), path(unmapped_read11) from mapped_out7085880000000000
		set val(sample), path(mapped_read12), path(unmapped_read12) from mapped_out7676370000000000
		set val(sample), path(mapped_read13), path(unmapped_read13) from mapped_out8266860000000000


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout590490000000000
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq  -s ${sample}_mapped_1.fq  -s ${sample}_mapped_2.fq  -s ${sample}_mapped_3.fq  -s ${sample}_mapped_4.fq  -s ${sample}_mapped_5.fq  -s ${sample}_mapped_6.fq  -s ${sample}_mapped_7.fq  -s ${sample}_mapped_8.fq  -s ${sample}_mapped_9.fq  -s ${sample}_mapped_10.fq  -s ${sample}_mapped_11.fq  -s ${sample}_mapped_12.fq  -s ${sample}_mapped_13.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate30 {

    publishDir "$params.output.folder30/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout590490000000000
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq590490000000000

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern30 {

    publishDir "$params.output.folder30/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq590490000000000
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads31 = params.input.fastq_path31+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads31,size:38)
params.genome31 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt31 {
    publishDir "$params.output.folder31/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out819628286980801
		set val(id),  path("${id}_trimmed_1.fq") into trim_out1639256573961602
		set val(id),  path("${id}_trimmed_2.fq") into trim_out2458884860942403
		set val(id),  path("${id}_trimmed_3.fq") into trim_out3278513147923204
		set val(id),  path("${id}_trimmed_4.fq") into trim_out4098141434904005
		set val(id),  path("${id}_trimmed_5.fq") into trim_out4917769721884806
		set val(id),  path("${id}_trimmed_6.fq") into trim_out5737398008865607
		set val(id),  path("${id}_trimmed_7.fq") into trim_out6557026295846408
		set val(id),  path("${id}_trimmed_8.fq") into trim_out7376654582827209
		set val(id),  path("${id}_trimmed_9.fq") into trim_out8196282869808010
		set val(id),  path("${id}_trimmed_10.fq") into trim_out9015911156788811
		set val(id),  path("${id}_trimmed_11.fq") into trim_out9835539443769612
		set val(id),  path("${id}_trimmed_12.fq") into trim_out10655167730750413
		set val(id),  path("${id}_trimmed_13.fq") into trim_out11474796017731214
		set val(id),  path("${id}_trimmed_14.fq") into trim_out12294424304712015
		set val(id),  path("${id}_trimmed_15.fq") into trim_out13114052591692816
		set val(id),  path("${id}_trimmed_16.fq") into trim_out13933680878673617
		set val(id),  path("${id}_trimmed_17.fq") into trim_out14753309165654418
		set val(id),  path("${id}_trimmed_18.fq") into trim_out15572937452635219
		set val(id),  path("${id}_trimmed_19.fq") into trim_out16392565739616020
		set val(id),  path("${id}_trimmed_20.fq") into trim_out17212194026596821
		set val(id),  path("${id}_trimmed_21.fq") into trim_out18031822313577622
		set val(id),  path("${id}_trimmed_22.fq") into trim_out18851450600558423
		set val(id),  path("${id}_trimmed_23.fq") into trim_out19671078887539224
		set val(id),  path("${id}_trimmed_24.fq") into trim_out20490707174520025
		set val(id),  path("${id}_trimmed_25.fq") into trim_out21310335461500826
		set val(id),  path("${id}_trimmed_26.fq") into trim_out22129963748481627
		set val(id),  path("${id}_trimmed_27.fq") into trim_out22949592035462428
		set val(id),  path("${id}_trimmed_28.fq") into trim_out23769220322443229
		set val(id),  path("${id}_trimmed_29.fq") into trim_out24588848609424030
		set val(id),  path("${id}_trimmed_30.fq") into trim_out25408476896404831
		set val(id),  path("${id}_trimmed_31.fq") into trim_out26228105183385632
		set val(id),  path("${id}_trimmed_32.fq") into trim_out27047733470366433
		set val(id),  path("${id}_trimmed_33.fq") into trim_out27867361757347234
		set val(id),  path("${id}_trimmed_34.fq") into trim_out28686990044328035
		set val(id),  path("${id}_trimmed_35.fq") into trim_out29506618331308836
		set val(id),  path("${id}_trimmed_36.fq") into trim_out30326246618289637
		set val(id),  path("${id}_trimmed_37.fq") into trim_out31145874905270438

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
		NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq
		NanoFilt ${id}_3.fastq -l 500 -q 10 > ${id}_trimmed_3.fq
		NanoFilt ${id}_4.fastq -l 500 -q 10 > ${id}_trimmed_4.fq
		NanoFilt ${id}_5.fastq -l 500 -q 10 > ${id}_trimmed_5.fq
		NanoFilt ${id}_6.fastq -l 500 -q 10 > ${id}_trimmed_6.fq
		NanoFilt ${id}_7.fastq -l 500 -q 10 > ${id}_trimmed_7.fq
		NanoFilt ${id}_8.fastq -l 500 -q 10 > ${id}_trimmed_8.fq
		NanoFilt ${id}_9.fastq -l 500 -q 10 > ${id}_trimmed_9.fq
		NanoFilt ${id}_10.fastq -l 500 -q 10 > ${id}_trimmed_10.fq
		NanoFilt ${id}_11.fastq -l 500 -q 10 > ${id}_trimmed_11.fq
		NanoFilt ${id}_12.fastq -l 500 -q 10 > ${id}_trimmed_12.fq
		NanoFilt ${id}_13.fastq -l 500 -q 10 > ${id}_trimmed_13.fq
		NanoFilt ${id}_14.fastq -l 500 -q 10 > ${id}_trimmed_14.fq
		NanoFilt ${id}_15.fastq -l 500 -q 10 > ${id}_trimmed_15.fq
		NanoFilt ${id}_16.fastq -l 500 -q 10 > ${id}_trimmed_16.fq
		NanoFilt ${id}_17.fastq -l 500 -q 10 > ${id}_trimmed_17.fq
		NanoFilt ${id}_18.fastq -l 500 -q 10 > ${id}_trimmed_18.fq
		NanoFilt ${id}_19.fastq -l 500 -q 10 > ${id}_trimmed_19.fq
		NanoFilt ${id}_20.fastq -l 500 -q 10 > ${id}_trimmed_20.fq
		NanoFilt ${id}_21.fastq -l 500 -q 10 > ${id}_trimmed_21.fq
		NanoFilt ${id}_22.fastq -l 500 -q 10 > ${id}_trimmed_22.fq
		NanoFilt ${id}_23.fastq -l 500 -q 10 > ${id}_trimmed_23.fq
		NanoFilt ${id}_24.fastq -l 500 -q 10 > ${id}_trimmed_24.fq
		NanoFilt ${id}_25.fastq -l 500 -q 10 > ${id}_trimmed_25.fq
		NanoFilt ${id}_26.fastq -l 500 -q 10 > ${id}_trimmed_26.fq
		NanoFilt ${id}_27.fastq -l 500 -q 10 > ${id}_trimmed_27.fq
		NanoFilt ${id}_28.fastq -l 500 -q 10 > ${id}_trimmed_28.fq
		NanoFilt ${id}_29.fastq -l 500 -q 10 > ${id}_trimmed_29.fq
		NanoFilt ${id}_30.fastq -l 500 -q 10 > ${id}_trimmed_30.fq
		NanoFilt ${id}_31.fastq -l 500 -q 10 > ${id}_trimmed_31.fq
		NanoFilt ${id}_32.fastq -l 500 -q 10 > ${id}_trimmed_32.fq
		NanoFilt ${id}_33.fastq -l 500 -q 10 > ${id}_trimmed_33.fq
		NanoFilt ${id}_34.fastq -l 500 -q 10 > ${id}_trimmed_34.fq
		NanoFilt ${id}_35.fastq -l 500 -q 10 > ${id}_trimmed_35.fq
		NanoFilt ${id}_36.fastq -l 500 -q 10 > ${id}_trimmed_36.fq
		NanoFilt ${id}_37.fastq -l 500 -q 10 > ${id}_trimmed_37.fq
		"""
        




}

process maptoreference31 {
    publishDir "$params.output.folder31/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder31/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out819628286980801
		set val(id), path(trim_read1) from trim_out1639256573961602
		set val(id), path(trim_read2) from trim_out2458884860942403
		set val(id), path(trim_read3) from trim_out3278513147923204
		set val(id), path(trim_read4) from trim_out4098141434904005
		set val(id), path(trim_read5) from trim_out4917769721884806
		set val(id), path(trim_read6) from trim_out5737398008865607
		set val(id), path(trim_read7) from trim_out6557026295846408
		set val(id), path(trim_read8) from trim_out7376654582827209
		set val(id), path(trim_read9) from trim_out8196282869808010
		set val(id), path(trim_read10) from trim_out9015911156788811
		set val(id), path(trim_read11) from trim_out9835539443769612
		set val(id), path(trim_read12) from trim_out10655167730750413
		set val(id), path(trim_read13) from trim_out11474796017731214
		set val(id), path(trim_read14) from trim_out12294424304712015
		set val(id), path(trim_read15) from trim_out13114052591692816
		set val(id), path(trim_read16) from trim_out13933680878673617
		set val(id), path(trim_read17) from trim_out14753309165654418
		set val(id), path(trim_read18) from trim_out15572937452635219
		set val(id), path(trim_read19) from trim_out16392565739616020
		set val(id), path(trim_read20) from trim_out17212194026596821
		set val(id), path(trim_read21) from trim_out18031822313577622
		set val(id), path(trim_read22) from trim_out18851450600558423
		set val(id), path(trim_read23) from trim_out19671078887539224
		set val(id), path(trim_read24) from trim_out20490707174520025
		set val(id), path(trim_read25) from trim_out21310335461500826
		set val(id), path(trim_read26) from trim_out22129963748481627
		set val(id), path(trim_read27) from trim_out22949592035462428
		set val(id), path(trim_read28) from trim_out23769220322443229
		set val(id), path(trim_read29) from trim_out24588848609424030
		set val(id), path(trim_read30) from trim_out25408476896404831
		set val(id), path(trim_read31) from trim_out26228105183385632
		set val(id), path(trim_read32) from trim_out27047733470366433
		set val(id), path(trim_read33) from trim_out27867361757347234
		set val(id), path(trim_read34) from trim_out28686990044328035
		set val(id), path(trim_read35) from trim_out29506618331308836
		set val(id), path(trim_read36) from trim_out30326246618289637
		set val(id), path(trim_read37) from trim_out31145874905270438

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out819628286980801
		set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out1639256573961602
		set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out2458884860942403
		set val(id), path("${id}_mapped_3.fq"), path("${id}_unmapped_3.fq") into mapped_out3278513147923204
		set val(id), path("${id}_mapped_4.fq"), path("${id}_unmapped_4.fq") into mapped_out4098141434904005
		set val(id), path("${id}_mapped_5.fq"), path("${id}_unmapped_5.fq") into mapped_out4917769721884806
		set val(id), path("${id}_mapped_6.fq"), path("${id}_unmapped_6.fq") into mapped_out5737398008865607
		set val(id), path("${id}_mapped_7.fq"), path("${id}_unmapped_7.fq") into mapped_out6557026295846408
		set val(id), path("${id}_mapped_8.fq"), path("${id}_unmapped_8.fq") into mapped_out7376654582827209
		set val(id), path("${id}_mapped_9.fq"), path("${id}_unmapped_9.fq") into mapped_out8196282869808010
		set val(id), path("${id}_mapped_10.fq"), path("${id}_unmapped_10.fq") into mapped_out9015911156788811
		set val(id), path("${id}_mapped_11.fq"), path("${id}_unmapped_11.fq") into mapped_out9835539443769612
		set val(id), path("${id}_mapped_12.fq"), path("${id}_unmapped_12.fq") into mapped_out10655167730750413
		set val(id), path("${id}_mapped_13.fq"), path("${id}_unmapped_13.fq") into mapped_out11474796017731214
		set val(id), path("${id}_mapped_14.fq"), path("${id}_unmapped_14.fq") into mapped_out12294424304712015
		set val(id), path("${id}_mapped_15.fq"), path("${id}_unmapped_15.fq") into mapped_out13114052591692816
		set val(id), path("${id}_mapped_16.fq"), path("${id}_unmapped_16.fq") into mapped_out13933680878673617
		set val(id), path("${id}_mapped_17.fq"), path("${id}_unmapped_17.fq") into mapped_out14753309165654418
		set val(id), path("${id}_mapped_18.fq"), path("${id}_unmapped_18.fq") into mapped_out15572937452635219
		set val(id), path("${id}_mapped_19.fq"), path("${id}_unmapped_19.fq") into mapped_out16392565739616020
		set val(id), path("${id}_mapped_20.fq"), path("${id}_unmapped_20.fq") into mapped_out17212194026596821
		set val(id), path("${id}_mapped_21.fq"), path("${id}_unmapped_21.fq") into mapped_out18031822313577622
		set val(id), path("${id}_mapped_22.fq"), path("${id}_unmapped_22.fq") into mapped_out18851450600558423
		set val(id), path("${id}_mapped_23.fq"), path("${id}_unmapped_23.fq") into mapped_out19671078887539224
		set val(id), path("${id}_mapped_24.fq"), path("${id}_unmapped_24.fq") into mapped_out20490707174520025
		set val(id), path("${id}_mapped_25.fq"), path("${id}_unmapped_25.fq") into mapped_out21310335461500826
		set val(id), path("${id}_mapped_26.fq"), path("${id}_unmapped_26.fq") into mapped_out22129963748481627
		set val(id), path("${id}_mapped_27.fq"), path("${id}_unmapped_27.fq") into mapped_out22949592035462428
		set val(id), path("${id}_mapped_28.fq"), path("${id}_unmapped_28.fq") into mapped_out23769220322443229
		set val(id), path("${id}_mapped_29.fq"), path("${id}_unmapped_29.fq") into mapped_out24588848609424030
		set val(id), path("${id}_mapped_30.fq"), path("${id}_unmapped_30.fq") into mapped_out25408476896404831
		set val(id), path("${id}_mapped_31.fq"), path("${id}_unmapped_31.fq") into mapped_out26228105183385632
		set val(id), path("${id}_mapped_32.fq"), path("${id}_unmapped_32.fq") into mapped_out27047733470366433
		set val(id), path("${id}_mapped_33.fq"), path("${id}_unmapped_33.fq") into mapped_out27867361757347234
		set val(id), path("${id}_mapped_34.fq"), path("${id}_unmapped_34.fq") into mapped_out28686990044328035
		set val(id), path("${id}_mapped_35.fq"), path("${id}_unmapped_35.fq") into mapped_out29506618331308836
		set val(id), path("${id}_mapped_36.fq"), path("${id}_unmapped_36.fq") into mapped_out30326246618289637
		set val(id), path("${id}_mapped_37.fq"), path("${id}_unmapped_37.fq") into mapped_out31145874905270438

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_3.fq outm=${id}_mapped_3.fq outu=${id}_unmapped_3.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_4.fq outm=${id}_mapped_4.fq outu=${id}_unmapped_4.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_5.fq outm=${id}_mapped_5.fq outu=${id}_unmapped_5.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_6.fq outm=${id}_mapped_6.fq outu=${id}_unmapped_6.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_7.fq outm=${id}_mapped_7.fq outu=${id}_unmapped_7.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_8.fq outm=${id}_mapped_8.fq outu=${id}_unmapped_8.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_9.fq outm=${id}_mapped_9.fq outu=${id}_unmapped_9.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_10.fq outm=${id}_mapped_10.fq outu=${id}_unmapped_10.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_11.fq outm=${id}_mapped_11.fq outu=${id}_unmapped_11.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_12.fq outm=${id}_mapped_12.fq outu=${id}_unmapped_12.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_13.fq outm=${id}_mapped_13.fq outu=${id}_unmapped_13.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_14.fq outm=${id}_mapped_14.fq outu=${id}_unmapped_14.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_15.fq outm=${id}_mapped_15.fq outu=${id}_unmapped_15.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_16.fq outm=${id}_mapped_16.fq outu=${id}_unmapped_16.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_17.fq outm=${id}_mapped_17.fq outu=${id}_unmapped_17.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_18.fq outm=${id}_mapped_18.fq outu=${id}_unmapped_18.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_19.fq outm=${id}_mapped_19.fq outu=${id}_unmapped_19.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_20.fq outm=${id}_mapped_20.fq outu=${id}_unmapped_20.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_21.fq outm=${id}_mapped_21.fq outu=${id}_unmapped_21.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_22.fq outm=${id}_mapped_22.fq outu=${id}_unmapped_22.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_23.fq outm=${id}_mapped_23.fq outu=${id}_unmapped_23.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_24.fq outm=${id}_mapped_24.fq outu=${id}_unmapped_24.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_25.fq outm=${id}_mapped_25.fq outu=${id}_unmapped_25.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_26.fq outm=${id}_mapped_26.fq outu=${id}_unmapped_26.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_27.fq outm=${id}_mapped_27.fq outu=${id}_unmapped_27.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_28.fq outm=${id}_mapped_28.fq outu=${id}_unmapped_28.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_29.fq outm=${id}_mapped_29.fq outu=${id}_unmapped_29.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_30.fq outm=${id}_mapped_30.fq outu=${id}_unmapped_30.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_31.fq outm=${id}_mapped_31.fq outu=${id}_unmapped_31.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_32.fq outm=${id}_mapped_32.fq outu=${id}_unmapped_32.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_33.fq outm=${id}_mapped_33.fq outu=${id}_unmapped_33.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_34.fq outm=${id}_mapped_34.fq outu=${id}_unmapped_34.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_35.fq outm=${id}_mapped_35.fq outu=${id}_unmapped_35.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_36.fq outm=${id}_mapped_36.fq outu=${id}_unmapped_36.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_37.fq outm=${id}_mapped_37.fq outu=${id}_unmapped_37.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble31 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder31/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out819628286980801
		set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out1639256573961602
		set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out2458884860942403
		set val(sample), path(mapped_read3), path(unmapped_read3) from mapped_out3278513147923204
		set val(sample), path(mapped_read4), path(unmapped_read4) from mapped_out4098141434904005
		set val(sample), path(mapped_read5), path(unmapped_read5) from mapped_out4917769721884806
		set val(sample), path(mapped_read6), path(unmapped_read6) from mapped_out5737398008865607
		set val(sample), path(mapped_read7), path(unmapped_read7) from mapped_out6557026295846408
		set val(sample), path(mapped_read8), path(unmapped_read8) from mapped_out7376654582827209
		set val(sample), path(mapped_read9), path(unmapped_read9) from mapped_out8196282869808010
		set val(sample), path(mapped_read10), path(unmapped_read10) from mapped_out9015911156788811
		set val(sample), path(mapped_read11), path(unmapped_read11) from mapped_out9835539443769612
		set val(sample), path(mapped_read12), path(unmapped_read12) from mapped_out10655167730750413
		set val(sample), path(mapped_read13), path(unmapped_read13) from mapped_out11474796017731214
		set val(sample), path(mapped_read14), path(unmapped_read14) from mapped_out12294424304712015
		set val(sample), path(mapped_read15), path(unmapped_read15) from mapped_out13114052591692816
		set val(sample), path(mapped_read16), path(unmapped_read16) from mapped_out13933680878673617
		set val(sample), path(mapped_read17), path(unmapped_read17) from mapped_out14753309165654418
		set val(sample), path(mapped_read18), path(unmapped_read18) from mapped_out15572937452635219
		set val(sample), path(mapped_read19), path(unmapped_read19) from mapped_out16392565739616020
		set val(sample), path(mapped_read20), path(unmapped_read20) from mapped_out17212194026596821
		set val(sample), path(mapped_read21), path(unmapped_read21) from mapped_out18031822313577622
		set val(sample), path(mapped_read22), path(unmapped_read22) from mapped_out18851450600558423
		set val(sample), path(mapped_read23), path(unmapped_read23) from mapped_out19671078887539224
		set val(sample), path(mapped_read24), path(unmapped_read24) from mapped_out20490707174520025
		set val(sample), path(mapped_read25), path(unmapped_read25) from mapped_out21310335461500826
		set val(sample), path(mapped_read26), path(unmapped_read26) from mapped_out22129963748481627
		set val(sample), path(mapped_read27), path(unmapped_read27) from mapped_out22949592035462428
		set val(sample), path(mapped_read28), path(unmapped_read28) from mapped_out23769220322443229
		set val(sample), path(mapped_read29), path(unmapped_read29) from mapped_out24588848609424030
		set val(sample), path(mapped_read30), path(unmapped_read30) from mapped_out25408476896404831
		set val(sample), path(mapped_read31), path(unmapped_read31) from mapped_out26228105183385632
		set val(sample), path(mapped_read32), path(unmapped_read32) from mapped_out27047733470366433
		set val(sample), path(mapped_read33), path(unmapped_read33) from mapped_out27867361757347234
		set val(sample), path(mapped_read34), path(unmapped_read34) from mapped_out28686990044328035
		set val(sample), path(mapped_read35), path(unmapped_read35) from mapped_out29506618331308836
		set val(sample), path(mapped_read36), path(unmapped_read36) from mapped_out30326246618289637
		set val(sample), path(mapped_read37), path(unmapped_read37) from mapped_out31145874905270438


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout819628286980801
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq  -s ${sample}_mapped_1.fq  -s ${sample}_mapped_2.fq  -s ${sample}_mapped_3.fq  -s ${sample}_mapped_4.fq  -s ${sample}_mapped_5.fq  -s ${sample}_mapped_6.fq  -s ${sample}_mapped_7.fq  -s ${sample}_mapped_8.fq  -s ${sample}_mapped_9.fq  -s ${sample}_mapped_10.fq  -s ${sample}_mapped_11.fq  -s ${sample}_mapped_12.fq  -s ${sample}_mapped_13.fq  -s ${sample}_mapped_14.fq  -s ${sample}_mapped_15.fq  -s ${sample}_mapped_16.fq  -s ${sample}_mapped_17.fq  -s ${sample}_mapped_18.fq  -s ${sample}_mapped_19.fq  -s ${sample}_mapped_20.fq  -s ${sample}_mapped_21.fq  -s ${sample}_mapped_22.fq  -s ${sample}_mapped_23.fq  -s ${sample}_mapped_24.fq  -s ${sample}_mapped_25.fq  -s ${sample}_mapped_26.fq  -s ${sample}_mapped_27.fq  -s ${sample}_mapped_28.fq  -s ${sample}_mapped_29.fq  -s ${sample}_mapped_30.fq  -s ${sample}_mapped_31.fq  -s ${sample}_mapped_32.fq  -s ${sample}_mapped_33.fq  -s ${sample}_mapped_34.fq  -s ${sample}_mapped_35.fq  -s ${sample}_mapped_36.fq  -s ${sample}_mapped_37.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate31 {

    publishDir "$params.output.folder31/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout819628286980801
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq819628286980801

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern31 {

    publishDir "$params.output.folder31/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq819628286980801
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}


params.reads32 = params.input.fastq_path32+'/*_{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199}.fastq'
fastq_path = Channel.fromFilePairs(params.reads32,size:183)
params.genome32 = "$baseDir/ref/XM_002808697.2.fasta"
pyscripts="$baseDir/pyscripts"
ref_dir="$baseDir/ref"

process trim_nanofilt32 {
    publishDir "$params.output.folder32/trimmedfastq/${id}", pattern: "*.fq", mode : "copy"
    
    input:
        set val(id), path(fastq_group) from fastq_path
        //set val(id), path(fastq_group2) from fastq_path2
        //set val(id), path(fastq_group3) from fastq_path3
       
    output:
		set val(id),  path("${id}_trimmed_0.fq") into trim_out1125899906842624
		set val(id),  path("${id}_trimmed_1.fq") into trim_out2251799813685248
		set val(id),  path("${id}_trimmed_2.fq") into trim_out3377699720527872
		set val(id),  path("${id}_trimmed_3.fq") into trim_out4503599627370496
		set val(id),  path("${id}_trimmed_4.fq") into trim_out5629499534213120
		set val(id),  path("${id}_trimmed_5.fq") into trim_out6755399441055744
		set val(id),  path("${id}_trimmed_6.fq") into trim_out7881299347898368
		set val(id),  path("${id}_trimmed_7.fq") into trim_out9007199254740992
		set val(id),  path("${id}_trimmed_8.fq") into trim_out10133099161583616
		set val(id),  path("${id}_trimmed_9.fq") into trim_out11258999068426240
		set val(id),  path("${id}_trimmed_10.fq") into trim_out12384898975268864
		set val(id),  path("${id}_trimmed_11.fq") into trim_out13510798882111488
		set val(id),  path("${id}_trimmed_12.fq") into trim_out14636698788954112
		set val(id),  path("${id}_trimmed_13.fq") into trim_out15762598695796736
		set val(id),  path("${id}_trimmed_14.fq") into trim_out16888498602639360
		set val(id),  path("${id}_trimmed_15.fq") into trim_out18014398509481984
		set val(id),  path("${id}_trimmed_16.fq") into trim_out19140298416324608
		set val(id),  path("${id}_trimmed_17.fq") into trim_out20266198323167232
		set val(id),  path("${id}_trimmed_18.fq") into trim_out21392098230009856
		set val(id),  path("${id}_trimmed_19.fq") into trim_out22517998136852480
		set val(id),  path("${id}_trimmed_20.fq") into trim_out23643898043695104
		set val(id),  path("${id}_trimmed_21.fq") into trim_out24769797950537728
		set val(id),  path("${id}_trimmed_22.fq") into trim_out25895697857380352
		set val(id),  path("${id}_trimmed_23.fq") into trim_out27021597764222976
		set val(id),  path("${id}_trimmed_24.fq") into trim_out28147497671065600
		set val(id),  path("${id}_trimmed_25.fq") into trim_out29273397577908224
		set val(id),  path("${id}_trimmed_26.fq") into trim_out30399297484750848
		set val(id),  path("${id}_trimmed_27.fq") into trim_out31525197391593472
		set val(id),  path("${id}_trimmed_28.fq") into trim_out32651097298436096
		set val(id),  path("${id}_trimmed_29.fq") into trim_out33776997205278720
		set val(id),  path("${id}_trimmed_30.fq") into trim_out34902897112121344
		set val(id),  path("${id}_trimmed_31.fq") into trim_out36028797018963968
		set val(id),  path("${id}_trimmed_32.fq") into trim_out37154696925806592
		set val(id),  path("${id}_trimmed_33.fq") into trim_out38280596832649216
		set val(id),  path("${id}_trimmed_34.fq") into trim_out39406496739491840
		set val(id),  path("${id}_trimmed_35.fq") into trim_out40532396646334464
		set val(id),  path("${id}_trimmed_36.fq") into trim_out41658296553177088
		set val(id),  path("${id}_trimmed_37.fq") into trim_out42784196460019712
		set val(id),  path("${id}_trimmed_38.fq") into trim_out43910096366862336
		set val(id),  path("${id}_trimmed_39.fq") into trim_out45035996273704960
		set val(id),  path("${id}_trimmed_40.fq") into trim_out46161896180547584
		set val(id),  path("${id}_trimmed_41.fq") into trim_out47287796087390208
		set val(id),  path("${id}_trimmed_42.fq") into trim_out48413695994232832
		set val(id),  path("${id}_trimmed_43.fq") into trim_out49539595901075456
		set val(id),  path("${id}_trimmed_44.fq") into trim_out50665495807918080
		set val(id),  path("${id}_trimmed_45.fq") into trim_out51791395714760704
		set val(id),  path("${id}_trimmed_46.fq") into trim_out52917295621603328
		set val(id),  path("${id}_trimmed_47.fq") into trim_out54043195528445952
		set val(id),  path("${id}_trimmed_48.fq") into trim_out55169095435288576
		set val(id),  path("${id}_trimmed_49.fq") into trim_out56294995342131200
		set val(id),  path("${id}_trimmed_50.fq") into trim_out57420895248973824
		set val(id),  path("${id}_trimmed_51.fq") into trim_out58546795155816448
		set val(id),  path("${id}_trimmed_52.fq") into trim_out59672695062659072
		set val(id),  path("${id}_trimmed_53.fq") into trim_out60798594969501696
		set val(id),  path("${id}_trimmed_54.fq") into trim_out61924494876344320
		set val(id),  path("${id}_trimmed_55.fq") into trim_out63050394783186944
		set val(id),  path("${id}_trimmed_56.fq") into trim_out64176294690029568
		set val(id),  path("${id}_trimmed_57.fq") into trim_out65302194596872192
		set val(id),  path("${id}_trimmed_58.fq") into trim_out66428094503714816
		set val(id),  path("${id}_trimmed_59.fq") into trim_out67553994410557440
		set val(id),  path("${id}_trimmed_60.fq") into trim_out68679894317400064
		set val(id),  path("${id}_trimmed_61.fq") into trim_out69805794224242688
		set val(id),  path("${id}_trimmed_62.fq") into trim_out70931694131085312
		set val(id),  path("${id}_trimmed_63.fq") into trim_out72057594037927936
		set val(id),  path("${id}_trimmed_64.fq") into trim_out73183493944770560
		set val(id),  path("${id}_trimmed_65.fq") into trim_out74309393851613184
		set val(id),  path("${id}_trimmed_66.fq") into trim_out75435293758455808
		set val(id),  path("${id}_trimmed_67.fq") into trim_out76561193665298432
		set val(id),  path("${id}_trimmed_68.fq") into trim_out77687093572141056
		set val(id),  path("${id}_trimmed_69.fq") into trim_out78812993478983680
		set val(id),  path("${id}_trimmed_70.fq") into trim_out79938893385826304
		set val(id),  path("${id}_trimmed_71.fq") into trim_out81064793292668928
		set val(id),  path("${id}_trimmed_72.fq") into trim_out82190693199511552
		set val(id),  path("${id}_trimmed_73.fq") into trim_out83316593106354176
		set val(id),  path("${id}_trimmed_74.fq") into trim_out84442493013196800
		set val(id),  path("${id}_trimmed_75.fq") into trim_out85568392920039424
		set val(id),  path("${id}_trimmed_76.fq") into trim_out86694292826882048
		set val(id),  path("${id}_trimmed_77.fq") into trim_out87820192733724672
		set val(id),  path("${id}_trimmed_78.fq") into trim_out88946092640567296
		set val(id),  path("${id}_trimmed_79.fq") into trim_out90071992547409920
		set val(id),  path("${id}_trimmed_80.fq") into trim_out91197892454252544
		set val(id),  path("${id}_trimmed_81.fq") into trim_out92323792361095168
		set val(id),  path("${id}_trimmed_82.fq") into trim_out93449692267937792
		set val(id),  path("${id}_trimmed_83.fq") into trim_out94575592174780416
		set val(id),  path("${id}_trimmed_84.fq") into trim_out95701492081623040
		set val(id),  path("${id}_trimmed_85.fq") into trim_out96827391988465664
		set val(id),  path("${id}_trimmed_86.fq") into trim_out97953291895308288
		set val(id),  path("${id}_trimmed_87.fq") into trim_out99079191802150912
		set val(id),  path("${id}_trimmed_88.fq") into trim_out100205091708993536
		set val(id),  path("${id}_trimmed_89.fq") into trim_out101330991615836160
		set val(id),  path("${id}_trimmed_90.fq") into trim_out102456891522678784
		set val(id),  path("${id}_trimmed_91.fq") into trim_out103582791429521408
		set val(id),  path("${id}_trimmed_92.fq") into trim_out104708691336364032
		set val(id),  path("${id}_trimmed_93.fq") into trim_out105834591243206656
		set val(id),  path("${id}_trimmed_94.fq") into trim_out106960491150049280
		set val(id),  path("${id}_trimmed_95.fq") into trim_out108086391056891904
		set val(id),  path("${id}_trimmed_96.fq") into trim_out109212290963734528
		set val(id),  path("${id}_trimmed_97.fq") into trim_out110338190870577152
		set val(id),  path("${id}_trimmed_98.fq") into trim_out111464090777419776
		set val(id),  path("${id}_trimmed_99.fq") into trim_out112589990684262400
		set val(id),  path("${id}_trimmed_100.fq") into trim_out113715890591105024
		set val(id),  path("${id}_trimmed_101.fq") into trim_out114841790497947648
		set val(id),  path("${id}_trimmed_102.fq") into trim_out115967690404790272
		set val(id),  path("${id}_trimmed_103.fq") into trim_out117093590311632896
		set val(id),  path("${id}_trimmed_104.fq") into trim_out118219490218475520
		set val(id),  path("${id}_trimmed_105.fq") into trim_out119345390125318144
		set val(id),  path("${id}_trimmed_106.fq") into trim_out120471290032160768
		set val(id),  path("${id}_trimmed_107.fq") into trim_out121597189939003392
		set val(id),  path("${id}_trimmed_108.fq") into trim_out122723089845846016
		set val(id),  path("${id}_trimmed_109.fq") into trim_out123848989752688640
		set val(id),  path("${id}_trimmed_110.fq") into trim_out124974889659531264
		set val(id),  path("${id}_trimmed_111.fq") into trim_out126100789566373888
		set val(id),  path("${id}_trimmed_112.fq") into trim_out127226689473216512
		set val(id),  path("${id}_trimmed_113.fq") into trim_out128352589380059136
		set val(id),  path("${id}_trimmed_114.fq") into trim_out129478489286901760
		set val(id),  path("${id}_trimmed_115.fq") into trim_out130604389193744384
		set val(id),  path("${id}_trimmed_116.fq") into trim_out131730289100587008
		set val(id),  path("${id}_trimmed_117.fq") into trim_out132856189007429632
		set val(id),  path("${id}_trimmed_118.fq") into trim_out133982088914272256
		set val(id),  path("${id}_trimmed_119.fq") into trim_out135107988821114880
		set val(id),  path("${id}_trimmed_120.fq") into trim_out136233888727957504
		set val(id),  path("${id}_trimmed_121.fq") into trim_out137359788634800128
		set val(id),  path("${id}_trimmed_122.fq") into trim_out138485688541642752
		set val(id),  path("${id}_trimmed_123.fq") into trim_out139611588448485376
		set val(id),  path("${id}_trimmed_124.fq") into trim_out140737488355328000
		set val(id),  path("${id}_trimmed_125.fq") into trim_out141863388262170624
		set val(id),  path("${id}_trimmed_126.fq") into trim_out142989288169013248
		set val(id),  path("${id}_trimmed_127.fq") into trim_out144115188075855872
		set val(id),  path("${id}_trimmed_128.fq") into trim_out145241087982698496
		set val(id),  path("${id}_trimmed_129.fq") into trim_out146366987889541120
		set val(id),  path("${id}_trimmed_130.fq") into trim_out147492887796383744
		set val(id),  path("${id}_trimmed_131.fq") into trim_out148618787703226368
		set val(id),  path("${id}_trimmed_132.fq") into trim_out149744687610068992
		set val(id),  path("${id}_trimmed_133.fq") into trim_out150870587516911616
		set val(id),  path("${id}_trimmed_134.fq") into trim_out151996487423754240
		set val(id),  path("${id}_trimmed_135.fq") into trim_out153122387330596864
		set val(id),  path("${id}_trimmed_136.fq") into trim_out154248287237439488
		set val(id),  path("${id}_trimmed_137.fq") into trim_out155374187144282112
		set val(id),  path("${id}_trimmed_138.fq") into trim_out156500087051124736
		set val(id),  path("${id}_trimmed_139.fq") into trim_out157625986957967360
		set val(id),  path("${id}_trimmed_140.fq") into trim_out158751886864809984
		set val(id),  path("${id}_trimmed_141.fq") into trim_out159877786771652608
		set val(id),  path("${id}_trimmed_142.fq") into trim_out161003686678495232
		set val(id),  path("${id}_trimmed_143.fq") into trim_out162129586585337856
		set val(id),  path("${id}_trimmed_144.fq") into trim_out163255486492180480
		set val(id),  path("${id}_trimmed_145.fq") into trim_out164381386399023104
		set val(id),  path("${id}_trimmed_146.fq") into trim_out165507286305865728
		set val(id),  path("${id}_trimmed_147.fq") into trim_out166633186212708352
		set val(id),  path("${id}_trimmed_148.fq") into trim_out167759086119550976
		set val(id),  path("${id}_trimmed_149.fq") into trim_out168884986026393600
		set val(id),  path("${id}_trimmed_150.fq") into trim_out170010885933236224
		set val(id),  path("${id}_trimmed_151.fq") into trim_out171136785840078848
		set val(id),  path("${id}_trimmed_152.fq") into trim_out172262685746921472
		set val(id),  path("${id}_trimmed_153.fq") into trim_out173388585653764096
		set val(id),  path("${id}_trimmed_154.fq") into trim_out174514485560606720
		set val(id),  path("${id}_trimmed_155.fq") into trim_out175640385467449344
		set val(id),  path("${id}_trimmed_156.fq") into trim_out176766285374291968
		set val(id),  path("${id}_trimmed_157.fq") into trim_out177892185281134592
		set val(id),  path("${id}_trimmed_158.fq") into trim_out179018085187977216
		set val(id),  path("${id}_trimmed_159.fq") into trim_out180143985094819840
		set val(id),  path("${id}_trimmed_160.fq") into trim_out181269885001662464
		set val(id),  path("${id}_trimmed_161.fq") into trim_out182395784908505088
		set val(id),  path("${id}_trimmed_162.fq") into trim_out183521684815347712
		set val(id),  path("${id}_trimmed_163.fq") into trim_out184647584722190336
		set val(id),  path("${id}_trimmed_164.fq") into trim_out185773484629032960
		set val(id),  path("${id}_trimmed_165.fq") into trim_out186899384535875584
		set val(id),  path("${id}_trimmed_166.fq") into trim_out188025284442718208
		set val(id),  path("${id}_trimmed_167.fq") into trim_out189151184349560832
		set val(id),  path("${id}_trimmed_168.fq") into trim_out190277084256403456
		set val(id),  path("${id}_trimmed_169.fq") into trim_out191402984163246080
		set val(id),  path("${id}_trimmed_170.fq") into trim_out192528884070088704
		set val(id),  path("${id}_trimmed_171.fq") into trim_out193654783976931328
		set val(id),  path("${id}_trimmed_172.fq") into trim_out194780683883773952
		set val(id),  path("${id}_trimmed_173.fq") into trim_out195906583790616576
		set val(id),  path("${id}_trimmed_174.fq") into trim_out197032483697459200
		set val(id),  path("${id}_trimmed_175.fq") into trim_out198158383604301824
		set val(id),  path("${id}_trimmed_176.fq") into trim_out199284283511144448
		set val(id),  path("${id}_trimmed_177.fq") into trim_out200410183417987072
		set val(id),  path("${id}_trimmed_178.fq") into trim_out201536083324829696
		set val(id),  path("${id}_trimmed_179.fq") into trim_out202661983231672320
		set val(id),  path("${id}_trimmed_180.fq") into trim_out203787883138514944
		set val(id),  path("${id}_trimmed_181.fq") into trim_out204913783045357568
		set val(id),  path("${id}_trimmed_182.fq") into trim_out206039682952200192

    script:
		"""
		NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq
		NanoFilt ${id}_1.fastq -l 500 -q 10 > ${id}_trimmed_1.fq
		NanoFilt ${id}_2.fastq -l 500 -q 10 > ${id}_trimmed_2.fq
		NanoFilt ${id}_3.fastq -l 500 -q 10 > ${id}_trimmed_3.fq
		NanoFilt ${id}_4.fastq -l 500 -q 10 > ${id}_trimmed_4.fq
		NanoFilt ${id}_5.fastq -l 500 -q 10 > ${id}_trimmed_5.fq
		NanoFilt ${id}_6.fastq -l 500 -q 10 > ${id}_trimmed_6.fq
		NanoFilt ${id}_7.fastq -l 500 -q 10 > ${id}_trimmed_7.fq
		NanoFilt ${id}_8.fastq -l 500 -q 10 > ${id}_trimmed_8.fq
		NanoFilt ${id}_9.fastq -l 500 -q 10 > ${id}_trimmed_9.fq
		NanoFilt ${id}_10.fastq -l 500 -q 10 > ${id}_trimmed_10.fq
		NanoFilt ${id}_11.fastq -l 500 -q 10 > ${id}_trimmed_11.fq
		NanoFilt ${id}_12.fastq -l 500 -q 10 > ${id}_trimmed_12.fq
		NanoFilt ${id}_13.fastq -l 500 -q 10 > ${id}_trimmed_13.fq
		NanoFilt ${id}_14.fastq -l 500 -q 10 > ${id}_trimmed_14.fq
		NanoFilt ${id}_15.fastq -l 500 -q 10 > ${id}_trimmed_15.fq
		NanoFilt ${id}_16.fastq -l 500 -q 10 > ${id}_trimmed_16.fq
		NanoFilt ${id}_17.fastq -l 500 -q 10 > ${id}_trimmed_17.fq
		NanoFilt ${id}_18.fastq -l 500 -q 10 > ${id}_trimmed_18.fq
		NanoFilt ${id}_19.fastq -l 500 -q 10 > ${id}_trimmed_19.fq
		NanoFilt ${id}_20.fastq -l 500 -q 10 > ${id}_trimmed_20.fq
		NanoFilt ${id}_21.fastq -l 500 -q 10 > ${id}_trimmed_21.fq
		NanoFilt ${id}_22.fastq -l 500 -q 10 > ${id}_trimmed_22.fq
		NanoFilt ${id}_23.fastq -l 500 -q 10 > ${id}_trimmed_23.fq
		NanoFilt ${id}_24.fastq -l 500 -q 10 > ${id}_trimmed_24.fq
		NanoFilt ${id}_25.fastq -l 500 -q 10 > ${id}_trimmed_25.fq
		NanoFilt ${id}_26.fastq -l 500 -q 10 > ${id}_trimmed_26.fq
		NanoFilt ${id}_27.fastq -l 500 -q 10 > ${id}_trimmed_27.fq
		NanoFilt ${id}_28.fastq -l 500 -q 10 > ${id}_trimmed_28.fq
		NanoFilt ${id}_29.fastq -l 500 -q 10 > ${id}_trimmed_29.fq
		NanoFilt ${id}_30.fastq -l 500 -q 10 > ${id}_trimmed_30.fq
		NanoFilt ${id}_31.fastq -l 500 -q 10 > ${id}_trimmed_31.fq
		NanoFilt ${id}_32.fastq -l 500 -q 10 > ${id}_trimmed_32.fq
		NanoFilt ${id}_33.fastq -l 500 -q 10 > ${id}_trimmed_33.fq
		NanoFilt ${id}_34.fastq -l 500 -q 10 > ${id}_trimmed_34.fq
		NanoFilt ${id}_35.fastq -l 500 -q 10 > ${id}_trimmed_35.fq
		NanoFilt ${id}_36.fastq -l 500 -q 10 > ${id}_trimmed_36.fq
		NanoFilt ${id}_37.fastq -l 500 -q 10 > ${id}_trimmed_37.fq
		NanoFilt ${id}_38.fastq -l 500 -q 10 > ${id}_trimmed_38.fq
		NanoFilt ${id}_39.fastq -l 500 -q 10 > ${id}_trimmed_39.fq
		NanoFilt ${id}_40.fastq -l 500 -q 10 > ${id}_trimmed_40.fq
		NanoFilt ${id}_41.fastq -l 500 -q 10 > ${id}_trimmed_41.fq
		NanoFilt ${id}_42.fastq -l 500 -q 10 > ${id}_trimmed_42.fq
		NanoFilt ${id}_43.fastq -l 500 -q 10 > ${id}_trimmed_43.fq
		NanoFilt ${id}_44.fastq -l 500 -q 10 > ${id}_trimmed_44.fq
		NanoFilt ${id}_45.fastq -l 500 -q 10 > ${id}_trimmed_45.fq
		NanoFilt ${id}_46.fastq -l 500 -q 10 > ${id}_trimmed_46.fq
		NanoFilt ${id}_47.fastq -l 500 -q 10 > ${id}_trimmed_47.fq
		NanoFilt ${id}_48.fastq -l 500 -q 10 > ${id}_trimmed_48.fq
		NanoFilt ${id}_49.fastq -l 500 -q 10 > ${id}_trimmed_49.fq
		NanoFilt ${id}_50.fastq -l 500 -q 10 > ${id}_trimmed_50.fq
		NanoFilt ${id}_51.fastq -l 500 -q 10 > ${id}_trimmed_51.fq
		NanoFilt ${id}_52.fastq -l 500 -q 10 > ${id}_trimmed_52.fq
		NanoFilt ${id}_53.fastq -l 500 -q 10 > ${id}_trimmed_53.fq
		NanoFilt ${id}_54.fastq -l 500 -q 10 > ${id}_trimmed_54.fq
		NanoFilt ${id}_55.fastq -l 500 -q 10 > ${id}_trimmed_55.fq
		NanoFilt ${id}_56.fastq -l 500 -q 10 > ${id}_trimmed_56.fq
		NanoFilt ${id}_57.fastq -l 500 -q 10 > ${id}_trimmed_57.fq
		NanoFilt ${id}_58.fastq -l 500 -q 10 > ${id}_trimmed_58.fq
		NanoFilt ${id}_59.fastq -l 500 -q 10 > ${id}_trimmed_59.fq
		NanoFilt ${id}_60.fastq -l 500 -q 10 > ${id}_trimmed_60.fq
		NanoFilt ${id}_61.fastq -l 500 -q 10 > ${id}_trimmed_61.fq
		NanoFilt ${id}_62.fastq -l 500 -q 10 > ${id}_trimmed_62.fq
		NanoFilt ${id}_63.fastq -l 500 -q 10 > ${id}_trimmed_63.fq
		NanoFilt ${id}_64.fastq -l 500 -q 10 > ${id}_trimmed_64.fq
		NanoFilt ${id}_65.fastq -l 500 -q 10 > ${id}_trimmed_65.fq
		NanoFilt ${id}_66.fastq -l 500 -q 10 > ${id}_trimmed_66.fq
		NanoFilt ${id}_67.fastq -l 500 -q 10 > ${id}_trimmed_67.fq
		NanoFilt ${id}_68.fastq -l 500 -q 10 > ${id}_trimmed_68.fq
		NanoFilt ${id}_69.fastq -l 500 -q 10 > ${id}_trimmed_69.fq
		NanoFilt ${id}_70.fastq -l 500 -q 10 > ${id}_trimmed_70.fq
		NanoFilt ${id}_71.fastq -l 500 -q 10 > ${id}_trimmed_71.fq
		NanoFilt ${id}_72.fastq -l 500 -q 10 > ${id}_trimmed_72.fq
		NanoFilt ${id}_73.fastq -l 500 -q 10 > ${id}_trimmed_73.fq
		NanoFilt ${id}_74.fastq -l 500 -q 10 > ${id}_trimmed_74.fq
		NanoFilt ${id}_75.fastq -l 500 -q 10 > ${id}_trimmed_75.fq
		NanoFilt ${id}_76.fastq -l 500 -q 10 > ${id}_trimmed_76.fq
		NanoFilt ${id}_77.fastq -l 500 -q 10 > ${id}_trimmed_77.fq
		NanoFilt ${id}_78.fastq -l 500 -q 10 > ${id}_trimmed_78.fq
		NanoFilt ${id}_79.fastq -l 500 -q 10 > ${id}_trimmed_79.fq
		NanoFilt ${id}_80.fastq -l 500 -q 10 > ${id}_trimmed_80.fq
		NanoFilt ${id}_81.fastq -l 500 -q 10 > ${id}_trimmed_81.fq
		NanoFilt ${id}_82.fastq -l 500 -q 10 > ${id}_trimmed_82.fq
		NanoFilt ${id}_83.fastq -l 500 -q 10 > ${id}_trimmed_83.fq
		NanoFilt ${id}_84.fastq -l 500 -q 10 > ${id}_trimmed_84.fq
		NanoFilt ${id}_85.fastq -l 500 -q 10 > ${id}_trimmed_85.fq
		NanoFilt ${id}_86.fastq -l 500 -q 10 > ${id}_trimmed_86.fq
		NanoFilt ${id}_87.fastq -l 500 -q 10 > ${id}_trimmed_87.fq
		NanoFilt ${id}_88.fastq -l 500 -q 10 > ${id}_trimmed_88.fq
		NanoFilt ${id}_89.fastq -l 500 -q 10 > ${id}_trimmed_89.fq
		NanoFilt ${id}_90.fastq -l 500 -q 10 > ${id}_trimmed_90.fq
		NanoFilt ${id}_91.fastq -l 500 -q 10 > ${id}_trimmed_91.fq
		NanoFilt ${id}_92.fastq -l 500 -q 10 > ${id}_trimmed_92.fq
		NanoFilt ${id}_93.fastq -l 500 -q 10 > ${id}_trimmed_93.fq
		NanoFilt ${id}_94.fastq -l 500 -q 10 > ${id}_trimmed_94.fq
		NanoFilt ${id}_95.fastq -l 500 -q 10 > ${id}_trimmed_95.fq
		NanoFilt ${id}_96.fastq -l 500 -q 10 > ${id}_trimmed_96.fq
		NanoFilt ${id}_97.fastq -l 500 -q 10 > ${id}_trimmed_97.fq
		NanoFilt ${id}_98.fastq -l 500 -q 10 > ${id}_trimmed_98.fq
		NanoFilt ${id}_99.fastq -l 500 -q 10 > ${id}_trimmed_99.fq
		NanoFilt ${id}_100.fastq -l 500 -q 10 > ${id}_trimmed_100.fq
		NanoFilt ${id}_101.fastq -l 500 -q 10 > ${id}_trimmed_101.fq
		NanoFilt ${id}_102.fastq -l 500 -q 10 > ${id}_trimmed_102.fq
		NanoFilt ${id}_103.fastq -l 500 -q 10 > ${id}_trimmed_103.fq
		NanoFilt ${id}_104.fastq -l 500 -q 10 > ${id}_trimmed_104.fq
		NanoFilt ${id}_105.fastq -l 500 -q 10 > ${id}_trimmed_105.fq
		NanoFilt ${id}_106.fastq -l 500 -q 10 > ${id}_trimmed_106.fq
		NanoFilt ${id}_107.fastq -l 500 -q 10 > ${id}_trimmed_107.fq
		NanoFilt ${id}_108.fastq -l 500 -q 10 > ${id}_trimmed_108.fq
		NanoFilt ${id}_109.fastq -l 500 -q 10 > ${id}_trimmed_109.fq
		NanoFilt ${id}_110.fastq -l 500 -q 10 > ${id}_trimmed_110.fq
		NanoFilt ${id}_111.fastq -l 500 -q 10 > ${id}_trimmed_111.fq
		NanoFilt ${id}_112.fastq -l 500 -q 10 > ${id}_trimmed_112.fq
		NanoFilt ${id}_113.fastq -l 500 -q 10 > ${id}_trimmed_113.fq
		NanoFilt ${id}_114.fastq -l 500 -q 10 > ${id}_trimmed_114.fq
		NanoFilt ${id}_115.fastq -l 500 -q 10 > ${id}_trimmed_115.fq
		NanoFilt ${id}_116.fastq -l 500 -q 10 > ${id}_trimmed_116.fq
		NanoFilt ${id}_117.fastq -l 500 -q 10 > ${id}_trimmed_117.fq
		NanoFilt ${id}_118.fastq -l 500 -q 10 > ${id}_trimmed_118.fq
		NanoFilt ${id}_119.fastq -l 500 -q 10 > ${id}_trimmed_119.fq
		NanoFilt ${id}_120.fastq -l 500 -q 10 > ${id}_trimmed_120.fq
		NanoFilt ${id}_121.fastq -l 500 -q 10 > ${id}_trimmed_121.fq
		NanoFilt ${id}_122.fastq -l 500 -q 10 > ${id}_trimmed_122.fq
		NanoFilt ${id}_123.fastq -l 500 -q 10 > ${id}_trimmed_123.fq
		NanoFilt ${id}_124.fastq -l 500 -q 10 > ${id}_trimmed_124.fq
		NanoFilt ${id}_125.fastq -l 500 -q 10 > ${id}_trimmed_125.fq
		NanoFilt ${id}_126.fastq -l 500 -q 10 > ${id}_trimmed_126.fq
		NanoFilt ${id}_127.fastq -l 500 -q 10 > ${id}_trimmed_127.fq
		NanoFilt ${id}_128.fastq -l 500 -q 10 > ${id}_trimmed_128.fq
		NanoFilt ${id}_129.fastq -l 500 -q 10 > ${id}_trimmed_129.fq
		NanoFilt ${id}_130.fastq -l 500 -q 10 > ${id}_trimmed_130.fq
		NanoFilt ${id}_131.fastq -l 500 -q 10 > ${id}_trimmed_131.fq
		NanoFilt ${id}_132.fastq -l 500 -q 10 > ${id}_trimmed_132.fq
		NanoFilt ${id}_133.fastq -l 500 -q 10 > ${id}_trimmed_133.fq
		NanoFilt ${id}_134.fastq -l 500 -q 10 > ${id}_trimmed_134.fq
		NanoFilt ${id}_135.fastq -l 500 -q 10 > ${id}_trimmed_135.fq
		NanoFilt ${id}_136.fastq -l 500 -q 10 > ${id}_trimmed_136.fq
		NanoFilt ${id}_137.fastq -l 500 -q 10 > ${id}_trimmed_137.fq
		NanoFilt ${id}_138.fastq -l 500 -q 10 > ${id}_trimmed_138.fq
		NanoFilt ${id}_139.fastq -l 500 -q 10 > ${id}_trimmed_139.fq
		NanoFilt ${id}_140.fastq -l 500 -q 10 > ${id}_trimmed_140.fq
		NanoFilt ${id}_141.fastq -l 500 -q 10 > ${id}_trimmed_141.fq
		NanoFilt ${id}_142.fastq -l 500 -q 10 > ${id}_trimmed_142.fq
		NanoFilt ${id}_143.fastq -l 500 -q 10 > ${id}_trimmed_143.fq
		NanoFilt ${id}_144.fastq -l 500 -q 10 > ${id}_trimmed_144.fq
		NanoFilt ${id}_145.fastq -l 500 -q 10 > ${id}_trimmed_145.fq
		NanoFilt ${id}_146.fastq -l 500 -q 10 > ${id}_trimmed_146.fq
		NanoFilt ${id}_147.fastq -l 500 -q 10 > ${id}_trimmed_147.fq
		NanoFilt ${id}_148.fastq -l 500 -q 10 > ${id}_trimmed_148.fq
		NanoFilt ${id}_149.fastq -l 500 -q 10 > ${id}_trimmed_149.fq
		NanoFilt ${id}_150.fastq -l 500 -q 10 > ${id}_trimmed_150.fq
		NanoFilt ${id}_151.fastq -l 500 -q 10 > ${id}_trimmed_151.fq
		NanoFilt ${id}_152.fastq -l 500 -q 10 > ${id}_trimmed_152.fq
		NanoFilt ${id}_153.fastq -l 500 -q 10 > ${id}_trimmed_153.fq
		NanoFilt ${id}_154.fastq -l 500 -q 10 > ${id}_trimmed_154.fq
		NanoFilt ${id}_155.fastq -l 500 -q 10 > ${id}_trimmed_155.fq
		NanoFilt ${id}_156.fastq -l 500 -q 10 > ${id}_trimmed_156.fq
		NanoFilt ${id}_157.fastq -l 500 -q 10 > ${id}_trimmed_157.fq
		NanoFilt ${id}_158.fastq -l 500 -q 10 > ${id}_trimmed_158.fq
		NanoFilt ${id}_159.fastq -l 500 -q 10 > ${id}_trimmed_159.fq
		NanoFilt ${id}_160.fastq -l 500 -q 10 > ${id}_trimmed_160.fq
		NanoFilt ${id}_161.fastq -l 500 -q 10 > ${id}_trimmed_161.fq
		NanoFilt ${id}_162.fastq -l 500 -q 10 > ${id}_trimmed_162.fq
		NanoFilt ${id}_163.fastq -l 500 -q 10 > ${id}_trimmed_163.fq
		NanoFilt ${id}_164.fastq -l 500 -q 10 > ${id}_trimmed_164.fq
		NanoFilt ${id}_165.fastq -l 500 -q 10 > ${id}_trimmed_165.fq
		NanoFilt ${id}_166.fastq -l 500 -q 10 > ${id}_trimmed_166.fq
		NanoFilt ${id}_167.fastq -l 500 -q 10 > ${id}_trimmed_167.fq
		NanoFilt ${id}_168.fastq -l 500 -q 10 > ${id}_trimmed_168.fq
		NanoFilt ${id}_169.fastq -l 500 -q 10 > ${id}_trimmed_169.fq
		NanoFilt ${id}_170.fastq -l 500 -q 10 > ${id}_trimmed_170.fq
		NanoFilt ${id}_171.fastq -l 500 -q 10 > ${id}_trimmed_171.fq
		NanoFilt ${id}_172.fastq -l 500 -q 10 > ${id}_trimmed_172.fq
		NanoFilt ${id}_173.fastq -l 500 -q 10 > ${id}_trimmed_173.fq
		NanoFilt ${id}_174.fastq -l 500 -q 10 > ${id}_trimmed_174.fq
		NanoFilt ${id}_175.fastq -l 500 -q 10 > ${id}_trimmed_175.fq
		NanoFilt ${id}_176.fastq -l 500 -q 10 > ${id}_trimmed_176.fq
		NanoFilt ${id}_177.fastq -l 500 -q 10 > ${id}_trimmed_177.fq
		NanoFilt ${id}_178.fastq -l 500 -q 10 > ${id}_trimmed_178.fq
		NanoFilt ${id}_179.fastq -l 500 -q 10 > ${id}_trimmed_179.fq
		NanoFilt ${id}_180.fastq -l 500 -q 10 > ${id}_trimmed_180.fq
		NanoFilt ${id}_181.fastq -l 500 -q 10 > ${id}_trimmed_181.fq
		NanoFilt ${id}_182.fastq -l 500 -q 10 > ${id}_trimmed_182.fq
		"""
        




}

process maptoreference32 {
    publishDir "$params.output.folder32/mappedreads/${id}", pattern: "*.fq", mode : "copy"
    publishDir "$params.output.folder32/unmappedreads/${id}", pattern: "*.fq", mode : "copy"

    input:
		set val(id), path(trim_read0) from trim_out1125899906842624
		set val(id), path(trim_read1) from trim_out2251799813685248
		set val(id), path(trim_read2) from trim_out3377699720527872
		set val(id), path(trim_read3) from trim_out4503599627370496
		set val(id), path(trim_read4) from trim_out5629499534213120
		set val(id), path(trim_read5) from trim_out6755399441055744
		set val(id), path(trim_read6) from trim_out7881299347898368
		set val(id), path(trim_read7) from trim_out9007199254740992
		set val(id), path(trim_read8) from trim_out10133099161583616
		set val(id), path(trim_read9) from trim_out11258999068426240
		set val(id), path(trim_read10) from trim_out12384898975268864
		set val(id), path(trim_read11) from trim_out13510798882111488
		set val(id), path(trim_read12) from trim_out14636698788954112
		set val(id), path(trim_read13) from trim_out15762598695796736
		set val(id), path(trim_read14) from trim_out16888498602639360
		set val(id), path(trim_read15) from trim_out18014398509481984
		set val(id), path(trim_read16) from trim_out19140298416324608
		set val(id), path(trim_read17) from trim_out20266198323167232
		set val(id), path(trim_read18) from trim_out21392098230009856
		set val(id), path(trim_read19) from trim_out22517998136852480
		set val(id), path(trim_read20) from trim_out23643898043695104
		set val(id), path(trim_read21) from trim_out24769797950537728
		set val(id), path(trim_read22) from trim_out25895697857380352
		set val(id), path(trim_read23) from trim_out27021597764222976
		set val(id), path(trim_read24) from trim_out28147497671065600
		set val(id), path(trim_read25) from trim_out29273397577908224
		set val(id), path(trim_read26) from trim_out30399297484750848
		set val(id), path(trim_read27) from trim_out31525197391593472
		set val(id), path(trim_read28) from trim_out32651097298436096
		set val(id), path(trim_read29) from trim_out33776997205278720
		set val(id), path(trim_read30) from trim_out34902897112121344
		set val(id), path(trim_read31) from trim_out36028797018963968
		set val(id), path(trim_read32) from trim_out37154696925806592
		set val(id), path(trim_read33) from trim_out38280596832649216
		set val(id), path(trim_read34) from trim_out39406496739491840
		set val(id), path(trim_read35) from trim_out40532396646334464
		set val(id), path(trim_read36) from trim_out41658296553177088
		set val(id), path(trim_read37) from trim_out42784196460019712
		set val(id), path(trim_read38) from trim_out43910096366862336
		set val(id), path(trim_read39) from trim_out45035996273704960
		set val(id), path(trim_read40) from trim_out46161896180547584
		set val(id), path(trim_read41) from trim_out47287796087390208
		set val(id), path(trim_read42) from trim_out48413695994232832
		set val(id), path(trim_read43) from trim_out49539595901075456
		set val(id), path(trim_read44) from trim_out50665495807918080
		set val(id), path(trim_read45) from trim_out51791395714760704
		set val(id), path(trim_read46) from trim_out52917295621603328
		set val(id), path(trim_read47) from trim_out54043195528445952
		set val(id), path(trim_read48) from trim_out55169095435288576
		set val(id), path(trim_read49) from trim_out56294995342131200
		set val(id), path(trim_read50) from trim_out57420895248973824
		set val(id), path(trim_read51) from trim_out58546795155816448
		set val(id), path(trim_read52) from trim_out59672695062659072
		set val(id), path(trim_read53) from trim_out60798594969501696
		set val(id), path(trim_read54) from trim_out61924494876344320
		set val(id), path(trim_read55) from trim_out63050394783186944
		set val(id), path(trim_read56) from trim_out64176294690029568
		set val(id), path(trim_read57) from trim_out65302194596872192
		set val(id), path(trim_read58) from trim_out66428094503714816
		set val(id), path(trim_read59) from trim_out67553994410557440
		set val(id), path(trim_read60) from trim_out68679894317400064
		set val(id), path(trim_read61) from trim_out69805794224242688
		set val(id), path(trim_read62) from trim_out70931694131085312
		set val(id), path(trim_read63) from trim_out72057594037927936
		set val(id), path(trim_read64) from trim_out73183493944770560
		set val(id), path(trim_read65) from trim_out74309393851613184
		set val(id), path(trim_read66) from trim_out75435293758455808
		set val(id), path(trim_read67) from trim_out76561193665298432
		set val(id), path(trim_read68) from trim_out77687093572141056
		set val(id), path(trim_read69) from trim_out78812993478983680
		set val(id), path(trim_read70) from trim_out79938893385826304
		set val(id), path(trim_read71) from trim_out81064793292668928
		set val(id), path(trim_read72) from trim_out82190693199511552
		set val(id), path(trim_read73) from trim_out83316593106354176
		set val(id), path(trim_read74) from trim_out84442493013196800
		set val(id), path(trim_read75) from trim_out85568392920039424
		set val(id), path(trim_read76) from trim_out86694292826882048
		set val(id), path(trim_read77) from trim_out87820192733724672
		set val(id), path(trim_read78) from trim_out88946092640567296
		set val(id), path(trim_read79) from trim_out90071992547409920
		set val(id), path(trim_read80) from trim_out91197892454252544
		set val(id), path(trim_read81) from trim_out92323792361095168
		set val(id), path(trim_read82) from trim_out93449692267937792
		set val(id), path(trim_read83) from trim_out94575592174780416
		set val(id), path(trim_read84) from trim_out95701492081623040
		set val(id), path(trim_read85) from trim_out96827391988465664
		set val(id), path(trim_read86) from trim_out97953291895308288
		set val(id), path(trim_read87) from trim_out99079191802150912
		set val(id), path(trim_read88) from trim_out100205091708993536
		set val(id), path(trim_read89) from trim_out101330991615836160
		set val(id), path(trim_read90) from trim_out102456891522678784
		set val(id), path(trim_read91) from trim_out103582791429521408
		set val(id), path(trim_read92) from trim_out104708691336364032
		set val(id), path(trim_read93) from trim_out105834591243206656
		set val(id), path(trim_read94) from trim_out106960491150049280
		set val(id), path(trim_read95) from trim_out108086391056891904
		set val(id), path(trim_read96) from trim_out109212290963734528
		set val(id), path(trim_read97) from trim_out110338190870577152
		set val(id), path(trim_read98) from trim_out111464090777419776
		set val(id), path(trim_read99) from trim_out112589990684262400
		set val(id), path(trim_read100) from trim_out113715890591105024
		set val(id), path(trim_read101) from trim_out114841790497947648
		set val(id), path(trim_read102) from trim_out115967690404790272
		set val(id), path(trim_read103) from trim_out117093590311632896
		set val(id), path(trim_read104) from trim_out118219490218475520
		set val(id), path(trim_read105) from trim_out119345390125318144
		set val(id), path(trim_read106) from trim_out120471290032160768
		set val(id), path(trim_read107) from trim_out121597189939003392
		set val(id), path(trim_read108) from trim_out122723089845846016
		set val(id), path(trim_read109) from trim_out123848989752688640
		set val(id), path(trim_read110) from trim_out124974889659531264
		set val(id), path(trim_read111) from trim_out126100789566373888
		set val(id), path(trim_read112) from trim_out127226689473216512
		set val(id), path(trim_read113) from trim_out128352589380059136
		set val(id), path(trim_read114) from trim_out129478489286901760
		set val(id), path(trim_read115) from trim_out130604389193744384
		set val(id), path(trim_read116) from trim_out131730289100587008
		set val(id), path(trim_read117) from trim_out132856189007429632
		set val(id), path(trim_read118) from trim_out133982088914272256
		set val(id), path(trim_read119) from trim_out135107988821114880
		set val(id), path(trim_read120) from trim_out136233888727957504
		set val(id), path(trim_read121) from trim_out137359788634800128
		set val(id), path(trim_read122) from trim_out138485688541642752
		set val(id), path(trim_read123) from trim_out139611588448485376
		set val(id), path(trim_read124) from trim_out140737488355328000
		set val(id), path(trim_read125) from trim_out141863388262170624
		set val(id), path(trim_read126) from trim_out142989288169013248
		set val(id), path(trim_read127) from trim_out144115188075855872
		set val(id), path(trim_read128) from trim_out145241087982698496
		set val(id), path(trim_read129) from trim_out146366987889541120
		set val(id), path(trim_read130) from trim_out147492887796383744
		set val(id), path(trim_read131) from trim_out148618787703226368
		set val(id), path(trim_read132) from trim_out149744687610068992
		set val(id), path(trim_read133) from trim_out150870587516911616
		set val(id), path(trim_read134) from trim_out151996487423754240
		set val(id), path(trim_read135) from trim_out153122387330596864
		set val(id), path(trim_read136) from trim_out154248287237439488
		set val(id), path(trim_read137) from trim_out155374187144282112
		set val(id), path(trim_read138) from trim_out156500087051124736
		set val(id), path(trim_read139) from trim_out157625986957967360
		set val(id), path(trim_read140) from trim_out158751886864809984
		set val(id), path(trim_read141) from trim_out159877786771652608
		set val(id), path(trim_read142) from trim_out161003686678495232
		set val(id), path(trim_read143) from trim_out162129586585337856
		set val(id), path(trim_read144) from trim_out163255486492180480
		set val(id), path(trim_read145) from trim_out164381386399023104
		set val(id), path(trim_read146) from trim_out165507286305865728
		set val(id), path(trim_read147) from trim_out166633186212708352
		set val(id), path(trim_read148) from trim_out167759086119550976
		set val(id), path(trim_read149) from trim_out168884986026393600
		set val(id), path(trim_read150) from trim_out170010885933236224
		set val(id), path(trim_read151) from trim_out171136785840078848
		set val(id), path(trim_read152) from trim_out172262685746921472
		set val(id), path(trim_read153) from trim_out173388585653764096
		set val(id), path(trim_read154) from trim_out174514485560606720
		set val(id), path(trim_read155) from trim_out175640385467449344
		set val(id), path(trim_read156) from trim_out176766285374291968
		set val(id), path(trim_read157) from trim_out177892185281134592
		set val(id), path(trim_read158) from trim_out179018085187977216
		set val(id), path(trim_read159) from trim_out180143985094819840
		set val(id), path(trim_read160) from trim_out181269885001662464
		set val(id), path(trim_read161) from trim_out182395784908505088
		set val(id), path(trim_read162) from trim_out183521684815347712
		set val(id), path(trim_read163) from trim_out184647584722190336
		set val(id), path(trim_read164) from trim_out185773484629032960
		set val(id), path(trim_read165) from trim_out186899384535875584
		set val(id), path(trim_read166) from trim_out188025284442718208
		set val(id), path(trim_read167) from trim_out189151184349560832
		set val(id), path(trim_read168) from trim_out190277084256403456
		set val(id), path(trim_read169) from trim_out191402984163246080
		set val(id), path(trim_read170) from trim_out192528884070088704
		set val(id), path(trim_read171) from trim_out193654783976931328
		set val(id), path(trim_read172) from trim_out194780683883773952
		set val(id), path(trim_read173) from trim_out195906583790616576
		set val(id), path(trim_read174) from trim_out197032483697459200
		set val(id), path(trim_read175) from trim_out198158383604301824
		set val(id), path(trim_read176) from trim_out199284283511144448
		set val(id), path(trim_read177) from trim_out200410183417987072
		set val(id), path(trim_read178) from trim_out201536083324829696
		set val(id), path(trim_read179) from trim_out202661983231672320
		set val(id), path(trim_read180) from trim_out203787883138514944
		set val(id), path(trim_read181) from trim_out204913783045357568
		set val(id), path(trim_read182) from trim_out206039682952200192

       
    output:
		set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1125899906842624
		set val(id), path("${id}_mapped_1.fq"), path("${id}_unmapped_1.fq") into mapped_out2251799813685248
		set val(id), path("${id}_mapped_2.fq"), path("${id}_unmapped_2.fq") into mapped_out3377699720527872
		set val(id), path("${id}_mapped_3.fq"), path("${id}_unmapped_3.fq") into mapped_out4503599627370496
		set val(id), path("${id}_mapped_4.fq"), path("${id}_unmapped_4.fq") into mapped_out5629499534213120
		set val(id), path("${id}_mapped_5.fq"), path("${id}_unmapped_5.fq") into mapped_out6755399441055744
		set val(id), path("${id}_mapped_6.fq"), path("${id}_unmapped_6.fq") into mapped_out7881299347898368
		set val(id), path("${id}_mapped_7.fq"), path("${id}_unmapped_7.fq") into mapped_out9007199254740992
		set val(id), path("${id}_mapped_8.fq"), path("${id}_unmapped_8.fq") into mapped_out10133099161583616
		set val(id), path("${id}_mapped_9.fq"), path("${id}_unmapped_9.fq") into mapped_out11258999068426240
		set val(id), path("${id}_mapped_10.fq"), path("${id}_unmapped_10.fq") into mapped_out12384898975268864
		set val(id), path("${id}_mapped_11.fq"), path("${id}_unmapped_11.fq") into mapped_out13510798882111488
		set val(id), path("${id}_mapped_12.fq"), path("${id}_unmapped_12.fq") into mapped_out14636698788954112
		set val(id), path("${id}_mapped_13.fq"), path("${id}_unmapped_13.fq") into mapped_out15762598695796736
		set val(id), path("${id}_mapped_14.fq"), path("${id}_unmapped_14.fq") into mapped_out16888498602639360
		set val(id), path("${id}_mapped_15.fq"), path("${id}_unmapped_15.fq") into mapped_out18014398509481984
		set val(id), path("${id}_mapped_16.fq"), path("${id}_unmapped_16.fq") into mapped_out19140298416324608
		set val(id), path("${id}_mapped_17.fq"), path("${id}_unmapped_17.fq") into mapped_out20266198323167232
		set val(id), path("${id}_mapped_18.fq"), path("${id}_unmapped_18.fq") into mapped_out21392098230009856
		set val(id), path("${id}_mapped_19.fq"), path("${id}_unmapped_19.fq") into mapped_out22517998136852480
		set val(id), path("${id}_mapped_20.fq"), path("${id}_unmapped_20.fq") into mapped_out23643898043695104
		set val(id), path("${id}_mapped_21.fq"), path("${id}_unmapped_21.fq") into mapped_out24769797950537728
		set val(id), path("${id}_mapped_22.fq"), path("${id}_unmapped_22.fq") into mapped_out25895697857380352
		set val(id), path("${id}_mapped_23.fq"), path("${id}_unmapped_23.fq") into mapped_out27021597764222976
		set val(id), path("${id}_mapped_24.fq"), path("${id}_unmapped_24.fq") into mapped_out28147497671065600
		set val(id), path("${id}_mapped_25.fq"), path("${id}_unmapped_25.fq") into mapped_out29273397577908224
		set val(id), path("${id}_mapped_26.fq"), path("${id}_unmapped_26.fq") into mapped_out30399297484750848
		set val(id), path("${id}_mapped_27.fq"), path("${id}_unmapped_27.fq") into mapped_out31525197391593472
		set val(id), path("${id}_mapped_28.fq"), path("${id}_unmapped_28.fq") into mapped_out32651097298436096
		set val(id), path("${id}_mapped_29.fq"), path("${id}_unmapped_29.fq") into mapped_out33776997205278720
		set val(id), path("${id}_mapped_30.fq"), path("${id}_unmapped_30.fq") into mapped_out34902897112121344
		set val(id), path("${id}_mapped_31.fq"), path("${id}_unmapped_31.fq") into mapped_out36028797018963968
		set val(id), path("${id}_mapped_32.fq"), path("${id}_unmapped_32.fq") into mapped_out37154696925806592
		set val(id), path("${id}_mapped_33.fq"), path("${id}_unmapped_33.fq") into mapped_out38280596832649216
		set val(id), path("${id}_mapped_34.fq"), path("${id}_unmapped_34.fq") into mapped_out39406496739491840
		set val(id), path("${id}_mapped_35.fq"), path("${id}_unmapped_35.fq") into mapped_out40532396646334464
		set val(id), path("${id}_mapped_36.fq"), path("${id}_unmapped_36.fq") into mapped_out41658296553177088
		set val(id), path("${id}_mapped_37.fq"), path("${id}_unmapped_37.fq") into mapped_out42784196460019712
		set val(id), path("${id}_mapped_38.fq"), path("${id}_unmapped_38.fq") into mapped_out43910096366862336
		set val(id), path("${id}_mapped_39.fq"), path("${id}_unmapped_39.fq") into mapped_out45035996273704960
		set val(id), path("${id}_mapped_40.fq"), path("${id}_unmapped_40.fq") into mapped_out46161896180547584
		set val(id), path("${id}_mapped_41.fq"), path("${id}_unmapped_41.fq") into mapped_out47287796087390208
		set val(id), path("${id}_mapped_42.fq"), path("${id}_unmapped_42.fq") into mapped_out48413695994232832
		set val(id), path("${id}_mapped_43.fq"), path("${id}_unmapped_43.fq") into mapped_out49539595901075456
		set val(id), path("${id}_mapped_44.fq"), path("${id}_unmapped_44.fq") into mapped_out50665495807918080
		set val(id), path("${id}_mapped_45.fq"), path("${id}_unmapped_45.fq") into mapped_out51791395714760704
		set val(id), path("${id}_mapped_46.fq"), path("${id}_unmapped_46.fq") into mapped_out52917295621603328
		set val(id), path("${id}_mapped_47.fq"), path("${id}_unmapped_47.fq") into mapped_out54043195528445952
		set val(id), path("${id}_mapped_48.fq"), path("${id}_unmapped_48.fq") into mapped_out55169095435288576
		set val(id), path("${id}_mapped_49.fq"), path("${id}_unmapped_49.fq") into mapped_out56294995342131200
		set val(id), path("${id}_mapped_50.fq"), path("${id}_unmapped_50.fq") into mapped_out57420895248973824
		set val(id), path("${id}_mapped_51.fq"), path("${id}_unmapped_51.fq") into mapped_out58546795155816448
		set val(id), path("${id}_mapped_52.fq"), path("${id}_unmapped_52.fq") into mapped_out59672695062659072
		set val(id), path("${id}_mapped_53.fq"), path("${id}_unmapped_53.fq") into mapped_out60798594969501696
		set val(id), path("${id}_mapped_54.fq"), path("${id}_unmapped_54.fq") into mapped_out61924494876344320
		set val(id), path("${id}_mapped_55.fq"), path("${id}_unmapped_55.fq") into mapped_out63050394783186944
		set val(id), path("${id}_mapped_56.fq"), path("${id}_unmapped_56.fq") into mapped_out64176294690029568
		set val(id), path("${id}_mapped_57.fq"), path("${id}_unmapped_57.fq") into mapped_out65302194596872192
		set val(id), path("${id}_mapped_58.fq"), path("${id}_unmapped_58.fq") into mapped_out66428094503714816
		set val(id), path("${id}_mapped_59.fq"), path("${id}_unmapped_59.fq") into mapped_out67553994410557440
		set val(id), path("${id}_mapped_60.fq"), path("${id}_unmapped_60.fq") into mapped_out68679894317400064
		set val(id), path("${id}_mapped_61.fq"), path("${id}_unmapped_61.fq") into mapped_out69805794224242688
		set val(id), path("${id}_mapped_62.fq"), path("${id}_unmapped_62.fq") into mapped_out70931694131085312
		set val(id), path("${id}_mapped_63.fq"), path("${id}_unmapped_63.fq") into mapped_out72057594037927936
		set val(id), path("${id}_mapped_64.fq"), path("${id}_unmapped_64.fq") into mapped_out73183493944770560
		set val(id), path("${id}_mapped_65.fq"), path("${id}_unmapped_65.fq") into mapped_out74309393851613184
		set val(id), path("${id}_mapped_66.fq"), path("${id}_unmapped_66.fq") into mapped_out75435293758455808
		set val(id), path("${id}_mapped_67.fq"), path("${id}_unmapped_67.fq") into mapped_out76561193665298432
		set val(id), path("${id}_mapped_68.fq"), path("${id}_unmapped_68.fq") into mapped_out77687093572141056
		set val(id), path("${id}_mapped_69.fq"), path("${id}_unmapped_69.fq") into mapped_out78812993478983680
		set val(id), path("${id}_mapped_70.fq"), path("${id}_unmapped_70.fq") into mapped_out79938893385826304
		set val(id), path("${id}_mapped_71.fq"), path("${id}_unmapped_71.fq") into mapped_out81064793292668928
		set val(id), path("${id}_mapped_72.fq"), path("${id}_unmapped_72.fq") into mapped_out82190693199511552
		set val(id), path("${id}_mapped_73.fq"), path("${id}_unmapped_73.fq") into mapped_out83316593106354176
		set val(id), path("${id}_mapped_74.fq"), path("${id}_unmapped_74.fq") into mapped_out84442493013196800
		set val(id), path("${id}_mapped_75.fq"), path("${id}_unmapped_75.fq") into mapped_out85568392920039424
		set val(id), path("${id}_mapped_76.fq"), path("${id}_unmapped_76.fq") into mapped_out86694292826882048
		set val(id), path("${id}_mapped_77.fq"), path("${id}_unmapped_77.fq") into mapped_out87820192733724672
		set val(id), path("${id}_mapped_78.fq"), path("${id}_unmapped_78.fq") into mapped_out88946092640567296
		set val(id), path("${id}_mapped_79.fq"), path("${id}_unmapped_79.fq") into mapped_out90071992547409920
		set val(id), path("${id}_mapped_80.fq"), path("${id}_unmapped_80.fq") into mapped_out91197892454252544
		set val(id), path("${id}_mapped_81.fq"), path("${id}_unmapped_81.fq") into mapped_out92323792361095168
		set val(id), path("${id}_mapped_82.fq"), path("${id}_unmapped_82.fq") into mapped_out93449692267937792
		set val(id), path("${id}_mapped_83.fq"), path("${id}_unmapped_83.fq") into mapped_out94575592174780416
		set val(id), path("${id}_mapped_84.fq"), path("${id}_unmapped_84.fq") into mapped_out95701492081623040
		set val(id), path("${id}_mapped_85.fq"), path("${id}_unmapped_85.fq") into mapped_out96827391988465664
		set val(id), path("${id}_mapped_86.fq"), path("${id}_unmapped_86.fq") into mapped_out97953291895308288
		set val(id), path("${id}_mapped_87.fq"), path("${id}_unmapped_87.fq") into mapped_out99079191802150912
		set val(id), path("${id}_mapped_88.fq"), path("${id}_unmapped_88.fq") into mapped_out100205091708993536
		set val(id), path("${id}_mapped_89.fq"), path("${id}_unmapped_89.fq") into mapped_out101330991615836160
		set val(id), path("${id}_mapped_90.fq"), path("${id}_unmapped_90.fq") into mapped_out102456891522678784
		set val(id), path("${id}_mapped_91.fq"), path("${id}_unmapped_91.fq") into mapped_out103582791429521408
		set val(id), path("${id}_mapped_92.fq"), path("${id}_unmapped_92.fq") into mapped_out104708691336364032
		set val(id), path("${id}_mapped_93.fq"), path("${id}_unmapped_93.fq") into mapped_out105834591243206656
		set val(id), path("${id}_mapped_94.fq"), path("${id}_unmapped_94.fq") into mapped_out106960491150049280
		set val(id), path("${id}_mapped_95.fq"), path("${id}_unmapped_95.fq") into mapped_out108086391056891904
		set val(id), path("${id}_mapped_96.fq"), path("${id}_unmapped_96.fq") into mapped_out109212290963734528
		set val(id), path("${id}_mapped_97.fq"), path("${id}_unmapped_97.fq") into mapped_out110338190870577152
		set val(id), path("${id}_mapped_98.fq"), path("${id}_unmapped_98.fq") into mapped_out111464090777419776
		set val(id), path("${id}_mapped_99.fq"), path("${id}_unmapped_99.fq") into mapped_out112589990684262400
		set val(id), path("${id}_mapped_100.fq"), path("${id}_unmapped_100.fq") into mapped_out113715890591105024
		set val(id), path("${id}_mapped_101.fq"), path("${id}_unmapped_101.fq") into mapped_out114841790497947648
		set val(id), path("${id}_mapped_102.fq"), path("${id}_unmapped_102.fq") into mapped_out115967690404790272
		set val(id), path("${id}_mapped_103.fq"), path("${id}_unmapped_103.fq") into mapped_out117093590311632896
		set val(id), path("${id}_mapped_104.fq"), path("${id}_unmapped_104.fq") into mapped_out118219490218475520
		set val(id), path("${id}_mapped_105.fq"), path("${id}_unmapped_105.fq") into mapped_out119345390125318144
		set val(id), path("${id}_mapped_106.fq"), path("${id}_unmapped_106.fq") into mapped_out120471290032160768
		set val(id), path("${id}_mapped_107.fq"), path("${id}_unmapped_107.fq") into mapped_out121597189939003392
		set val(id), path("${id}_mapped_108.fq"), path("${id}_unmapped_108.fq") into mapped_out122723089845846016
		set val(id), path("${id}_mapped_109.fq"), path("${id}_unmapped_109.fq") into mapped_out123848989752688640
		set val(id), path("${id}_mapped_110.fq"), path("${id}_unmapped_110.fq") into mapped_out124974889659531264
		set val(id), path("${id}_mapped_111.fq"), path("${id}_unmapped_111.fq") into mapped_out126100789566373888
		set val(id), path("${id}_mapped_112.fq"), path("${id}_unmapped_112.fq") into mapped_out127226689473216512
		set val(id), path("${id}_mapped_113.fq"), path("${id}_unmapped_113.fq") into mapped_out128352589380059136
		set val(id), path("${id}_mapped_114.fq"), path("${id}_unmapped_114.fq") into mapped_out129478489286901760
		set val(id), path("${id}_mapped_115.fq"), path("${id}_unmapped_115.fq") into mapped_out130604389193744384
		set val(id), path("${id}_mapped_116.fq"), path("${id}_unmapped_116.fq") into mapped_out131730289100587008
		set val(id), path("${id}_mapped_117.fq"), path("${id}_unmapped_117.fq") into mapped_out132856189007429632
		set val(id), path("${id}_mapped_118.fq"), path("${id}_unmapped_118.fq") into mapped_out133982088914272256
		set val(id), path("${id}_mapped_119.fq"), path("${id}_unmapped_119.fq") into mapped_out135107988821114880
		set val(id), path("${id}_mapped_120.fq"), path("${id}_unmapped_120.fq") into mapped_out136233888727957504
		set val(id), path("${id}_mapped_121.fq"), path("${id}_unmapped_121.fq") into mapped_out137359788634800128
		set val(id), path("${id}_mapped_122.fq"), path("${id}_unmapped_122.fq") into mapped_out138485688541642752
		set val(id), path("${id}_mapped_123.fq"), path("${id}_unmapped_123.fq") into mapped_out139611588448485376
		set val(id), path("${id}_mapped_124.fq"), path("${id}_unmapped_124.fq") into mapped_out140737488355328000
		set val(id), path("${id}_mapped_125.fq"), path("${id}_unmapped_125.fq") into mapped_out141863388262170624
		set val(id), path("${id}_mapped_126.fq"), path("${id}_unmapped_126.fq") into mapped_out142989288169013248
		set val(id), path("${id}_mapped_127.fq"), path("${id}_unmapped_127.fq") into mapped_out144115188075855872
		set val(id), path("${id}_mapped_128.fq"), path("${id}_unmapped_128.fq") into mapped_out145241087982698496
		set val(id), path("${id}_mapped_129.fq"), path("${id}_unmapped_129.fq") into mapped_out146366987889541120
		set val(id), path("${id}_mapped_130.fq"), path("${id}_unmapped_130.fq") into mapped_out147492887796383744
		set val(id), path("${id}_mapped_131.fq"), path("${id}_unmapped_131.fq") into mapped_out148618787703226368
		set val(id), path("${id}_mapped_132.fq"), path("${id}_unmapped_132.fq") into mapped_out149744687610068992
		set val(id), path("${id}_mapped_133.fq"), path("${id}_unmapped_133.fq") into mapped_out150870587516911616
		set val(id), path("${id}_mapped_134.fq"), path("${id}_unmapped_134.fq") into mapped_out151996487423754240
		set val(id), path("${id}_mapped_135.fq"), path("${id}_unmapped_135.fq") into mapped_out153122387330596864
		set val(id), path("${id}_mapped_136.fq"), path("${id}_unmapped_136.fq") into mapped_out154248287237439488
		set val(id), path("${id}_mapped_137.fq"), path("${id}_unmapped_137.fq") into mapped_out155374187144282112
		set val(id), path("${id}_mapped_138.fq"), path("${id}_unmapped_138.fq") into mapped_out156500087051124736
		set val(id), path("${id}_mapped_139.fq"), path("${id}_unmapped_139.fq") into mapped_out157625986957967360
		set val(id), path("${id}_mapped_140.fq"), path("${id}_unmapped_140.fq") into mapped_out158751886864809984
		set val(id), path("${id}_mapped_141.fq"), path("${id}_unmapped_141.fq") into mapped_out159877786771652608
		set val(id), path("${id}_mapped_142.fq"), path("${id}_unmapped_142.fq") into mapped_out161003686678495232
		set val(id), path("${id}_mapped_143.fq"), path("${id}_unmapped_143.fq") into mapped_out162129586585337856
		set val(id), path("${id}_mapped_144.fq"), path("${id}_unmapped_144.fq") into mapped_out163255486492180480
		set val(id), path("${id}_mapped_145.fq"), path("${id}_unmapped_145.fq") into mapped_out164381386399023104
		set val(id), path("${id}_mapped_146.fq"), path("${id}_unmapped_146.fq") into mapped_out165507286305865728
		set val(id), path("${id}_mapped_147.fq"), path("${id}_unmapped_147.fq") into mapped_out166633186212708352
		set val(id), path("${id}_mapped_148.fq"), path("${id}_unmapped_148.fq") into mapped_out167759086119550976
		set val(id), path("${id}_mapped_149.fq"), path("${id}_unmapped_149.fq") into mapped_out168884986026393600
		set val(id), path("${id}_mapped_150.fq"), path("${id}_unmapped_150.fq") into mapped_out170010885933236224
		set val(id), path("${id}_mapped_151.fq"), path("${id}_unmapped_151.fq") into mapped_out171136785840078848
		set val(id), path("${id}_mapped_152.fq"), path("${id}_unmapped_152.fq") into mapped_out172262685746921472
		set val(id), path("${id}_mapped_153.fq"), path("${id}_unmapped_153.fq") into mapped_out173388585653764096
		set val(id), path("${id}_mapped_154.fq"), path("${id}_unmapped_154.fq") into mapped_out174514485560606720
		set val(id), path("${id}_mapped_155.fq"), path("${id}_unmapped_155.fq") into mapped_out175640385467449344
		set val(id), path("${id}_mapped_156.fq"), path("${id}_unmapped_156.fq") into mapped_out176766285374291968
		set val(id), path("${id}_mapped_157.fq"), path("${id}_unmapped_157.fq") into mapped_out177892185281134592
		set val(id), path("${id}_mapped_158.fq"), path("${id}_unmapped_158.fq") into mapped_out179018085187977216
		set val(id), path("${id}_mapped_159.fq"), path("${id}_unmapped_159.fq") into mapped_out180143985094819840
		set val(id), path("${id}_mapped_160.fq"), path("${id}_unmapped_160.fq") into mapped_out181269885001662464
		set val(id), path("${id}_mapped_161.fq"), path("${id}_unmapped_161.fq") into mapped_out182395784908505088
		set val(id), path("${id}_mapped_162.fq"), path("${id}_unmapped_162.fq") into mapped_out183521684815347712
		set val(id), path("${id}_mapped_163.fq"), path("${id}_unmapped_163.fq") into mapped_out184647584722190336
		set val(id), path("${id}_mapped_164.fq"), path("${id}_unmapped_164.fq") into mapped_out185773484629032960
		set val(id), path("${id}_mapped_165.fq"), path("${id}_unmapped_165.fq") into mapped_out186899384535875584
		set val(id), path("${id}_mapped_166.fq"), path("${id}_unmapped_166.fq") into mapped_out188025284442718208
		set val(id), path("${id}_mapped_167.fq"), path("${id}_unmapped_167.fq") into mapped_out189151184349560832
		set val(id), path("${id}_mapped_168.fq"), path("${id}_unmapped_168.fq") into mapped_out190277084256403456
		set val(id), path("${id}_mapped_169.fq"), path("${id}_unmapped_169.fq") into mapped_out191402984163246080
		set val(id), path("${id}_mapped_170.fq"), path("${id}_unmapped_170.fq") into mapped_out192528884070088704
		set val(id), path("${id}_mapped_171.fq"), path("${id}_unmapped_171.fq") into mapped_out193654783976931328
		set val(id), path("${id}_mapped_172.fq"), path("${id}_unmapped_172.fq") into mapped_out194780683883773952
		set val(id), path("${id}_mapped_173.fq"), path("${id}_unmapped_173.fq") into mapped_out195906583790616576
		set val(id), path("${id}_mapped_174.fq"), path("${id}_unmapped_174.fq") into mapped_out197032483697459200
		set val(id), path("${id}_mapped_175.fq"), path("${id}_unmapped_175.fq") into mapped_out198158383604301824
		set val(id), path("${id}_mapped_176.fq"), path("${id}_unmapped_176.fq") into mapped_out199284283511144448
		set val(id), path("${id}_mapped_177.fq"), path("${id}_unmapped_177.fq") into mapped_out200410183417987072
		set val(id), path("${id}_mapped_178.fq"), path("${id}_unmapped_178.fq") into mapped_out201536083324829696
		set val(id), path("${id}_mapped_179.fq"), path("${id}_unmapped_179.fq") into mapped_out202661983231672320
		set val(id), path("${id}_mapped_180.fq"), path("${id}_unmapped_180.fq") into mapped_out203787883138514944
		set val(id), path("${id}_mapped_181.fq"), path("${id}_unmapped_181.fq") into mapped_out204913783045357568
		set val(id), path("${id}_mapped_182.fq"), path("${id}_unmapped_182.fq") into mapped_out206039682952200192

    script:
		"""
		bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_1.fq outm=${id}_mapped_1.fq outu=${id}_unmapped_1.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_2.fq outm=${id}_mapped_2.fq outu=${id}_unmapped_2.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_3.fq outm=${id}_mapped_3.fq outu=${id}_unmapped_3.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_4.fq outm=${id}_mapped_4.fq outu=${id}_unmapped_4.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_5.fq outm=${id}_mapped_5.fq outu=${id}_unmapped_5.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_6.fq outm=${id}_mapped_6.fq outu=${id}_unmapped_6.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_7.fq outm=${id}_mapped_7.fq outu=${id}_unmapped_7.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_8.fq outm=${id}_mapped_8.fq outu=${id}_unmapped_8.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_9.fq outm=${id}_mapped_9.fq outu=${id}_unmapped_9.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_10.fq outm=${id}_mapped_10.fq outu=${id}_unmapped_10.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_11.fq outm=${id}_mapped_11.fq outu=${id}_unmapped_11.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_12.fq outm=${id}_mapped_12.fq outu=${id}_unmapped_12.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_13.fq outm=${id}_mapped_13.fq outu=${id}_unmapped_13.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_14.fq outm=${id}_mapped_14.fq outu=${id}_unmapped_14.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_15.fq outm=${id}_mapped_15.fq outu=${id}_unmapped_15.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_16.fq outm=${id}_mapped_16.fq outu=${id}_unmapped_16.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_17.fq outm=${id}_mapped_17.fq outu=${id}_unmapped_17.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_18.fq outm=${id}_mapped_18.fq outu=${id}_unmapped_18.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_19.fq outm=${id}_mapped_19.fq outu=${id}_unmapped_19.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_20.fq outm=${id}_mapped_20.fq outu=${id}_unmapped_20.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_21.fq outm=${id}_mapped_21.fq outu=${id}_unmapped_21.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_22.fq outm=${id}_mapped_22.fq outu=${id}_unmapped_22.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_23.fq outm=${id}_mapped_23.fq outu=${id}_unmapped_23.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_24.fq outm=${id}_mapped_24.fq outu=${id}_unmapped_24.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_25.fq outm=${id}_mapped_25.fq outu=${id}_unmapped_25.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_26.fq outm=${id}_mapped_26.fq outu=${id}_unmapped_26.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_27.fq outm=${id}_mapped_27.fq outu=${id}_unmapped_27.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_28.fq outm=${id}_mapped_28.fq outu=${id}_unmapped_28.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_29.fq outm=${id}_mapped_29.fq outu=${id}_unmapped_29.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_30.fq outm=${id}_mapped_30.fq outu=${id}_unmapped_30.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_31.fq outm=${id}_mapped_31.fq outu=${id}_unmapped_31.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_32.fq outm=${id}_mapped_32.fq outu=${id}_unmapped_32.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_33.fq outm=${id}_mapped_33.fq outu=${id}_unmapped_33.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_34.fq outm=${id}_mapped_34.fq outu=${id}_unmapped_34.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_35.fq outm=${id}_mapped_35.fq outu=${id}_unmapped_35.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_36.fq outm=${id}_mapped_36.fq outu=${id}_unmapped_36.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_37.fq outm=${id}_mapped_37.fq outu=${id}_unmapped_37.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_38.fq outm=${id}_mapped_38.fq outu=${id}_unmapped_38.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_39.fq outm=${id}_mapped_39.fq outu=${id}_unmapped_39.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_40.fq outm=${id}_mapped_40.fq outu=${id}_unmapped_40.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_41.fq outm=${id}_mapped_41.fq outu=${id}_unmapped_41.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_42.fq outm=${id}_mapped_42.fq outu=${id}_unmapped_42.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_43.fq outm=${id}_mapped_43.fq outu=${id}_unmapped_43.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_44.fq outm=${id}_mapped_44.fq outu=${id}_unmapped_44.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_45.fq outm=${id}_mapped_45.fq outu=${id}_unmapped_45.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_46.fq outm=${id}_mapped_46.fq outu=${id}_unmapped_46.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_47.fq outm=${id}_mapped_47.fq outu=${id}_unmapped_47.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_48.fq outm=${id}_mapped_48.fq outu=${id}_unmapped_48.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_49.fq outm=${id}_mapped_49.fq outu=${id}_unmapped_49.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_50.fq outm=${id}_mapped_50.fq outu=${id}_unmapped_50.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_51.fq outm=${id}_mapped_51.fq outu=${id}_unmapped_51.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_52.fq outm=${id}_mapped_52.fq outu=${id}_unmapped_52.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_53.fq outm=${id}_mapped_53.fq outu=${id}_unmapped_53.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_54.fq outm=${id}_mapped_54.fq outu=${id}_unmapped_54.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_55.fq outm=${id}_mapped_55.fq outu=${id}_unmapped_55.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_56.fq outm=${id}_mapped_56.fq outu=${id}_unmapped_56.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_57.fq outm=${id}_mapped_57.fq outu=${id}_unmapped_57.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_58.fq outm=${id}_mapped_58.fq outu=${id}_unmapped_58.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_59.fq outm=${id}_mapped_59.fq outu=${id}_unmapped_59.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_60.fq outm=${id}_mapped_60.fq outu=${id}_unmapped_60.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_61.fq outm=${id}_mapped_61.fq outu=${id}_unmapped_61.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_62.fq outm=${id}_mapped_62.fq outu=${id}_unmapped_62.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_63.fq outm=${id}_mapped_63.fq outu=${id}_unmapped_63.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_64.fq outm=${id}_mapped_64.fq outu=${id}_unmapped_64.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_65.fq outm=${id}_mapped_65.fq outu=${id}_unmapped_65.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_66.fq outm=${id}_mapped_66.fq outu=${id}_unmapped_66.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_67.fq outm=${id}_mapped_67.fq outu=${id}_unmapped_67.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_68.fq outm=${id}_mapped_68.fq outu=${id}_unmapped_68.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_69.fq outm=${id}_mapped_69.fq outu=${id}_unmapped_69.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_70.fq outm=${id}_mapped_70.fq outu=${id}_unmapped_70.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_71.fq outm=${id}_mapped_71.fq outu=${id}_unmapped_71.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_72.fq outm=${id}_mapped_72.fq outu=${id}_unmapped_72.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_73.fq outm=${id}_mapped_73.fq outu=${id}_unmapped_73.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_74.fq outm=${id}_mapped_74.fq outu=${id}_unmapped_74.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_75.fq outm=${id}_mapped_75.fq outu=${id}_unmapped_75.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_76.fq outm=${id}_mapped_76.fq outu=${id}_unmapped_76.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_77.fq outm=${id}_mapped_77.fq outu=${id}_unmapped_77.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_78.fq outm=${id}_mapped_78.fq outu=${id}_unmapped_78.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_79.fq outm=${id}_mapped_79.fq outu=${id}_unmapped_79.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_80.fq outm=${id}_mapped_80.fq outu=${id}_unmapped_80.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_81.fq outm=${id}_mapped_81.fq outu=${id}_unmapped_81.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_82.fq outm=${id}_mapped_82.fq outu=${id}_unmapped_82.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_83.fq outm=${id}_mapped_83.fq outu=${id}_unmapped_83.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_84.fq outm=${id}_mapped_84.fq outu=${id}_unmapped_84.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_85.fq outm=${id}_mapped_85.fq outu=${id}_unmapped_85.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_86.fq outm=${id}_mapped_86.fq outu=${id}_unmapped_86.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_87.fq outm=${id}_mapped_87.fq outu=${id}_unmapped_87.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_88.fq outm=${id}_mapped_88.fq outu=${id}_unmapped_88.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_89.fq outm=${id}_mapped_89.fq outu=${id}_unmapped_89.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_90.fq outm=${id}_mapped_90.fq outu=${id}_unmapped_90.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_91.fq outm=${id}_mapped_91.fq outu=${id}_unmapped_91.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_92.fq outm=${id}_mapped_92.fq outu=${id}_unmapped_92.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_93.fq outm=${id}_mapped_93.fq outu=${id}_unmapped_93.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_94.fq outm=${id}_mapped_94.fq outu=${id}_unmapped_94.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_95.fq outm=${id}_mapped_95.fq outu=${id}_unmapped_95.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_96.fq outm=${id}_mapped_96.fq outu=${id}_unmapped_96.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_97.fq outm=${id}_mapped_97.fq outu=${id}_unmapped_97.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_98.fq outm=${id}_mapped_98.fq outu=${id}_unmapped_98.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_99.fq outm=${id}_mapped_99.fq outu=${id}_unmapped_99.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_100.fq outm=${id}_mapped_100.fq outu=${id}_unmapped_100.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_101.fq outm=${id}_mapped_101.fq outu=${id}_unmapped_101.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_102.fq outm=${id}_mapped_102.fq outu=${id}_unmapped_102.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_103.fq outm=${id}_mapped_103.fq outu=${id}_unmapped_103.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_104.fq outm=${id}_mapped_104.fq outu=${id}_unmapped_104.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_105.fq outm=${id}_mapped_105.fq outu=${id}_unmapped_105.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_106.fq outm=${id}_mapped_106.fq outu=${id}_unmapped_106.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_107.fq outm=${id}_mapped_107.fq outu=${id}_unmapped_107.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_108.fq outm=${id}_mapped_108.fq outu=${id}_unmapped_108.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_109.fq outm=${id}_mapped_109.fq outu=${id}_unmapped_109.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_110.fq outm=${id}_mapped_110.fq outu=${id}_unmapped_110.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_111.fq outm=${id}_mapped_111.fq outu=${id}_unmapped_111.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_112.fq outm=${id}_mapped_112.fq outu=${id}_unmapped_112.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_113.fq outm=${id}_mapped_113.fq outu=${id}_unmapped_113.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_114.fq outm=${id}_mapped_114.fq outu=${id}_unmapped_114.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_115.fq outm=${id}_mapped_115.fq outu=${id}_unmapped_115.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_116.fq outm=${id}_mapped_116.fq outu=${id}_unmapped_116.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_117.fq outm=${id}_mapped_117.fq outu=${id}_unmapped_117.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_118.fq outm=${id}_mapped_118.fq outu=${id}_unmapped_118.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_119.fq outm=${id}_mapped_119.fq outu=${id}_unmapped_119.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_120.fq outm=${id}_mapped_120.fq outu=${id}_unmapped_120.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_121.fq outm=${id}_mapped_121.fq outu=${id}_unmapped_121.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_122.fq outm=${id}_mapped_122.fq outu=${id}_unmapped_122.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_123.fq outm=${id}_mapped_123.fq outu=${id}_unmapped_123.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_124.fq outm=${id}_mapped_124.fq outu=${id}_unmapped_124.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_125.fq outm=${id}_mapped_125.fq outu=${id}_unmapped_125.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_126.fq outm=${id}_mapped_126.fq outu=${id}_unmapped_126.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_127.fq outm=${id}_mapped_127.fq outu=${id}_unmapped_127.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_128.fq outm=${id}_mapped_128.fq outu=${id}_unmapped_128.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_129.fq outm=${id}_mapped_129.fq outu=${id}_unmapped_129.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_130.fq outm=${id}_mapped_130.fq outu=${id}_unmapped_130.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_131.fq outm=${id}_mapped_131.fq outu=${id}_unmapped_131.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_132.fq outm=${id}_mapped_132.fq outu=${id}_unmapped_132.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_133.fq outm=${id}_mapped_133.fq outu=${id}_unmapped_133.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_134.fq outm=${id}_mapped_134.fq outu=${id}_unmapped_134.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_135.fq outm=${id}_mapped_135.fq outu=${id}_unmapped_135.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_136.fq outm=${id}_mapped_136.fq outu=${id}_unmapped_136.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_137.fq outm=${id}_mapped_137.fq outu=${id}_unmapped_137.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_138.fq outm=${id}_mapped_138.fq outu=${id}_unmapped_138.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_139.fq outm=${id}_mapped_139.fq outu=${id}_unmapped_139.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_140.fq outm=${id}_mapped_140.fq outu=${id}_unmapped_140.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_141.fq outm=${id}_mapped_141.fq outu=${id}_unmapped_141.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_142.fq outm=${id}_mapped_142.fq outu=${id}_unmapped_142.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_143.fq outm=${id}_mapped_143.fq outu=${id}_unmapped_143.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_144.fq outm=${id}_mapped_144.fq outu=${id}_unmapped_144.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_145.fq outm=${id}_mapped_145.fq outu=${id}_unmapped_145.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_146.fq outm=${id}_mapped_146.fq outu=${id}_unmapped_146.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_147.fq outm=${id}_mapped_147.fq outu=${id}_unmapped_147.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_148.fq outm=${id}_mapped_148.fq outu=${id}_unmapped_148.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_149.fq outm=${id}_mapped_149.fq outu=${id}_unmapped_149.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_150.fq outm=${id}_mapped_150.fq outu=${id}_unmapped_150.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_151.fq outm=${id}_mapped_151.fq outu=${id}_unmapped_151.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_152.fq outm=${id}_mapped_152.fq outu=${id}_unmapped_152.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_153.fq outm=${id}_mapped_153.fq outu=${id}_unmapped_153.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_154.fq outm=${id}_mapped_154.fq outu=${id}_unmapped_154.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_155.fq outm=${id}_mapped_155.fq outu=${id}_unmapped_155.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_156.fq outm=${id}_mapped_156.fq outu=${id}_unmapped_156.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_157.fq outm=${id}_mapped_157.fq outu=${id}_unmapped_157.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_158.fq outm=${id}_mapped_158.fq outu=${id}_unmapped_158.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_159.fq outm=${id}_mapped_159.fq outu=${id}_unmapped_159.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_160.fq outm=${id}_mapped_160.fq outu=${id}_unmapped_160.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_161.fq outm=${id}_mapped_161.fq outu=${id}_unmapped_161.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_162.fq outm=${id}_mapped_162.fq outu=${id}_unmapped_162.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_163.fq outm=${id}_mapped_163.fq outu=${id}_unmapped_163.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_164.fq outm=${id}_mapped_164.fq outu=${id}_unmapped_164.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_165.fq outm=${id}_mapped_165.fq outu=${id}_unmapped_165.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_166.fq outm=${id}_mapped_166.fq outu=${id}_unmapped_166.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_167.fq outm=${id}_mapped_167.fq outu=${id}_unmapped_167.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_168.fq outm=${id}_mapped_168.fq outu=${id}_unmapped_168.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_169.fq outm=${id}_mapped_169.fq outu=${id}_unmapped_169.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_170.fq outm=${id}_mapped_170.fq outu=${id}_unmapped_170.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_171.fq outm=${id}_mapped_171.fq outu=${id}_unmapped_171.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_172.fq outm=${id}_mapped_172.fq outu=${id}_unmapped_172.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_173.fq outm=${id}_mapped_173.fq outu=${id}_unmapped_173.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_174.fq outm=${id}_mapped_174.fq outu=${id}_unmapped_174.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_175.fq outm=${id}_mapped_175.fq outu=${id}_unmapped_175.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_176.fq outm=${id}_mapped_176.fq outu=${id}_unmapped_176.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_177.fq outm=${id}_mapped_177.fq outu=${id}_unmapped_177.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_178.fq outm=${id}_mapped_178.fq outu=${id}_unmapped_178.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_179.fq outm=${id}_mapped_179.fq outu=${id}_unmapped_179.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_180.fq outm=${id}_mapped_180.fq outu=${id}_unmapped_180.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_181.fq outm=${id}_mapped_181.fq outu=${id}_unmapped_181.fq ref=$ref_dir/XM_002808697.2.fasta
		bbduk.sh  -Xmx1g -in=${id}_trimmed_182.fq outm=${id}_mapped_182.fq outu=${id}_unmapped_182.fq ref=$ref_dir/XM_002808697.2.fasta
		"""


}


process assemble32 {
	errorStrategy 'ignore'

    publishDir "$params.output.folder32/contigs/${sample}", mode : "copy"
    input:
		set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1125899906842624
		set val(sample), path(mapped_read1), path(unmapped_read1) from mapped_out2251799813685248
		set val(sample), path(mapped_read2), path(unmapped_read2) from mapped_out3377699720527872
		set val(sample), path(mapped_read3), path(unmapped_read3) from mapped_out4503599627370496
		set val(sample), path(mapped_read4), path(unmapped_read4) from mapped_out5629499534213120
		set val(sample), path(mapped_read5), path(unmapped_read5) from mapped_out6755399441055744
		set val(sample), path(mapped_read6), path(unmapped_read6) from mapped_out7881299347898368
		set val(sample), path(mapped_read7), path(unmapped_read7) from mapped_out9007199254740992
		set val(sample), path(mapped_read8), path(unmapped_read8) from mapped_out10133099161583616
		set val(sample), path(mapped_read9), path(unmapped_read9) from mapped_out11258999068426240
		set val(sample), path(mapped_read10), path(unmapped_read10) from mapped_out12384898975268864
		set val(sample), path(mapped_read11), path(unmapped_read11) from mapped_out13510798882111488
		set val(sample), path(mapped_read12), path(unmapped_read12) from mapped_out14636698788954112
		set val(sample), path(mapped_read13), path(unmapped_read13) from mapped_out15762598695796736
		set val(sample), path(mapped_read14), path(unmapped_read14) from mapped_out16888498602639360
		set val(sample), path(mapped_read15), path(unmapped_read15) from mapped_out18014398509481984
		set val(sample), path(mapped_read16), path(unmapped_read16) from mapped_out19140298416324608
		set val(sample), path(mapped_read17), path(unmapped_read17) from mapped_out20266198323167232
		set val(sample), path(mapped_read18), path(unmapped_read18) from mapped_out21392098230009856
		set val(sample), path(mapped_read19), path(unmapped_read19) from mapped_out22517998136852480
		set val(sample), path(mapped_read20), path(unmapped_read20) from mapped_out23643898043695104
		set val(sample), path(mapped_read21), path(unmapped_read21) from mapped_out24769797950537728
		set val(sample), path(mapped_read22), path(unmapped_read22) from mapped_out25895697857380352
		set val(sample), path(mapped_read23), path(unmapped_read23) from mapped_out27021597764222976
		set val(sample), path(mapped_read24), path(unmapped_read24) from mapped_out28147497671065600
		set val(sample), path(mapped_read25), path(unmapped_read25) from mapped_out29273397577908224
		set val(sample), path(mapped_read26), path(unmapped_read26) from mapped_out30399297484750848
		set val(sample), path(mapped_read27), path(unmapped_read27) from mapped_out31525197391593472
		set val(sample), path(mapped_read28), path(unmapped_read28) from mapped_out32651097298436096
		set val(sample), path(mapped_read29), path(unmapped_read29) from mapped_out33776997205278720
		set val(sample), path(mapped_read30), path(unmapped_read30) from mapped_out34902897112121344
		set val(sample), path(mapped_read31), path(unmapped_read31) from mapped_out36028797018963968
		set val(sample), path(mapped_read32), path(unmapped_read32) from mapped_out37154696925806592
		set val(sample), path(mapped_read33), path(unmapped_read33) from mapped_out38280596832649216
		set val(sample), path(mapped_read34), path(unmapped_read34) from mapped_out39406496739491840
		set val(sample), path(mapped_read35), path(unmapped_read35) from mapped_out40532396646334464
		set val(sample), path(mapped_read36), path(unmapped_read36) from mapped_out41658296553177088
		set val(sample), path(mapped_read37), path(unmapped_read37) from mapped_out42784196460019712
		set val(sample), path(mapped_read38), path(unmapped_read38) from mapped_out43910096366862336
		set val(sample), path(mapped_read39), path(unmapped_read39) from mapped_out45035996273704960
		set val(sample), path(mapped_read40), path(unmapped_read40) from mapped_out46161896180547584
		set val(sample), path(mapped_read41), path(unmapped_read41) from mapped_out47287796087390208
		set val(sample), path(mapped_read42), path(unmapped_read42) from mapped_out48413695994232832
		set val(sample), path(mapped_read43), path(unmapped_read43) from mapped_out49539595901075456
		set val(sample), path(mapped_read44), path(unmapped_read44) from mapped_out50665495807918080
		set val(sample), path(mapped_read45), path(unmapped_read45) from mapped_out51791395714760704
		set val(sample), path(mapped_read46), path(unmapped_read46) from mapped_out52917295621603328
		set val(sample), path(mapped_read47), path(unmapped_read47) from mapped_out54043195528445952
		set val(sample), path(mapped_read48), path(unmapped_read48) from mapped_out55169095435288576
		set val(sample), path(mapped_read49), path(unmapped_read49) from mapped_out56294995342131200
		set val(sample), path(mapped_read50), path(unmapped_read50) from mapped_out57420895248973824
		set val(sample), path(mapped_read51), path(unmapped_read51) from mapped_out58546795155816448
		set val(sample), path(mapped_read52), path(unmapped_read52) from mapped_out59672695062659072
		set val(sample), path(mapped_read53), path(unmapped_read53) from mapped_out60798594969501696
		set val(sample), path(mapped_read54), path(unmapped_read54) from mapped_out61924494876344320
		set val(sample), path(mapped_read55), path(unmapped_read55) from mapped_out63050394783186944
		set val(sample), path(mapped_read56), path(unmapped_read56) from mapped_out64176294690029568
		set val(sample), path(mapped_read57), path(unmapped_read57) from mapped_out65302194596872192
		set val(sample), path(mapped_read58), path(unmapped_read58) from mapped_out66428094503714816
		set val(sample), path(mapped_read59), path(unmapped_read59) from mapped_out67553994410557440
		set val(sample), path(mapped_read60), path(unmapped_read60) from mapped_out68679894317400064
		set val(sample), path(mapped_read61), path(unmapped_read61) from mapped_out69805794224242688
		set val(sample), path(mapped_read62), path(unmapped_read62) from mapped_out70931694131085312
		set val(sample), path(mapped_read63), path(unmapped_read63) from mapped_out72057594037927936
		set val(sample), path(mapped_read64), path(unmapped_read64) from mapped_out73183493944770560
		set val(sample), path(mapped_read65), path(unmapped_read65) from mapped_out74309393851613184
		set val(sample), path(mapped_read66), path(unmapped_read66) from mapped_out75435293758455808
		set val(sample), path(mapped_read67), path(unmapped_read67) from mapped_out76561193665298432
		set val(sample), path(mapped_read68), path(unmapped_read68) from mapped_out77687093572141056
		set val(sample), path(mapped_read69), path(unmapped_read69) from mapped_out78812993478983680
		set val(sample), path(mapped_read70), path(unmapped_read70) from mapped_out79938893385826304
		set val(sample), path(mapped_read71), path(unmapped_read71) from mapped_out81064793292668928
		set val(sample), path(mapped_read72), path(unmapped_read72) from mapped_out82190693199511552
		set val(sample), path(mapped_read73), path(unmapped_read73) from mapped_out83316593106354176
		set val(sample), path(mapped_read74), path(unmapped_read74) from mapped_out84442493013196800
		set val(sample), path(mapped_read75), path(unmapped_read75) from mapped_out85568392920039424
		set val(sample), path(mapped_read76), path(unmapped_read76) from mapped_out86694292826882048
		set val(sample), path(mapped_read77), path(unmapped_read77) from mapped_out87820192733724672
		set val(sample), path(mapped_read78), path(unmapped_read78) from mapped_out88946092640567296
		set val(sample), path(mapped_read79), path(unmapped_read79) from mapped_out90071992547409920
		set val(sample), path(mapped_read80), path(unmapped_read80) from mapped_out91197892454252544
		set val(sample), path(mapped_read81), path(unmapped_read81) from mapped_out92323792361095168
		set val(sample), path(mapped_read82), path(unmapped_read82) from mapped_out93449692267937792
		set val(sample), path(mapped_read83), path(unmapped_read83) from mapped_out94575592174780416
		set val(sample), path(mapped_read84), path(unmapped_read84) from mapped_out95701492081623040
		set val(sample), path(mapped_read85), path(unmapped_read85) from mapped_out96827391988465664
		set val(sample), path(mapped_read86), path(unmapped_read86) from mapped_out97953291895308288
		set val(sample), path(mapped_read87), path(unmapped_read87) from mapped_out99079191802150912
		set val(sample), path(mapped_read88), path(unmapped_read88) from mapped_out100205091708993536
		set val(sample), path(mapped_read89), path(unmapped_read89) from mapped_out101330991615836160
		set val(sample), path(mapped_read90), path(unmapped_read90) from mapped_out102456891522678784
		set val(sample), path(mapped_read91), path(unmapped_read91) from mapped_out103582791429521408
		set val(sample), path(mapped_read92), path(unmapped_read92) from mapped_out104708691336364032
		set val(sample), path(mapped_read93), path(unmapped_read93) from mapped_out105834591243206656
		set val(sample), path(mapped_read94), path(unmapped_read94) from mapped_out106960491150049280
		set val(sample), path(mapped_read95), path(unmapped_read95) from mapped_out108086391056891904
		set val(sample), path(mapped_read96), path(unmapped_read96) from mapped_out109212290963734528
		set val(sample), path(mapped_read97), path(unmapped_read97) from mapped_out110338190870577152
		set val(sample), path(mapped_read98), path(unmapped_read98) from mapped_out111464090777419776
		set val(sample), path(mapped_read99), path(unmapped_read99) from mapped_out112589990684262400
		set val(sample), path(mapped_read100), path(unmapped_read100) from mapped_out113715890591105024
		set val(sample), path(mapped_read101), path(unmapped_read101) from mapped_out114841790497947648
		set val(sample), path(mapped_read102), path(unmapped_read102) from mapped_out115967690404790272
		set val(sample), path(mapped_read103), path(unmapped_read103) from mapped_out117093590311632896
		set val(sample), path(mapped_read104), path(unmapped_read104) from mapped_out118219490218475520
		set val(sample), path(mapped_read105), path(unmapped_read105) from mapped_out119345390125318144
		set val(sample), path(mapped_read106), path(unmapped_read106) from mapped_out120471290032160768
		set val(sample), path(mapped_read107), path(unmapped_read107) from mapped_out121597189939003392
		set val(sample), path(mapped_read108), path(unmapped_read108) from mapped_out122723089845846016
		set val(sample), path(mapped_read109), path(unmapped_read109) from mapped_out123848989752688640
		set val(sample), path(mapped_read110), path(unmapped_read110) from mapped_out124974889659531264
		set val(sample), path(mapped_read111), path(unmapped_read111) from mapped_out126100789566373888
		set val(sample), path(mapped_read112), path(unmapped_read112) from mapped_out127226689473216512
		set val(sample), path(mapped_read113), path(unmapped_read113) from mapped_out128352589380059136
		set val(sample), path(mapped_read114), path(unmapped_read114) from mapped_out129478489286901760
		set val(sample), path(mapped_read115), path(unmapped_read115) from mapped_out130604389193744384
		set val(sample), path(mapped_read116), path(unmapped_read116) from mapped_out131730289100587008
		set val(sample), path(mapped_read117), path(unmapped_read117) from mapped_out132856189007429632
		set val(sample), path(mapped_read118), path(unmapped_read118) from mapped_out133982088914272256
		set val(sample), path(mapped_read119), path(unmapped_read119) from mapped_out135107988821114880
		set val(sample), path(mapped_read120), path(unmapped_read120) from mapped_out136233888727957504
		set val(sample), path(mapped_read121), path(unmapped_read121) from mapped_out137359788634800128
		set val(sample), path(mapped_read122), path(unmapped_read122) from mapped_out138485688541642752
		set val(sample), path(mapped_read123), path(unmapped_read123) from mapped_out139611588448485376
		set val(sample), path(mapped_read124), path(unmapped_read124) from mapped_out140737488355328000
		set val(sample), path(mapped_read125), path(unmapped_read125) from mapped_out141863388262170624
		set val(sample), path(mapped_read126), path(unmapped_read126) from mapped_out142989288169013248
		set val(sample), path(mapped_read127), path(unmapped_read127) from mapped_out144115188075855872
		set val(sample), path(mapped_read128), path(unmapped_read128) from mapped_out145241087982698496
		set val(sample), path(mapped_read129), path(unmapped_read129) from mapped_out146366987889541120
		set val(sample), path(mapped_read130), path(unmapped_read130) from mapped_out147492887796383744
		set val(sample), path(mapped_read131), path(unmapped_read131) from mapped_out148618787703226368
		set val(sample), path(mapped_read132), path(unmapped_read132) from mapped_out149744687610068992
		set val(sample), path(mapped_read133), path(unmapped_read133) from mapped_out150870587516911616
		set val(sample), path(mapped_read134), path(unmapped_read134) from mapped_out151996487423754240
		set val(sample), path(mapped_read135), path(unmapped_read135) from mapped_out153122387330596864
		set val(sample), path(mapped_read136), path(unmapped_read136) from mapped_out154248287237439488
		set val(sample), path(mapped_read137), path(unmapped_read137) from mapped_out155374187144282112
		set val(sample), path(mapped_read138), path(unmapped_read138) from mapped_out156500087051124736
		set val(sample), path(mapped_read139), path(unmapped_read139) from mapped_out157625986957967360
		set val(sample), path(mapped_read140), path(unmapped_read140) from mapped_out158751886864809984
		set val(sample), path(mapped_read141), path(unmapped_read141) from mapped_out159877786771652608
		set val(sample), path(mapped_read142), path(unmapped_read142) from mapped_out161003686678495232
		set val(sample), path(mapped_read143), path(unmapped_read143) from mapped_out162129586585337856
		set val(sample), path(mapped_read144), path(unmapped_read144) from mapped_out163255486492180480
		set val(sample), path(mapped_read145), path(unmapped_read145) from mapped_out164381386399023104
		set val(sample), path(mapped_read146), path(unmapped_read146) from mapped_out165507286305865728
		set val(sample), path(mapped_read147), path(unmapped_read147) from mapped_out166633186212708352
		set val(sample), path(mapped_read148), path(unmapped_read148) from mapped_out167759086119550976
		set val(sample), path(mapped_read149), path(unmapped_read149) from mapped_out168884986026393600
		set val(sample), path(mapped_read150), path(unmapped_read150) from mapped_out170010885933236224
		set val(sample), path(mapped_read151), path(unmapped_read151) from mapped_out171136785840078848
		set val(sample), path(mapped_read152), path(unmapped_read152) from mapped_out172262685746921472
		set val(sample), path(mapped_read153), path(unmapped_read153) from mapped_out173388585653764096
		set val(sample), path(mapped_read154), path(unmapped_read154) from mapped_out174514485560606720
		set val(sample), path(mapped_read155), path(unmapped_read155) from mapped_out175640385467449344
		set val(sample), path(mapped_read156), path(unmapped_read156) from mapped_out176766285374291968
		set val(sample), path(mapped_read157), path(unmapped_read157) from mapped_out177892185281134592
		set val(sample), path(mapped_read158), path(unmapped_read158) from mapped_out179018085187977216
		set val(sample), path(mapped_read159), path(unmapped_read159) from mapped_out180143985094819840
		set val(sample), path(mapped_read160), path(unmapped_read160) from mapped_out181269885001662464
		set val(sample), path(mapped_read161), path(unmapped_read161) from mapped_out182395784908505088
		set val(sample), path(mapped_read162), path(unmapped_read162) from mapped_out183521684815347712
		set val(sample), path(mapped_read163), path(unmapped_read163) from mapped_out184647584722190336
		set val(sample), path(mapped_read164), path(unmapped_read164) from mapped_out185773484629032960
		set val(sample), path(mapped_read165), path(unmapped_read165) from mapped_out186899384535875584
		set val(sample), path(mapped_read166), path(unmapped_read166) from mapped_out188025284442718208
		set val(sample), path(mapped_read167), path(unmapped_read167) from mapped_out189151184349560832
		set val(sample), path(mapped_read168), path(unmapped_read168) from mapped_out190277084256403456
		set val(sample), path(mapped_read169), path(unmapped_read169) from mapped_out191402984163246080
		set val(sample), path(mapped_read170), path(unmapped_read170) from mapped_out192528884070088704
		set val(sample), path(mapped_read171), path(unmapped_read171) from mapped_out193654783976931328
		set val(sample), path(mapped_read172), path(unmapped_read172) from mapped_out194780683883773952
		set val(sample), path(mapped_read173), path(unmapped_read173) from mapped_out195906583790616576
		set val(sample), path(mapped_read174), path(unmapped_read174) from mapped_out197032483697459200
		set val(sample), path(mapped_read175), path(unmapped_read175) from mapped_out198158383604301824
		set val(sample), path(mapped_read176), path(unmapped_read176) from mapped_out199284283511144448
		set val(sample), path(mapped_read177), path(unmapped_read177) from mapped_out200410183417987072
		set val(sample), path(mapped_read178), path(unmapped_read178) from mapped_out201536083324829696
		set val(sample), path(mapped_read179), path(unmapped_read179) from mapped_out202661983231672320
		set val(sample), path(mapped_read180), path(unmapped_read180) from mapped_out203787883138514944
		set val(sample), path(mapped_read181), path(unmapped_read181) from mapped_out204913783045357568
		set val(sample), path(mapped_read182), path(unmapped_read182) from mapped_out206039682952200192


    output:
        tuple val(sample), path("scaffolds.fasta") into assemblyout1125899906842624
    script:
		"""
		$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py  -s ${sample}_mapped_0.fq  -s ${sample}_mapped_1.fq  -s ${sample}_mapped_2.fq  -s ${sample}_mapped_3.fq  -s ${sample}_mapped_4.fq  -s ${sample}_mapped_5.fq  -s ${sample}_mapped_6.fq  -s ${sample}_mapped_7.fq  -s ${sample}_mapped_8.fq  -s ${sample}_mapped_9.fq  -s ${sample}_mapped_10.fq  -s ${sample}_mapped_11.fq  -s ${sample}_mapped_12.fq  -s ${sample}_mapped_13.fq  -s ${sample}_mapped_14.fq  -s ${sample}_mapped_15.fq  -s ${sample}_mapped_16.fq  -s ${sample}_mapped_17.fq  -s ${sample}_mapped_18.fq  -s ${sample}_mapped_19.fq  -s ${sample}_mapped_20.fq  -s ${sample}_mapped_21.fq  -s ${sample}_mapped_22.fq  -s ${sample}_mapped_23.fq  -s ${sample}_mapped_24.fq  -s ${sample}_mapped_25.fq  -s ${sample}_mapped_26.fq  -s ${sample}_mapped_27.fq  -s ${sample}_mapped_28.fq  -s ${sample}_mapped_29.fq  -s ${sample}_mapped_30.fq  -s ${sample}_mapped_31.fq  -s ${sample}_mapped_32.fq  -s ${sample}_mapped_33.fq  -s ${sample}_mapped_34.fq  -s ${sample}_mapped_35.fq  -s ${sample}_mapped_36.fq  -s ${sample}_mapped_37.fq  -s ${sample}_mapped_38.fq  -s ${sample}_mapped_39.fq  -s ${sample}_mapped_40.fq  -s ${sample}_mapped_41.fq  -s ${sample}_mapped_42.fq  -s ${sample}_mapped_43.fq  -s ${sample}_mapped_44.fq  -s ${sample}_mapped_45.fq  -s ${sample}_mapped_46.fq  -s ${sample}_mapped_47.fq  -s ${sample}_mapped_48.fq  -s ${sample}_mapped_49.fq  -s ${sample}_mapped_50.fq  -s ${sample}_mapped_51.fq  -s ${sample}_mapped_52.fq  -s ${sample}_mapped_53.fq  -s ${sample}_mapped_54.fq  -s ${sample}_mapped_55.fq  -s ${sample}_mapped_56.fq  -s ${sample}_mapped_57.fq  -s ${sample}_mapped_58.fq  -s ${sample}_mapped_59.fq  -s ${sample}_mapped_60.fq  -s ${sample}_mapped_61.fq  -s ${sample}_mapped_62.fq  -s ${sample}_mapped_63.fq  -s ${sample}_mapped_64.fq  -s ${sample}_mapped_65.fq  -s ${sample}_mapped_66.fq  -s ${sample}_mapped_67.fq  -s ${sample}_mapped_68.fq  -s ${sample}_mapped_69.fq  -s ${sample}_mapped_70.fq  -s ${sample}_mapped_71.fq  -s ${sample}_mapped_72.fq  -s ${sample}_mapped_73.fq  -s ${sample}_mapped_74.fq  -s ${sample}_mapped_75.fq  -s ${sample}_mapped_76.fq  -s ${sample}_mapped_77.fq  -s ${sample}_mapped_78.fq  -s ${sample}_mapped_79.fq  -s ${sample}_mapped_80.fq  -s ${sample}_mapped_81.fq  -s ${sample}_mapped_82.fq  -s ${sample}_mapped_83.fq  -s ${sample}_mapped_84.fq  -s ${sample}_mapped_85.fq  -s ${sample}_mapped_86.fq  -s ${sample}_mapped_87.fq  -s ${sample}_mapped_88.fq  -s ${sample}_mapped_89.fq  -s ${sample}_mapped_90.fq  -s ${sample}_mapped_91.fq  -s ${sample}_mapped_92.fq  -s ${sample}_mapped_93.fq  -s ${sample}_mapped_94.fq  -s ${sample}_mapped_95.fq  -s ${sample}_mapped_96.fq  -s ${sample}_mapped_97.fq  -s ${sample}_mapped_98.fq  -s ${sample}_mapped_99.fq  -s ${sample}_mapped_100.fq  -s ${sample}_mapped_101.fq  -s ${sample}_mapped_102.fq  -s ${sample}_mapped_103.fq  -s ${sample}_mapped_104.fq  -s ${sample}_mapped_105.fq  -s ${sample}_mapped_106.fq  -s ${sample}_mapped_107.fq  -s ${sample}_mapped_108.fq  -s ${sample}_mapped_109.fq  -s ${sample}_mapped_110.fq  -s ${sample}_mapped_111.fq  -s ${sample}_mapped_112.fq  -s ${sample}_mapped_113.fq  -s ${sample}_mapped_114.fq  -s ${sample}_mapped_115.fq  -s ${sample}_mapped_116.fq  -s ${sample}_mapped_117.fq  -s ${sample}_mapped_118.fq  -s ${sample}_mapped_119.fq  -s ${sample}_mapped_120.fq  -s ${sample}_mapped_121.fq  -s ${sample}_mapped_122.fq  -s ${sample}_mapped_123.fq  -s ${sample}_mapped_124.fq  -s ${sample}_mapped_125.fq  -s ${sample}_mapped_126.fq  -s ${sample}_mapped_127.fq  -s ${sample}_mapped_128.fq  -s ${sample}_mapped_129.fq  -s ${sample}_mapped_130.fq  -s ${sample}_mapped_131.fq  -s ${sample}_mapped_132.fq  -s ${sample}_mapped_133.fq  -s ${sample}_mapped_134.fq  -s ${sample}_mapped_135.fq  -s ${sample}_mapped_136.fq  -s ${sample}_mapped_137.fq  -s ${sample}_mapped_138.fq  -s ${sample}_mapped_139.fq  -s ${sample}_mapped_140.fq  -s ${sample}_mapped_141.fq  -s ${sample}_mapped_142.fq  -s ${sample}_mapped_143.fq  -s ${sample}_mapped_144.fq  -s ${sample}_mapped_145.fq  -s ${sample}_mapped_146.fq  -s ${sample}_mapped_147.fq  -s ${sample}_mapped_148.fq  -s ${sample}_mapped_149.fq  -s ${sample}_mapped_150.fq  -s ${sample}_mapped_151.fq  -s ${sample}_mapped_152.fq  -s ${sample}_mapped_153.fq  -s ${sample}_mapped_154.fq  -s ${sample}_mapped_155.fq  -s ${sample}_mapped_156.fq  -s ${sample}_mapped_157.fq  -s ${sample}_mapped_158.fq  -s ${sample}_mapped_159.fq  -s ${sample}_mapped_160.fq  -s ${sample}_mapped_161.fq  -s ${sample}_mapped_162.fq  -s ${sample}_mapped_163.fq  -s ${sample}_mapped_164.fq  -s ${sample}_mapped_165.fq  -s ${sample}_mapped_166.fq  -s ${sample}_mapped_167.fq  -s ${sample}_mapped_168.fq  -s ${sample}_mapped_169.fq  -s ${sample}_mapped_170.fq  -s ${sample}_mapped_171.fq  -s ${sample}_mapped_172.fq  -s ${sample}_mapped_173.fq  -s ${sample}_mapped_174.fq  -s ${sample}_mapped_175.fq  -s ${sample}_mapped_176.fq  -s ${sample}_mapped_177.fq  -s ${sample}_mapped_178.fq  -s ${sample}_mapped_179.fq  -s ${sample}_mapped_180.fq  -s ${sample}_mapped_181.fq  -s ${sample}_mapped_182.fq -o . 
		"""
       
        
}
        //skesa --reads ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq > ${sample}_assembled.fasta --memory 256
        //spades.py  -k 21,33 -s ${sample}_mapped_0.fq -s ${sample}_mapped_1.fq -s ${sample}_mapped_2.fq -o .
        // velveth ${sample} 77 -separate -fastq ${sample}_mapped_0.fq ${sample}_mapped_1.fq ${sample}_mapped_2.fq
        // velvetg ${sample} -cov_cutoff auto -exp_cov auto
        // cp ${sample}/contigs.fa ${sample}_assembled.fasta
process translate32 {

    publishDir "$params.output.folder32/translated_sequences/${sample}/", mode : "copy"
    input:
        tuple val(sample), path(assembly) from assemblyout1125899906842624
    
    output:
        tuple val(sample), path("${sample}_translated.fasta") into translatedseq1125899906842624

    script:
        """
       transeq ${assembly} ${sample}_translated.fasta
        """
}


process countpattern32 {

    publishDir "$params.output.folder32/resultfiles/${sample}", mode : "copy"
    input:
        tuple val(sample), path(transfasta) from translatedseq1125899906842624
        path(pyscripts_path) from pyscripts

    output:
        file("${sample}_results.txt")
    
    script:
        """
        python ${pyscripts_path}/hrp2_correctedpath.py ${transfasta} > ${sample}_results.txt
        """
}



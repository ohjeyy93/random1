import subprocess
import csv
import os

fastq_path1 = "fq/minion/FLA/fastq_fail/barcode01"
fastq_path2 = "fq/minion/FLA/fastq_fail/barcode02"
fastq_path3 = "fq/minion/FLA/fastq_fail/barcode11"
fastq_path4 = "fq/minion/FLA/fastq_fail/unclassified"
fastq_path5 = "fq/minion/FLA/fastq_pass/barcode01"
fastq_path6 = "fq/minion/FLA/fastq_pass/barcode02"
fastq_path7 = "fq/minion/FLA/fastq_pass/barcode11"
fastq_path8 = "fq/minion/FLA/fastq_pass/unclassified"
fastq_path9 = "fq/minion/REPA/fastq_fail/barcode01"
fastq_path10 = "fq/minion/REPA/fastq_fail/barcode02"
fastq_path11 = "fq/minion/REPA/fastq_fail/barcode07"
fastq_path12 = "fq/minion/REPA/fastq_fail/barcode09"
fastq_path13 = "fq/minion/REPA/fastq_fail/barcode11"
fastq_path14 = "fq/minion/REPA/fastq_fail/unclassified"
fastq_path15 = "fq/minion/REPA/fastq_pass/barcode01"
fastq_path16 = "fq/minion/REPA/fastq_pass/barcode02"
fastq_path17 = "fq/minion/REPA/fastq_pass/barcode07"
fastq_path18 = "fq/minion/REPA/fastq_pass/barcode09"
fastq_path19 = "fq/minion/REPA/fastq_pass/barcode11"
fastq_path20 = "fq/minion/REPA/fastq_pass/unclassified"
fastq_path21 = "fq/minion/REPB/fastq_fail/barcode01"
fastq_path22 = "fq/minion/REPB/fastq_fail/barcode02"
fastq_path23 = "fq/minion/REPB/fastq_fail/barcode07"
fastq_path24 = "fq/minion/REPB/fastq_fail/barcode09"
fastq_path25 = "fq/minion/REPB/fastq_fail/barcode11"
fastq_path26 = "fq/minion/REPB/fastq_fail/unclassified"
fastq_path27 = "fq/minion/REPB/fastq_pass/barcode01"
fastq_path28 = "fq/minion/REPB/fastq_pass/barcode02"
fastq_path29 = "fq/minion/REPB/fastq_pass/barcode07"
fastq_path30 = "fq/minion/REPB/fastq_pass/barcode09"
fastq_path31 = "fq/minion/REPB/fastq_pass/barcode11"
fastq_path32 = "fq/minion/REPB/fastq_pass/unclassified"

wholedir1=[fastq_path1,fastq_path2,fastq_path3,fastq_path4,fastq_path5,fastq_path6,fastq_path7,
fastq_path8,fastq_path9,fastq_path10,fastq_path11,fastq_path12,fastq_path13,fastq_path14,
fastq_path15,fastq_path16,fastq_path17,fastq_path18,fastq_path19,fastq_path20,fastq_path21,
fastq_path22,fastq_path23,fastq_path24,fastq_path25,fastq_path26,fastq_path27,fastq_path28,
fastq_path29,fastq_path30,fastq_path31,fastq_path32]

#print(wholedir1[31])
#for x in wholedir1:
#    print(len([name for name in os.listdir(x) if os.path.isfile(os.path.join(x, name))]))

with open("test_hrp2multi.nf", "r") as th1:
    with open("test_hrp4multi.nf", "w") as th2:
        count=0
        currentprocess=""
        for line in th1:
            if count==len(wholedir1):
                #print("TRUE")
                break
            if currentprocess=="" and line.find("trim_nanofilt")==-1:
                th2.write(line)
            if line.find("trim_nanofilt")!=-1:
                currentprocess="trim_nanofilt"
            if line.find("maptoreference")!=-1:
                currentprocess="maptoreference"
            if line.find("assemble")!=-1:
                currentprocess="assemble"
            if line.find("translate")!=-1:
                currentprocess="translate"
            if line.find("countpattern")!=-1:
                currentprocess="countpattern"
            if currentprocess=="translate" or currentprocess=="countpattern":
                th2.write(line)
            ###################trim_nanofilt############################
            if currentprocess=="trim_nanofilt" and line.find("${id}_trimmed")==-1 and line.find('\"\"\"')==-1:
                th2.write(line)
            if line.find("output:")!=-1 and currentprocess=="trim_nanofilt":
                #print(count)
                for x in range(len([name for name in os.listdir(wholedir1[count]) if os.path.isfile(os.path.join(wholedir1[count], name))])):
                    th2.write((("\t"+"\t"+"""set val(id),  path("${id}_trimmed_0.fq") into trim_out1"""+"\n").replace("0",str(x))).replace("trim_out1","trim_out"+str(x+1)))
            if line.find("script:")!=-1 and currentprocess=="trim_nanofilt":
                th2.write("\t"+"\t"+'\"\"\"'+"\n")
                for x in range(len([name for name in os.listdir(wholedir1[count]) if os.path.isfile(os.path.join(wholedir1[count], name))])):
                    th2.write((("\t"+"\t"+"""NanoFilt ${id}_0.fastq -l 500 -q 10 > ${id}_trimmed_0.fq"""+"\n").replace("${id}_0","${id}_"+str(x))).replace("${id}_trimmed_0.fq","${id}_trimmed_"+str(x)+".fq"))
                th2.write("\t"+"\t"+'\"\"\"'+"\n")
            #######################maptoreference########################
            if currentprocess=="maptoreference" and line.find("${id}_mapped_")==-1 and line.find('\"\"\"')==-1 and line.find("ref_dir")==-1 and line.find("trim_read")==-1:
                th2.write(line)
            if line.find("input:")!=-1 and currentprocess=="maptoreference":
                for x in range(len([name for name in os.listdir(wholedir1[count]) if os.path.isfile(os.path.join(wholedir1[count], name))])):
                    if x!=len([name for name in os.listdir(wholedir1[count]) if os.path.isfile(os.path.join(wholedir1[count], name))]):
                        th2.write((("\t"+"\t"+"""set val(id), path(trim_read0) from trim_out1"""+"\n").replace("trim_read0","trim_read"+str(x))).replace("trim_out1","trim_out"+str(x+1)))
                    if x==len([name for name in os.listdir(wholedir1[count]) if os.path.isfile(os.path.join(wholedir1[count], name))]):
                        th2.write((("\t"+"\t"+"""path(ref_dir) from ref_dir"""+"\n")))
            if line.find("output:")!=-1 and currentprocess=="maptoreference":
                for x in range(len([name for name in os.listdir(wholedir1[count]) if os.path.isfile(os.path.join(wholedir1[count], name))])):
                    th2.write(((("\t"+"\t"+"""set val(id), path("${id}_mapped_0.fq"), path("${id}_unmapped_0.fq") into mapped_out1"""+"\n").replace("${id}_mapped_0","${id}_mapped_"+str(x))).replace("${id}_unmapped_0","${id}_unmapped_"+str(x))).replace("mapped_out1","mapped_out"+str(x+1)))
            if line.find("script:")!=-1 and currentprocess=="maptoreference":
                th2.write("\t"+"\t"+'\"\"\"'+"\n")
                for x in range(len([name for name in os.listdir(wholedir1[count]) if os.path.isfile(os.path.join(wholedir1[count], name))])):
                    th2.write(((("\t"+"\t"+"""bbduk.sh  -Xmx1g -in=${id}_trimmed_0.fq outm=${id}_mapped_0.fq outu=${id}_unmapped_0.fq ref=$ref_dir/XM_002808697.2.fasta"""+"\n").replace("${id}_trimmed_0","${id}_trimmed_"+str(x))).replace("${id}_mapped_0","${id}_mapped_"+str(x))).replace("${id}_unmapped_0","${id}_unmapped_"+str(x)))
                th2.write("\t"+"\t"+'\"\"\"'+"\n")
            #######################assemble########################
            if currentprocess=="assemble" and line.find("process")!=-1 and line.find("mapped_out")==-1 and line.find('\"\"\"')==-1 and line.find("$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py")==-1:
                th2.write(line)
                th2.write("\t"+"errorStrategy 'ignore'"+"\n")
            if currentprocess=="assemble" and line.find("process")==-1 and line.find("mapped_out")==-1 and line.find('\"\"\"')==-1 and line.find("$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py")==-1:
                th2.write(line)
            if line.find("input:")!=-1 and currentprocess=="assemble":
                for x in range(len([name for name in os.listdir(wholedir1[count]) if os.path.isfile(os.path.join(wholedir1[count], name))])):
                    th2.write(((("\t"+"\t"+"""set val(sample), path(mapped_read0), path(unmapped_read0) from mapped_out1"""+"\n").replace("mapped_read0","mapped_read"+str(x))).replace("unmapped_read0","unmapped_read"+str(x))).replace("mapped_out1","mapped_out"+str(x+1)))
            if line.find("script:")!=-1 and currentprocess=="assemble":
                th2.write("\t"+"\t"+'\"\"\"'+"\n")
                th2.write(("\t"+"\t"+"""$baseDir/SPAdes-3.15.2-Darwin/bin/spades.py"""+" "))
                for x in range(len([name for name in os.listdir(wholedir1[count]) if os.path.isfile(os.path.join(wholedir1[count], name))])):
                    th2.write((""" -s ${sample}_mapped_0.fq """).replace("0",str(x)))
                    if x==len([name for name in os.listdir(wholedir1[count]) if os.path.isfile(os.path.join(wholedir1[count], name))])-1:
                        th2.write("-o . ")
                        th2.write("\n")
                th2.write("\t"+"\t"+'\"\"\"'+"\n")
            #######################################################

            if line.find("params.reads")!=-1 and line.find("params.input.fastq_path")!=-1 and currentprocess!="":
                count+=1

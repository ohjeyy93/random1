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

for x in wholedir1:
    print(len([name for name in os.listdir(x) if os.path.isfile(os.path.join(x, name))]))

with open("test_hrp2multi.nf", "r") as th1:
    with open("test_hrp3multi.nf", "w") as th2:
        count=0
        for line in th1:
            templine=""
            if line.find("params.reads")!=-1 and line.find("${id}_trimmed_")!=-1:
                count+=1
                tempcount=0
                for word in line:
                    templine+=word
                    tempcount+=1
                    if tempcount==line.find("params.reads")+12:
                        templine+=str(count)
                    if tempcount==line.find("params.input.fastq_path")+23:
                        templine+=str(count)
            if line.find("$params.output.folder")!=-1:
                tempcount=0
                for word in line:
                    templine+=word
                    tempcount+=1
                    if tempcount==line.find("$params.output.folder")+21:
                        templine+=str(count)
            if line.find("params.reads")!=-1 and line.find("params.input.fastq_path")==-1:
                tempcount=0
                for word in line:
                    templine+=word
                    tempcount+=1
                    if tempcount==line.find("params.reads")+12:
                        templine+=str(count)
            if line.find("params.genome")!=-1:
                tempcount=0
                for word in line:
                    templine+=word
                    tempcount+=1
                    if tempcount==line.find("params.genome")+13:
                        templine+=str(count)
            if line.find("process trim_nanofilt")!=-1:
                tempcount=0
                for word in line:
                    templine+=word
                    tempcount+=1
                    if tempcount==line.find("process trim_nanofilt")+21:
                        templine+=str(count)
            if line.find("process maptoreference")!=-1:
                tempcount=0
                for word in line:
                    templine+=word
                    tempcount+=1
                    if tempcount==line.find("process maptoreference")+22:
                        templine+=str(count)
            if line.find("process assemble")!=-1:
                tempcount=0
                for word in line:
                    templine+=word
                    tempcount+=1
                    if tempcount==line.find("process assemble")+16:
                        templine+=str(count)
            if line.find("process translate")!=-1:
                tempcount=0
                for word in line:
                    templine+=word
                    tempcount+=1
                    if tempcount==line.find("process translate")+17:
                        templine+=str(count)
            if line.find("process countpattern")!=-1:
                tempcount=0
                for word in line:
                    templine+=word
                    tempcount+=1
                    if tempcount==line.find("process countpattern")+20:
                        templine+=str(count)
            if line.find("trim_out")!=-1:
                tempcount=0
                for word in line:
                    templine+=word
                    tempcount+=1
                    if tempcount==len(line)-1:
                        templine=templine[0:-1]
                        templine+=str(int(word)*count**2)
            if line.find("mapped_out")!=-1:
                tempcount=0
                for word in line:
                    templine+=word
                    tempcount+=1
                    if tempcount==len(line)-1:
                        templine=templine[0:-1]
                        templine+=str(int(word)*count**2)
            if line.find("assemblyout")!=-1:
                tempcount=0
                for word in line:
                    templine+=word
                    tempcount+=1
                    if tempcount==len(line)-1:
                        templine=templine[0:-1]
                        templine+=str(int(word)*count**2)
            if line.find("translatedseq")!=-1:
                tempcount=0
                for word in line:
                    templine+=word
                    tempcount+=1
                    if tempcount==len(line)-1:
                        templine=templine[0:-1]
                        templine+=str(int(word)*count**2)
            if templine=="":
                templine=line
            th2.write(templine)

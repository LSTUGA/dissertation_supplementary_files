#command: python run_tree_pipeline_cfsan_ONT.py raw_reads_illumina/ raw_reads_ONT/ reference_genome.fasta OUTPUT_NAME 8 yes
import os,sys

filepath=os.path.abspath(sys.argv[1]) #directory for Illumina raw reads (e.g. *_1.fastq.gz, *_2.fastq.gz; *1.fq.gz, *2.fq.gz).
filepath_ONT=os.path.abspath(sys.argv[2]) #directory for ONT raw reads (e.g. *_ont.fastq; *porechop.fastq).
reference_file=os.path.abspath(sys.argv[3]) #e.g. Salmonella_enterica_serovar_Enteritidis_P125109_uid30687.fna
project_name=sys.argv[4] #name of the project. e.g. Salmonella_Enteritidis
threads=sys.argv[5] #number of threads, e.g. 10
ref_type=sys.argv[6] #"yes" or "no", if "yes", the reference genome will be included in the final trees.
dirpath = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

# run cfsan pipeline
samples  = [x for x in os.listdir(filepath) if x.endswith("_1.fastq.gz") or x.endswith("1.fq.gz")]
samples_ONT  = [x for x in os.listdir(filepath_ONT) if x.endswith(".fastq")]
cfsan_path = project_name+'_cfsan'
os.system("mkdir -p "+cfsan_path)
os.chdir(cfsan_path)
os.system("mkdir samples")
handle=open('sampleFullPathNames.txt','w')
for sample_r1 in samples:
    sample=sample_r1.split('_')[0]
    os.system("mkdir samples/"+sample)
    if sample_r1.endswith("_1.fastq.gz"):
        sample_r2=sample_r1.replace("_1.fastq.gz","_2.fastq.gz")
    	os.system("ln -sf "+filepath+'/'+sample_r1+" samples/"+sample)
    	os.system("ln -sf "+filepath+'/'+sample_r2+" samples/"+sample)
    	handle.write("samples/"+sample+"/"+sample_r1+" "+"samples/"+sample+"/"+sample_r2+"\n")
    if sample_r1.endswith("1.fq.gz"):
        sample_r2=sample_r1.replace("1.fq.gz","2.fq.gz")
        os.system("ln -sf "+filepath+'/'+sample_r1+" samples/"+sample)
        os.system("ln -sf "+filepath+'/'+sample_r2+" samples/"+sample)
        handle.write("samples/"+sample+"/"+sample_r1+" "+"samples/"+sample+"/"+sample_r2+"\n")
os.system("bwa index "+reference_file)
for sample_ONT in samples_ONT:
    sample=sample_ONT.split('_')[0]
    os.system("mkdir samples/"+sample)
    os.system("ln -sf "+filepath_ONT+'/'+sample_ONT+" samples/"+sample)
    os.system("bwa mem -x ont2d "+reference_file+" samples/"+sample+"/"+sample_ONT+" -t "+threads+" > samples/"+sample+"/reads.sam")
    os.system("samtools view -S -b -F 4 -q 30 --threads "+threads+" -o samples/"+sample+"/reads.unsorted.bam samples/"+sample+"/reads.sam")
    os.system("samtools sort --threads "+threads+" -o samples/"+sample+"/reads.sorted.bam samples/"+sample+"/reads.unsorted.bam")
    os.system("java -jar /home/programs/picard-2.9.2/picard.jar MarkDuplicates INPUT=samples/"+sample+"/reads.sorted.bam OUTPUT=samples/"+sample+"/reads.sorted.deduped.bam METRICS_FILE=samples/"+sample+"/duplicate_reads_metrics.txt")
    os.system("samtools index -@ 8 samples/"+sample+"/reads.sorted.deduped.bam samples/"+sample+"/reads.sorted.deduped.bai")
handle.close()

os.system("ls -d samples/* > sampleDirectories.txt")
os.system("cfsan_snp_pipeline index_ref "+reference_file)
os.environ["Bowtie2Align_ExtraParams"] = "--reorder -X 1000"
os.environ["SamtoolsSamFilter_ExtraParams"] = "-F 4 -q 30"
os.environ["EnableLocalRealignment"] = "false" # disable GATK realignment
os.system('cat sampleFullPathNames.txt | xargs -n 2 -L 1 cfsan_snp_pipeline map_reads --threads '+threads+' '+reference_file)
os.environ["SamtoolsMpileup_ExtraParams"] = "-q 0 -Q 13 -A -B"
os.environ["VarscanMpileup2snp_ExtraParams"] = "--min-var-freq 0.90" # cfsan default
os.system('cat sampleDirectories.txt | xargs -n 1 -P '+threads+' cfsan_snp_pipeline call_sites '+reference_file)

os.system("cfsan_snp_pipeline filter_regions --window_size 1000 125 15 --max_snp 3 2 1 -n var.flt.vcf sampleDirectories.txt "+reference_file)
os.system("cfsan_snp_pipeline merge_sites -f -n var.flt.vcf -o snplist.txt sampleDirectories.txt sampleDirectories.txt.OrigVCF.filtered")
os.system("cfsan_snp_pipeline merge_sites -f -n var.flt_preserved.vcf -o snplist_preserved.txt sampleDirectories.txt sampleDirectories.txt.PresVCF.filtered")
os.system("cat sampleDirectories.txt | xargs -n 1 -P "+threads+" -I XX cfsan_snp_pipeline call_consensus -f -l snplist.txt --minConsDpth 3 --vcfFileName consensus.vcf -o XX/consensus.fasta XX/reads.all.pileup")
os.system("cat sampleDirectories.txt | xargs -n 1 -P "+threads+" -I XX cfsan_snp_pipeline call_consensus -f -l snplist_preserved.txt --minConsDpth 3 --vcfFileName consensus_preserved.vcf -o XX/consensus_preserved.fasta -e XX/var.flt_removed.vcf XX/reads.all.pileup")
os.system("cfsan_snp_pipeline snp_matrix -f -c consensus.fasta -o snpma.fasta sampleDirectories.txt.OrigVCF.filtered")
os.system("cfsan_snp_pipeline snp_matrix -f -c consensus_preserved.fasta -o snpma_preserved.fasta sampleDirectories.txt.PresVCF.filtered")
os.system("cfsan_snp_pipeline snp_reference -f -l snplist.txt -o referenceSNP.fasta "+reference_file)
os.system("cfsan_snp_pipeline snp_reference -f -l snplist_preserved.txt -o referenceSNP_preserved.fasta "+reference_file)
os.system("cat sampleDirectories.txt | xargs -n 1 -P "+threads+" -I XX cfsan_snp_pipeline collect_metrics -f -o XX/metrics XX "+reference_file)
os.system("cfsan_snp_pipeline combine_metrics -f -n metrics -o metrics.tsv sampleDirectories.txt")
os.system("cfsan_snp_pipeline merge_vcfs -f -n consensus.vcf -o snpma.vcf sampleDirectories.txt.OrigVCF.filtered")
os.system("cfsan_snp_pipeline merge_vcfs -f -n consensus_preserved.vcf -o snpma_preserved.vcf sampleDirectories.txt.PresVCF.filtered")
os.system("cfsan_snp_pipeline distance -f -p snp_distance_pairwise.tsv -m snp_distance_matrix.tsv snpma.fasta")
os.system("cfsan_snp_pipeline distance -f -p snp_distance_pairwise_preserved.tsv -m snp_distance_matrix_preserved.tsv snpma_preserved.fasta")
if ref_type=="yes":
    os.system('cat snpma_preserved.fasta referenceSNP_preserved.fasta > '+project_name+'_cfsan.fasta')
elif ref_type=="no":
    os.system('cp snpma_preserved.fasta '+project_name+'_cfsan.fasta')
os.system('python '+dirpath+'/extract_un-gaps.py '+project_name+'_cfsan.fasta '+project_name+'_cfsan -')
os.system(dirpath+'/Fasta2Phylip.pl '+project_name+'_cfsan_concate_snp.fasta '+project_name+'_cfsan_concate_snp.fasta.phylip')
os.system('phyml -q -i '+project_name+'_cfsan_concate_snp.fasta.phylip')

os.system(dirpath+'/Fasta2Phylip.pl snpma_preserved.fasta '+project_name+'_snpma_preserved.fasta.phylip')
os.system('phyml -q -i '+project_name+'_snpma_preserved.fasta.phylip')

#!/usr/bin/env python


# Name:         assembly.py
# Author:       Tom van Wijk - RIVM Bilthoven
# Date:         26-09-2018
# Licence:      GNU General Public License v3.0 (copy provided in directory)

# For detailed information and instruction on how to install and use this software
# please view the included "READNE.md" file in this repository


# import python libraries
from argparse import ArgumentParser
import logging
import logging.handlers
import os
import sys
import random


# Function to parse the command-line arguments
# Returns a namespace with argument keys and values
def parse_arguments(args, log):
	log.info("Parsing command line arguments...")
	parser = ArgumentParser(prog="assembly.py")
	parser.add_argument("-i", "--indir", dest = "input_dir",
		action = "store", default = None, type = str,
		help = "Location of input directory (required)",
		required = True)
	parser.add_argument("-o", "--outdir", dest = "output_dir",
		action = "store", default = "inputdir", type = str,
		help = "Location of output directory (default=inputdir)")
	parser.add_argument("-t", "--threads", dest = "threads",
		action = "store", default = 4, type = int,
		help = "Number of threads to be used (default=4)")
	parser.add_argument("-m", "--memory", dest = "ram",
		action = "store", default = 13, type = int,
		help = "Maximum amount of RAM (GB) to be used (default=13)")
        parser.add_argument("-x", "--savetemp", dest = "save_temp",
                action = "store", default = "false", type = str,
                help = "Option to save temporary files (default=false)")
	return parser.parse_args()


# Function creates logger with handlers for both logfile and console output
# Returns logger
def create_logger(logid):
        # create logger
        log = logging.getLogger()
        log.setLevel(logging.INFO)
        # create file handler
        fh = logging.FileHandler(str(logid)+'_assembly.log')
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter('%(message)s'))
        log.addHandler(fh)
        # create console handler
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(logging.Formatter('%(message)s'))
        log.addHandler(ch)
        return log


# Function for trimming reads at q=25 using erne-filter
# Writes erne-filter output files to output_dir
def trim_reads(R1_file, R2_file, output_dir, prefix, threads, log):
	# Running the actual trimming for each sample using enre-filter
	log.info("Trimming reads...")
	os.system("erne-filter --query1 %s --query2 %s --output-prefix %s --q 25 --threads %s"
			  % (R1_file, R2_file, output_dir+"/"+prefix, str(threads)))


# Function creates a list of files or directories in <inputdir>
# on the specified directory depth
def list_directory(input_dir, obj_type, depth):
        dir_depth = 1
        for root, dirs, files in os.walk(input_dir):
                if dir_depth == depth:
                        if obj_type ==  'files':
                                return files
                        elif obj_type == 'dirs':
                                return dirs
                dir_depth += 1


# Function to perform no-novo assembly of reads to contigs using SPAdes
# Writes spades output to temp_dir and moves relevant files to output_dir
def assemble_reads(prefix, input_dir, temp_dir, output_dir, logs_dir, threads, ram, log):
	log.info("Assembling reads using SPAdes...")
	# Make list of files in input_dir
	files = list_directory(input_dir, 'files', 1)
	print files
	if prefix+"_1.fastq" in files and prefix+"_2.fastq" in files:
		# building SPAdes command
		spades_command = ("spades.py --pe1-1 %s --pe1-2 %s"
			% (input_dir+"/"+prefix+"_1.fastq", input_dir+"/"+prefix+"_2.fastq"))
		# include orphaned sequences if available
		if prefix+"_unpaired.fastq" in files:
			spades_command += (" --pe1-s "+input_dir+"/"+prefix+"_unpaired.fastq")
		spades_command += (" -o %s --threads %s --memory %s" % (temp_dir, str(threads), str(ram)))
		# run SPAdes
		log.info("running spades with command:\n"+spades_command)
		os.system(spades_command)
		# Move relevant output from temp_dir to output_dir
		os.system("mv "+temp_dir+"/scaffolds.fasta "+output_dir+"/"+prefix+'_assembly.fasta')
		os.system("mv "+temp_dir+"/contigs.paths "+output_dir+"/metadata/"+prefix+'_contigs.paths')
		os.system("mv "+temp_dir+"/scaffolds.paths "+output_dir+"/metadata/"+prefix+'_scaffolds.paths')
		os.system("mv "+temp_dir+"/spades.log "+logs_dir+"/"+prefix+'_spades.log')
		os.system("mv "+temp_dir+"/params.txt "+output_dir+"/metadata/"+prefix+'_params.txt')
	else:
		log.warning("Could not find at least one of the following files:\n"
			+prefix+"_1.fastq\n"+prefix+"_2.fastq\nSkip sample...")


# Function to perform a QUAST quality assessment on the fasta files in the given directory
# Writes QUAST output do output_dir
def assess_assemblies(input_dir, output_dir, threads, log):
	log.info("Creating quality rapport of assemblies using QUAST...")
	# Iterate over files in directory
	file_string = ""
	for file in list_directory(input_dir, 'files', 1):
		file_string+=(input_dir+"/"+file+" ")
	# Running  QUAST
	os.system("quast.py -t "+str(threads)+" -o "+output_dir+" "+file_string)
	# Cleaning up QUAST output directory


# Function closes logger handlers
def close_logger(log):
        for handler in log.handlers:
                handler.close()
                log.removeFilter(handler)


# MAIN function
def main():
        # create logger
        logid = random.randint(99999, 9999999)
        log = create_logger(logid)
        # parse command line arguments
        args = parse_arguments(sys.argv, log)
        # creating output directory
        if args.output_dir == 'inputdir':
                outdir = os.path.abspath(args.input_dir+"/assembly_pipeline_output")
        else:
                outdir = os.path.abspath(args.output_dir)
        log.info("output directory: "+outdir)
	# create output dir and sub-dirs
	os.system("mkdir -p "+outdir+"/temp/trimmed_reads "
		+outdir+"/temp/spades "
		+outdir+"/logfiles/spades "
		+outdir+"/quality_control_reads/fastqc "
		+outdir+"/quality_control_assemblies "
		+outdir+"/assemblies/metadata ")
		# iterate over filepairs
	filelist = list_directory(args.input_dir, "files", 1)
	for file in filelist:
		if '_R1' in file and file in filelist and file.replace('_R1', '_R2') in filelist:
			R1_file = args.input_dir+"/"+file
			R2_file = args.input_dir+"/"+file.replace('_R1', '_R2')
			prefix = file.split("_R1")[0].replace(".","_").replace(" ","_")+"_trimmed"
			log.info("processing filepair:\t"+R1_file+",\t"+R2_file)
			# create a fastCQ quality report
			os.system("fastqc "+R1_file+" "+R2_file+" -o "+outdir+"/quality_control_reads/fastqc")
			# Perform quality trimming the raw sequences
			erne_out = outdir+"/temp/trimmed_reads"
			trim_reads(R1_file, R2_file, erne_out, prefix, args.threads, log)
			# Perform de novo assembly of reads
			spades_out = outdir+"/temp/spades/"+prefix
			os.system("mkdir -p "+spades_out)
			assemble_reads(prefix, erne_out, spades_out, outdir+"/assemblies",
				outdir+"/logfiles/spades", args.threads, args.ram, log)
			# empty /outdir/temp if param save_temp = 'false'
#			if args.save_temp == "false":
#				log.info("Emptying "+outdir+"/temp to save storage space...")
#				os.system("rm -rf "+outdir+"/temp")
#				os.system("mkdir -p "+outdir+"/temp/spades")
	# quality assessment of assemblies
	quast_out = outdir+"/quality_control_assemblies"
	assess_assemblies(outdir+"/assemblies", quast_out, args.threads, log)
	# unzipping fastqc reports
	files = list_directory(outdir+"/quality_control_reads/fastqc", "files", 1)
	for file in files:
		if file.endswith(".zip"):
			os.system("unzip "+outdir+"/quality_control_reads/fastqc/"+file
				+" -d "+outdir+"/quality_control_reads/fastqc/"+file.replace('.zip', ''))
	# create multiqc report
	os.system("multiqc "+outdir+"/quality_control_reads/fastqc/ -o "+outdir+"/quality_control_reads")
	# Remove outdir/temp directory if temp parameter is not set to anything but false
	if args.save_temp == "false":
		log.info("Emptying "+outdir+"/temp to save storage space...")
		os.system("rm -rf "+outdir+"/temp")
	else:
		log.info("Keeping temporary files and directories...")
	# close logger handlers
        log.info("\nClosing logger and finalising bacterial_genome_assembler.py")
        close_logger(log)
	# move logfile to output directory
        os.system("mv "+str(logid)+"_assembly.log "+outdir+"/logfiles/assembly.log")


main()

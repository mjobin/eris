#!/usr/bin/python

#####################################
#####        HPG Lab            #####
#####    updated March 2018     #####
#####       MJJ                 #####
#####################################

# Script for simulating fragmentation and damage to DNA
# including aDNA alterations. It then packs the fragments
# into randomly sized FASTA files.

import argparse
from argparse import RawTextHelpFormatter
import os
import sys
import progressbar
import datetime
import gzip
import random
import shutil
import fnmatch
import subprocess
from subprocess import Popen, PIPE
import numpy as np
from Bio import SeqIO

def bash_command(cmd):
    # print cmd
    cmdfile.write(cmd)
    cmdfile.write("\n\n")
    subp = subprocess.Popen(['/bin/bash', '-c', cmd], stdout=PIPE, stderr=PIPE)
    stdout, stderr = subp.communicate()
    return stdout

def bwa_align_se(infile, ref, refbase, outdir):
    infilebase = os.path.basename(infile).split(".")[0]
    bwaalignline = "bwa aln -l " + seed_disable + " -n " + bwamaxedit + " -t " + threads + " " + ref + " " + infile
    bwasline = "bwa samse " + ref + " - " + infile
    samviewline = "samtools view " + " -F 4 -buh - "
    samsortline = "samtools sort -@ " + threads + " -m 4G -o " + outdir + "/" + infilebase + "-" + refbase + ".bam"
    seline = bwaalignline + " | " + bwasline + " | " + samviewline + " | " + samsortline
    bash_command(seline)

def blastn():
    blastline = "blastn -word_size " + wordsize + " -perc_identity 90 -evalue 0.0000001 -num_threads " + threads + " -db nt -query $i -out ${i//.fasta/.st-nt.blast}"


def chk_mkdir(dname):
    if overwrite:
        if os.path.exists(dname):
            shutil.rmtree(dname)
        os.mkdir(dname)
    else:
        if os.path.exists(dname):
            print "ERROR: Directory " + dname + " exists and overwrite not set to true. Exiting."
            exit(1)
        else:
            os.mkdir(dname)



if __name__ == "__main__":

    print "\n****************\nERIS\n****************\n"

    parser = argparse.ArgumentParser(description="# This script:\n"
                                                    "1. loads a dir full of genomes\n"
                                                    "2. cuts them in stepped sizes or normally distrubuted sizes\n"
                                                    "3. applies simulated aDNA damage\n"
                                                    "4. packs them into random FASTA files\n"
                                                 "- ", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-genomedir', metavar='<genomedir>', help='location of genomes.', default='/data/genomes/eris/')
    parser.add_argument('-wd', metavar='<wd>', help='Working directory. Defaults to current.', default='.')
    parser.add_argument('-verbose', dest='verbose', help='Print stdout and stderr to console.',
                        action='store_true')
    parser.set_defaults(verbose=False)
    parser.add_argument('-nfiles', metavar='<nfiles>', help='Number of FASTA outputfiles',
                        default=100)
    parser.add_argument('-nlmu', metavar='<nlmu>', help='Mean numbers of lines in any one FASTA file',
                        default=100000)
    parser.add_argument('-nlsigma', metavar='<nlsigma>', help='Variance in numbers of lines in any one FASTA file',
                        default=5)
    parser.add_argument('-fmode', metavar='<fmode>', help='fixed, normal',
                        default='fixed')
    parser.add_argument('-fixedmin', metavar='<fixedmin>', help='Fixed sequence minimum length.',
                        default=20)
    parser.add_argument('-fixedgap', metavar='<fixedgap>', help='Fixed sequence length gap.',
                        default=20)
    parser.add_argument('-fixednum', metavar='<fixednum>', help='Fixed sequence number of stepped sequence lengths.',
                        default=10)
    parser.add_argument('-nmu', metavar='<nmu>', help='Mean size of fragments',
                        default=30)
    parser.add_argument('-nsigma', metavar='<nsigma>', help='Variance of sizes of fragments',
                        default=10)
    parser.add_argument('-sigmas', metavar='<sigmas>', help='cytosine de- amination rate for single stranded',
                        default=0.05)
    parser.add_argument('-pmdp', metavar='<pmdp>', help='PMDTools p factor',
                        default=0.3)
    parser.add_argument('-pmdc', metavar='<pmdc>', help='PMDTools C copnstant',
                        default=0.01)
    parser.add_argument('-bwamaxedit', metavar='<bwamaxedit>', help='Maximum edit distance if the value is INT, or the fraction of missing alignments given 2 percent uniform base error rate if FLOAT. In the latter case, the maximum edit distance is automatically chosen for different read lengths.',
                        default="0.01")
    parser.add_argument('-seed_disable', metavar='<seed_disable>', help='Following ancient DNA data processing protocols',
                        default="1024")
    parser.add_argument('-threads', metavar='<threads>', help='Number of concerrent threads to use when invoking software.',
                        default="23")
    parser.add_argument('-q', metavar='<q>', help='BWA min quality. 20 provides a fairly low cutoff',
                        default="20")
    parser.add_argument('-refs', dest='refs', nargs='+', default=[],
                        help='List of reference sequences to map against.')
    parser.add_argument('-genomes', dest='genomes', nargs='+', default=[],
                        help='List of source genomes.')
    parser.add_argument('-wordsize', metavar='<wordsize>', help='BlastN wordsize.',
                        default="15")
    parser.add_argument('-perc_identity', metavar='<perc_identity>', help='BlastN perc_identity.',
                        default="90")
    parser.add_argument('-adnadmg', dest='adnadmg', help='Turn aDNA damage on.',
                        action='store_true')
    parser.set_defaults(adnadmg=False)
    parser.add_argument('-overwrite', dest='overwrite', help='Overwrite existing files and directories.',
                        action='store_true')
    parser.set_defaults(overwrite=False)
    parser.add_argument('-bwaindex', dest='bwaindex', help='Need to index if never used the reference genome before.',
                        action='store_true')
    parser.set_defaults(bwaindex=False)
    parser.add_argument('-index_algorithm', metavar='<index_algorithm>',
                        help='If reference is <2Gb use is, if >2Gb use bwtsw',
                        default='is')
    parser.add_argument('-glm', dest='glm', help='Override nlmean with 4.6*mean genome length / mean fragment length.',
                        action='store_true')
    parser.set_defaults(glm=False)




    args = parser.parse_args()
    genomedir = args.genomedir
    wd = args.wd
    nfiles = int(args.nfiles)
    nlmu = int(args.nlmu)
    nlsigma = int(args.nlsigma)
    fmode = args.fmode
    fixedmin = int(args.fixedmin)
    fixedgap = int(args.fixedgap)
    fixednum = int(args.fixednum)
    nmu = int(args.nmu)
    nsigma = int(args.nsigma)
    sigmas = float(args.sigmas)
    pmdp = float(args.pmdp)
    pmdc = float(args.pmdc)
    seed_disable = args.seed_disable
    threads = args.threads
    bwamaxedit = args.bwamaxedit
    verbose = bool(args.verbose)
    refs = args.refs
    genomes = args.genomes
    wordsize = args.wordsize
    perc_identity = args.perc_identity
    adnadmg = bool(args.adnadmg)
    overwrite = bool(args.overwrite)
    bwaindex = bool(args.bwaindex)
    index_algorithm = args.index_algorithm
    glm = bool(args.glm)

    os.chdir(wd)
    cwd = os.getcwd()
    print "Working in: ", cwd

    rng = random.SystemRandom()  # Uses /dev/urandom

    cmdfile = open("eris_cmds", 'w')

    today = datetime.date.today()
    logfilename = wd + "/out.eris." + str(today) + ".log"
    print "Logging to: ", logfilename
    logfile = open(logfilename, 'w')

    gdict = {}
    outdict = {}
    infilelist = []

    logfile.write("\n****************\nERIS\n****************\n")
    logfile.write("Run started: " + str(today) + "\n")

    logfile.write("Number of arguments:" + str(len(sys.argv)) + "arguments.")
    logfile.write("Argument List:")
    logfile.write(str(sys.argv))

    logfile.write("aDNA damage: ")
    if adnadmg:
        logfile.write("ON\n")
        logfile.write("aDNA damage parametqqers:")
        logfile.write("PMDTools p factor:" + str(pmdp))
        logfile.write(("PMDTools c constant: " + str(pmdc)))
    else:
        logfile.write("OFF\n")
    logfile.write("-------------------\n\n")

    logfile.write("Original genomes used:\n")

    if genomes:
        for genome in genomes:
            logfile.write(genome)
            logfile.write("\n")
            infilelist.append(genome)
    else:
        for gfilename in os.listdir(genomedir):
            if gfilename.endswith(".fa") or gfilename.endswith(".fasta") or gfilename.endswith(".fas"):
                if " " in gfilename:
                    print "ERROR: spaces not allowed in filename: " + gfilename
                logfile.write(os.path.join(genomedir, gfilename))
                logfile.write("\n")
                infilelist.append(os.path.join(genomedir, gfilename))

    reflist = []
    if refs:
        for ref in refs:
            reflist.append(ref)
    else:
        for infile in infilelist:
            reflist.append(infile)


    logfile.write("-------------------\n\n")
    logfile.write("Reference sequences mapped: \n")
    for ref in reflist:
        refname = os.path.basename(ref)
        filebase, fileext = os.path.splitext(refname)
        logfile.write(ref + "\n")
    logfile.write("-------------------\n\n")

    # logfile.write("-------------------\n\n")
    # refdic = {}
    # logfile.write("Reference sequences mapped: \n")
    # logfile.write("Self (always on)\n")
    # for ref in arefs:
    #     refname = os.path.basename(ref)
    #     filebase, fileext = os.path.splitext(refname)
    #     refdic[filebase] = ref
    #     logfile.write(ref + "\n")
    # logfile.write("-------------------\n\n")


    ########### BWA - aligning reads to a reference sequence ############
    if bwaindex:
        print "\nIndexing genomes..."
        for infilename in infilelist:
            print "Indexing " + infilename + "...\n"
            bash_command("bwa index -p " + infilename + " -a " + index_algorithm + " " + infilename)

    glengths = []
    print "\nLoading input genomes..."
    for infilename in infilelist:
        filename = os.path.basename(infilename)
        filebase, fileext = os.path.splitext(filename)

        fasta_sequences = SeqIO.parse(open(infilename), 'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq.upper())
            fullname = filebase
            print fullname + "    " +  str(len(sequence))
            glengths.append(len(sequence))
            gdict[fullname] = sequence


    meanglength = np.mean(glengths)

    if glm:
        logfile.write("Mean numbers of lines in any one FASTA file: " + str(4.6 * meanglength) +" / mean length of fragment\n")
        logfile.write("Variance in numbers of lines in any one FASTA file: 0.0 \n")
    else:
        logfile.write("Mean numbers of lines in any one FASTA file: " + str(nlmu)+"\n")
        logfile.write("Variance in numbers of lines in any one FASTA file: " + str(nlsigma)+"\n")

    # for infilename in infilelist:
    #     filename = os.path.basename(infilename)
    #     filebase, fileext = os.path.splitext(filename)
    #     filebase = filebase.replace(' ', '_')
    #
    #
    #     chk_mkdir(filebase)
    #     print "Processing " + filebase + "...\n"




    sublist = []

    if fmode == 'normal':

       sublist.append(-1)
    elif fmode == 'fixed':
        fixedsize = fixedmin
        for i in xrange(fixednum):
            sublist.append(fixedsize)
            fixedsize += fixedgap

    for sub in sublist:
        slen = sub
        newdir = str(sub) + "bp"
        if sub < 1: #normal
            newdir = str(nmu) + "-mu-"+str(nsigma)+"-sig-"
            slen = np.random.normal(nmu, nsigma)
        print "\nRandomizing files in: " + newdir

        chk_mkdir(newdir)

        if glm:
            nlmu = int(4.6 * meanglength / slen) #Lander-waterman 99% coverage
            nlsigma = 0.0

        onames = []
        for ofile in range(nfiles):
            ofilename = newdir + "/erisout-" + str(ofile) + ".fasta"
            onames.append(ofilename)
            outfile = open(ofilename, 'w')
            outs = int(np.random.normal(nlmu, nlsigma))
            for i in range(outs):
                oname = rng.choice(gdict.keys())
                oseq = gdict[oname]
                slen = sub
                if fmode == 'normal':
                    slen = np.random.normal(nmu, nsigma)

                seqloc = rng.randint(0, (len(oseq) - 1))

                five = 0
                three = 0
                if rng.randint(0, 1) == 0:  # put down minumum of range first
                    five = seqloc
                    three = len(oseq) - 1
                    if five + slen < oseq[seqloc]:
                        three = five + slen
                else:  # put down max of range first
                    three = seqloc
                    five = 0
                    if three - slen > 0:
                        five = three - slen
                sample = oseq[int(five):int(three)]

                # And now apply damage
                slist = list(sample)
                if adnadmg:

                    for j in range(len(slist)):
                        z = -1
                        if slist[j] == 'G':  # distance from 3' end
                            z = len(slist) - j - 1
                        elif slist[j] == 'C':  # distance from 5' end
                            z = j
                        if z >= 0:
                            dz = (pow((1 - pmdp), (z - 1)) * pmdp) + pmdc
                            rando = rng.random()
                            if rando < dz:
                                if slist[j] == 'G':
                                    slist[j] = 'A'
                                elif slist[j] == 'C':
                                    slist[j] = 'T'

                newsample = ''.join(slist)

                outfile.write(">")
                outfile.write(oname)
                outfile.write("\n")
                outfile.write(newsample)
                outfile.write("\n")
            outfile.close()

        for ref in reflist:
            refbase = os.path.basename(ref).split(".")[0]
            refbase = refbase.replace(' ', '_')
            refdir = newdir + "/" + refbase
            chk_mkdir(refdir)
            print "Mapping to: " + ref
            bar = progressbar.ProgressBar()
            for i in bar(range(len(onames))):
                bwa_align_se(onames[i], ref, refbase, refdir)
                infilebase = os.path.basename(onames[i]).split(".")[0]
                namefinalbase = refdir + "/" + infilebase + "-" + refbase
                bash_command(
                    "bedtools bamtofastq -i " + namefinalbase + ".bam -fq " + namefinalbase + ".aln.fastq")
                bash_command(
                    "seqtk seq -a " + namefinalbase + ".aln.fastq  > " + namefinalbase + ".aln.fasta")

    logfile.close()
    cmdfile.close()
    print "eris.py complete."
    exit(0)
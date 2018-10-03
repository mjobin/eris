#!/usr/bin/python

#####################################
#####        HPG Lab            #####
#####    updated March 2018     #####
#####       MJJ                 #####
#####################################

# Script for s

import argparse
import os
import progressbar
import gzip
import shutil
import subprocess
from subprocess import PIPE
import datetime
import numpy as np
import fnmatch

def bash_command(cmd):
    # print cmd
    cmdfile.write(cmd)
    cmdfile.write("\n\n")
    subp = subprocess.Popen(['/bin/bash', '-c', cmd], stdout=PIPE, stderr=PIPE)
    stdout, stderr = subp.communicate()
    if verbose:
        print stdout
    logfile.write(stdout)
    if verbose:
        print stderr
    logfile.write(stderr)
    return stdout


def nomosnew(statsfile, subdir):
    sizedirs = next(os.walk(wd))[1]
    for sizedir in sizedirs:
        print "\nSize : " + sizedir

        fragsize = sizedir
        if fragsize.endswith('bp'):
            fragsize = fragsize[:-2]

        cursizedir = os.getcwd()

        os.chdir(cursizedir + "/" + sizedir)

        origfiles = []

        curdir = os.getcwd()
        rawfiles = os.listdir(curdir)
        for rawfile in rawfiles:
            if rawfile.endswith(".fasta"):
                origfiles.append(rawfile)
            elif rawfile.endswith(".fasta.gz"):
                fastaname = os.path.splitext(rawfile)[0]
                with gzip.open(rawfile, 'rb') as f_in, open(fastaname, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                flist.append(fastaname)

        refdirs = next(os.walk(curdir))[1]

        for refdir in refdirs:
            print "\nReference: " + refdir
            countfile = open(refdir + "/" + refdir + "counts.csv", 'w')
            countfile.write("Orig file, Orig Count,  Ref File, Ref Count, Correctly Mapped Count\n")
            origcounts = []
            origmatches = []
            outgcounts = []
            refmatchcounts = []

            origmeanlengths = []
            origmeanmatchlengths = []
            outgmeanlengths = []
            corrmeanlengths = []
            incorrmeanlengths = []


            print "Counting..."
            bar = progressbar.ProgressBar()
            for i in bar(range(len(origfiles))):
            # for i in range(len(origfiles)):
                origfilename = origfiles[i]
                origbase, origext = os.path.splitext(origfilename)

                outgfilename = refdir + "/" + origbase + "-" + refdir + ".aln.fasta"

                if os.path.isfile(outgfilename):
                    origlengths = []
                    origmatchlengths = []
                    corrlengths = []
                    incorrlengths = []
                    origrefmatchcount = 0

                    # oldorigout = bash_command("wc -l " + origfilename)
                    # oldorigcols = oldorigout.split()
                    # oldorigcount = int(oldorigcols[0]) / 2

                    origcount = 0
                    origfile = open(origfilename, 'r')
                    for origline in origfile:
                        origline = origline.strip()
                        if origline.startswith(">"):
                            origcount = origcount + 1
                            seqline = origfile.next()
                            origlengths.append(len(seqline))
                            if origline[1:] == refdir:
                                origrefmatchcount = origrefmatchcount + 1
                                origmatchlengths.append(len(seqline))

                    origcounts.append(origcount)
                    # print str(oldorigcount) + "vs" + str(origcount)
                    origmatches.append(origrefmatchcount)


                    outgcount = 0
                    refmatches = 0
                    oglengths = []
                    outgfile = open(outgfilename, 'r')
                    for outgline in outgfile:
                        ogline = outgline.strip()
                        if ogline.startswith(">"):
                            outgcount = outgcount + 1
                            seqline = outgfile.next()
                            oglengths.append(len(seqline))
                            if ogline[1:] == refdir:
                                refmatches = refmatches + 1
                                corrlengths.append(len(seqline))
                            else:
                                incorrlengths.append(len(seqline))
                    outgcounts.append(outgcount)
                    refmatchcounts.append(refmatches)


                    countfile.write(origfilename)
                    countfile.write(",")
                    countfile.write(str(origcount))
                    countfile.write(",")
                    countfile.write(outgfilename)
                    countfile.write(",")
                    countfile.write(str(outgcount))
                    countfile.write(",")
                    countfile.write(str(refmatches))
                    countfile.write("\n")
                    origmeanlengths.append(np.mean(origlengths))
                    origmeanmatchlengths.append(np.mean(origmatchlengths))
                    outgmeanlengths.append(np.mean(oglengths))
                    corrmeanlengths.append(np.mean(corrlengths))
                    incorrmeanlengths.append(np.mean(incorrlengths))
            countfile.close()

            statsfile.write(subdir)
            statsfile.write(",")
            statsfile.write(fragsize)
            statsfile.write(",")
            statsfile.write(refdir)
            statsfile.write(",")
            statsfile.write(str(np.mean(origcounts)))
            statsfile.write(",")
            statsfile.write(str(np.std(origcounts)))
            statsfile.write(",")
            statsfile.write(str(np.mean(origmeanlengths)))
            statsfile.write(",")
            statsfile.write(str(np.mean(origmatches)))
            statsfile.write(",")
            statsfile.write(str(np.std(origmatches)))
            statsfile.write(",")
            statsfile.write(str(np.mean(origmeanmatchlengths)))
            statsfile.write(",")
            statsfile.write(str(np.mean(outgcounts)))
            statsfile.write(",")
            statsfile.write(str(np.std(outgcounts)))
            statsfile.write(",")
            statsfile.write(str(np.mean(outgmeanlengths)))
            statsfile.write(",")
            statsfile.write(str(np.mean(refmatchcounts)))
            statsfile.write(",")
            statsfile.write(str(np.std(refmatchcounts)))
            statsfile.write(",")
            statsfile.write(str(np.mean(corrmeanlengths)))
            statsfile.write(",")
            statsfile.write(str(np.mean(incorrmeanlengths)))
            statsfile.write("\n")


        if kraken:
            print "RELEASE THE KRAKEN!!!!"
            bar = progressbar.ProgressBar()
            for i in bar(range(len(origfiles))):
                origfile = origfiles[i]

                bash_command(
                    "kraken --threads " + threads + " --db " + krakendb + " " + origfile + " > " + origfile + ".kraken")

                bash_command("kraken-translate --db " + krakendb + " " + origfile + ".kraken > " + origfile + ".labels")

                bash_command(
                    "kraken-translate --db " + krakendb + " " + origfile + ".kraken| cut -f2 | python /data/scripts/make_counts.py > " + origfile + ".kraken4krona")

                bash_command("ktImportText " + origfile + ".kraken4krona -o " + origfile + ".kraken.krona.html")

        if diamond:
            print "\nDIAMOND metagenome analysis, mapping fasta reads to NR NCBI database..."


            bar = progressbar.ProgressBar()
            for i in bar(range(len(origfiles))):
                origfile = origfiles[i]

                bash_command("diamond blastx -p " + threads + " -t " + diamondtmp + " -b " + diamondblock + " -q " + origfile + " -d " + blastdir + "/nr.dmnd -o " + origfile + ".dmnd.matches.txt")

                # 	### turn diamond results into rma files viewable in MEGAN
                bash_command(megandir + "/tools/blast2rma -i " + origfile + ".dmnd.matches.txt -f BlastTab -r " + origfile + " -o " + origfile + ".dmnd.matches.rma -g2t " + blastdir + "/gi2tax-July2016.bin")

                # 	### turn diamond output into Krona html file for visualization
                bash_command("ktImportBLAST " + origfile + ".dmnd.matches.txt -o " + origfile + ".dmnd.krona.html")

        os.chdir("..")
    os.chdir("..")


if __name__ == "__main__":

    print "\n****************\nNOMOS\n****************\n"

    parser = argparse.ArgumentParser(description="# This script:\n"

                                                 "- ", formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument('-wd', metavar='<wd>', help='Working directory. Defaults to current.', default='.')
    parser.add_argument('-krakendb', metavar='<krakendb>', help='KrakenDB location',
                        default='/data/db/krakenDB')
    parser.add_argument('-blastdir', metavar='<blastdir>', help='BLAST db directory',
                        default='/data/db/BLAST')
    parser.add_argument('-scriptsdir', metavar='<scriptsdir>', help='Default scripts directory',
                        default='/data/scripts')
    parser.add_argument('-threads', metavar='<threads>',
                        help='Number of concerrent threads to use when invoking software.',
                        default="23")
    parser.add_argument('-verbose', dest='verbose', help='Print stdout and stderr to console.',
                        action='store_true')
    parser.set_defaults(verbose=False)
    parser.add_argument('-dam', dest='dam', help='Search for files with the .dam. extension. Otherwise search for .wodam.',
                        action='store_true')
    parser.set_defaults(dam=False)
    parser.add_argument('-oldstyle', dest='oldstyle', help='Eris was run under old style',
                        action='store_true')
    parser.set_defaults(oldstyle=False)

    parser.add_argument('-kraken', dest='kraken', help='Run Kraken',
                        action='store_true')
    parser.set_defaults(kraken=False)
    parser.add_argument('-multisub', dest='multisub', help='True if multiple sub directories, each their ther own subdirs for the fragment sizes (i.e. 20bp, 40bp...)',
                        action='store_true')
    parser.set_defaults(multisub=False)
    parser.add_argument('-diamond', dest='diamond', help='Run Diamond',
                        action='store_true')
    parser.set_defaults(diamond=False)
    parser.add_argument('-diamondblock', metavar='<diamondblock>', help="Diamond block size",
                        default='48.0')
    parser.add_argument('-diamondtmp', metavar='<diamondtmp>', help="Temp directory for Diamond",
                        default='/dev/shm')
    parser.add_argument('-megandir', metavar='<megandir>', help='MEGAN6 directory',
                        default='/opt/megan')





    args = parser.parse_args()
    wd = args.wd
    krakendb = args.krakendb
    blastdir = args.blastdir
    scriptsdir = args.scriptsdir
    threads = args.threads
    verbose = bool(args.verbose)
    dam = args.dam
    oldstyle = bool(args.oldstyle)
    kraken = bool(args.kraken)
    diamond = bool(args.diamond)
    diamondblock = args.diamondblock
    diamondtmp = args.diamondtmp
    megandir = args.megandir
    multisub = bool(args.multisub)

    cmdfile = open("nomos_cmds", 'w')

    today = datetime.date.today()
    logfilename = wd + "/out.nomos." + str(today) + ".log"
    print "Logging to: ", logfilename
    logfile = open(logfilename, 'w')



    if oldstyle:


        subdirs = next(os.walk(wd))[1]
        for subdir in subdirs:

            os.chdir(wd + "/" + subdir)
            print "**********Now working in " + os.getcwd()

            statsfile = open(subdir + "-stats.csv", 'w')
            statsfile.write("Subdir, Ref,Orig. Mean.,Orig. StD., Ref Mean, Ref StD\n")

            subsubdirs = next(os.walk(wd))[1]
            for subsubdir in subsubdirs:
                os.chdir("./" + subsubdir)
                print "Now working in " + os.getcwd()

                flist = []
                print "\nUncompressing..."
                wdfiles = os.listdir(wd)
                bar = progressbar.ProgressBar()
                for i in bar(range(len(wdfiles))):
                    wdfile = wdfiles[i]
                    if dam:
                        if wdfile.endswith(".dam.fasta"):
                            flist.append(wdfile)
                        elif wdfile.endswith(".dam.fasta.gz"):
                            fastaname = os.path.splitext(wdfile)[0]
                            with gzip.open(wdfile, 'rb') as f_in, open(fastaname, 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                            flist.append(fastaname)
                    else:
                        if wdfile.endswith(".wodam.fasta"):
                            flist.append(wdfile)
                        elif wdfile.endswith(".wodam.fasta.gz"):
                            fastaname = os.path.splitext(wdfile)[0]
                            with gzip.open(wdfile, 'rb') as f_in, open(fastaname, 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                            flist.append(fastaname)




                print "Counting..."

                origcounts = []
                outgcounts = []
                bar = progressbar.ProgressBar()
                for i in bar(range(len(flist))):
                    ffile = flist[i]

                    fname, fext = os.path.splitext(ffile)



                    file = None
                    outgpattern = ""
                    if dam:
                        files = os.listdir("bwa_dam")
                    else:
                        files = os.listdir("bwa")

                    outgpattern += fname
                    outgpattern += ".MTb.fasta*"

                    # if dam:
                    #     outgpattern += ".dam"
                    # else:
                    #     outgpattern += ".wodam"

                    # outgpattern += ".fasta*"


                    outgname = None

                    for name in files:

                        if (fnmatch.fnmatch(name, outgpattern)):
                            outgname = name
                            continue

                    outgfile = ""
                    if dam:
                        outgfile += "bwa_dam/"
                    else:
                        outgfile += "bwa/"

                    outgfile += outgname


                    # print fname + " " + outgfile
                    if not os.path.isfile(ffile):
                        print "ERROR: Cannot find " + ffile
                        exit(1)
                    if not os.path.isfile(outgfile):
                        print "ERROR: Cannot find " + outgfile
                        continue


                    if outgfile.endswith(".fasta.gz"):
                        fastaname = os.path.splitext(outgfile)[0]
                        with gzip.open(outgfile, 'rb') as f_in, open(fastaname, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                        outgfile = fastaname

                    origout = bash_command("wc -l " + ffile)

                    origcols = origout.split()
                    origcount = int(origcols[0]) / 2
                    # print "Wat " +str(origcount)

                    outgout = bash_command("wc -l " + outgfile)
                    outgcols = outgout.split()
                    outgcount = int(outgcols[0]) / 2

                    # print "Waot " + str(outgcount)

                    origcounts.append(origcount)
                    outgcounts.append(outgcount)





                    with open(outgfile, 'rb') as f_in, gzip.open(outgfile + ".gz", 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)




                statsfile.write(subsubdir)
                statsfile.write(",")
                statsfile.write("MTb")
                statsfile.write(",")
                statsfile.write(str(np.mean(origcounts)))
                statsfile.write(",")
                statsfile.write(str(np.std(origcounts)))
                statsfile.write(",")
                statsfile.write(str(np.mean(outgcounts)))
                statsfile.write(",")
                statsfile.write(str(np.std(outgcounts)))
                statsfile.write("\n")






                if kraken:
                    print "RELEASE THE KRAKEN!!!!"
                    bar = progressbar.ProgressBar()
                    for i in bar(range(len(flist))):
                        ffile = flist[i]

                        bash_command(
                            "kraken --threads " + threads + " --db " + krakendb + " " + ffile + " > " + ffile + ".kraken")

                        bash_command("kraken-translate --db " + krakendb + " " + ffile + ".kraken > " + ffile + ".labels")

                        bash_command(
                            "kraken-translate --db " + krakendb + " " + ffile + ".kraken| cut -f2 | python /data/scripts/make_counts.py > " + ffile + ".kraken4krona")

                        bash_command("ktImportText " + ffile + ".kraken4krona -o " + ffile + ".kraken.krona.html")


                #
                # print "\nCompressing..."
                # bar = progressbar.ProgressBar()
                # for i in bar(range(len(flist))):
                #     ftrz = flist[i]
                #     if os.path.isfile(ftrz) and os.path.isfile(ftrz + ".gz"):
                #         os.remove(ftrz + ".gz")
                #     if os.path.isfile(ftrz):
                #         with open(ftrz, 'rb') as f_in, gzip.open(ftrz + ".gz", 'wb') as f_out:
                #             shutil.copyfileobj(f_in, f_out)
                #         os.remove(ftrz)



                os.chdir("..")

            statsfile.close()
            os.chdir("..")




    ###new style
    else:


        statsfile = open(os.path.dirname(wd) + "nomos-stats.csv", 'w')
        statsfile.write("SubDir, Frag. Length, Ref,Orig Mean Reads per File, Orig StD. Reads, Orig Mean Frag. Length, Orig. Mean Reads Matching Ref. per File,Orig. StD Reads Matching Ref. per File, Orig Mean Matching Ref. Frag. Length, Mapped Mean Reads, Mapped StD Reads, Mapped Mean Frag. Length, Correctly Mapped Mean Reads, Correctly Mapped StD Reads, Correctly Mapped Mean Frag. Length, Incorrectly Mapped Mean Frag. Length\n")

        if multisub:
            subdirs = next(os.walk(wd))[1]
            for subdir in subdirs:
                print "\n*******\nNow working in subdirectory " + subdir
                os.chdir(wd + "/" + subdir)

                nomosnew(statsfile, subdir)

        else:
            nomosnew(statsfile, "N/A")



        statsfile.close()












    cmdfile.close()
    logfile.close()
    exit(0)



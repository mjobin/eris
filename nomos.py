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
    parser.add_argument('-outg', metavar='<outg>', help='Name of genome.',
                        default='')
    parser.add_argument('-dam', dest='dam', help='Search for files with the .dam. extension. Otherwise search for .wodam.',
                        action='store_true')
    parser.set_defaults(dam=False)
    parser.add_argument('-oldstyle', dest='oldstyle', help='Eris was run under old style',
                        action='store_true')
    parser.set_defaults(oldstyle=False)




    args = parser.parse_args()
    wd = args.wd
    krakendb = args.krakendb
    blastdir = args.blastdir
    scriptsdir = args.scriptsdir
    threads = args.threads
    verbose = bool(args.verbose)
    outg = args.outg
    dam = args.dam
    oldstyle = bool(args.oldstyle)

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



                print "\nCompressing..."
                bar = progressbar.ProgressBar()
                for i in bar(range(len(flist))):
                    ftrz = flist[i]
                    if os.path.isfile(ftrz) and os.path.isfile(ftrz + ".gz"):
                        os.remove(ftrz + ".gz")
                    if os.path.isfile(ftrz):
                        with open(ftrz, 'rb') as f_in, gzip.open(ftrz + ".gz", 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                        os.remove(ftrz)



                os.chdir("..")

            statsfile.close()
            os.chdir("..")












    ###new style
    else:
        subdirs = next(os.walk(wd))[1]

        for subdir in subdirs:
            print "Now working in " + subdir
            os.chdir(wd + "/" + subdir)


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


            statsfile = open(subdir+"-stats.csv", 'w')
            statsfile.write("Ref,Orig. Mean.,Orig. StD., Ref Mean, Ref StD\n")

            refdirs = next(os.walk(curdir))[1]

            for refdir in refdirs:
                countfile = open(refdir+"/"+refdir+"counts.csv", 'w')
                countfile.write("Orig file, Orig Count, Ref File, Ref Count\n")
                origcounts = []
                outgcounts = []
                for origfile in origfiles:
                    origbase, origext = os.path.splitext(origfile)

                    outgfile = refdir + "/" + origbase + "-" + refdir + ".aln.fasta"
                    if os.path.isfile(outgfile):
                        origout = bash_command("wc -l " + origfile)
                        origcols = origout.split()
                        origcount = int(origcols[0]) / 2

                        outgout = bash_command("wc -l " + outgfile)
                        outgcols = outgout.split()
                        outgcount = int(outgcols[0]) / 2

                        origcounts.append(origcount)
                        outgcounts.append(outgcount)

                        countfile.write(origfile)
                        countfile.write(",")
                        countfile.write(str(origcount))
                        countfile.write(",")
                        countfile.write(outgfile)
                        countfile.write(",")
                        countfile.write(str(outgcount))
                        countfile.write("\n")
                countfile.close()

                statsfile.write(refdir)
                statsfile.write(",")
                statsfile.write(str(np.mean(origcounts)))
                statsfile.write(",")
                statsfile.write(str(np.std(origcounts)))
                statsfile.write(",")
                statsfile.write(str(np.mean(outgcounts)))
                statsfile.write(",")
                statsfile.write(str(np.std(outgcounts)))
                statsfile.write("\n")

            statsfile.close()

            print "RELEASE THE KRAKEN!!!!"
            bar = progressbar.ProgressBar()
            for i in bar(range(len(origfiles))):
                origfile = origfiles[i]

                bash_command("kraken --threads " + threads + " --db " + krakendb + " " + origfile + " > " + origfile + ".kraken")

                bash_command("kraken-translate --db " + krakendb + " " + origfile + ".kraken > " + origfile + ".labels")

                bash_command("kraken-translate --db " + krakendb + " " + origfile + ".kraken| cut -f2 | python /data/scripts/make_counts.py > " + origfile + ".kraken4krona")

                bash_command("ktImportText " + origfile + ".kraken4krona -o " + origfile + ".kraken.krona.html")


            os.chdir("..")





    cmdfile.close()
    logfile.close()
    exit(0)



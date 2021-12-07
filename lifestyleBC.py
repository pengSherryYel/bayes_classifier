#!/usr/bin/env python3
# coding: utf-8
# author: sherry peng
# mail: xue.peng@helmholtz-muenchen.de
# date: 2021.12.6


import re
from collections import defaultdict
import os
from subprocess import Popen
import math
from utility import mkdirs, checkEnv


# ------ function ------
# def runProdigal(inputseq,prefix,wd):
#    cmd = "prodigal -a {2}/{1}.prodigal.gene.faa -d {2}/{1}.prodigal.gene.fna -g 11  -i {0} -o {2}/{1}.prodigal.output -s {2}/{1}.prodigal.gene.score".format(inputseq,prefix,wd)
#    print("RUN command: %s\n"%cmd)
#    obj = Popen(cmd,shell=True)
#    obj.wait()
#    print("prodigal done!")
#    return "%s/%s.prodigal.gene.faa"%(wd,prefix)


def runHmmsearch(inputfile, prefix, wd, hmmModel, otherPara="--cpu 1"):
    '''
    Aim: run hmmer search for a pfam hmm database

    Usage: runHmmsearch(inputfile,prefix,wd,hmmModel,otherPara="--cpu 1")
        inputfile: a protein set from a metabin or a genome 
        prefix: sample ientifier, will be used as a prefix to the output.
        wd: work path where put the result
        hmmModel: a pfam *.hmm file
        otherPara: parameter for run the hmmer search.
            default: "--cpu 1"

    Return: output file path (*.tblout)
    '''

    checkEnv("hmmsearch")
    mkdirs(wd)
    cmd = "hmmsearch {4} -o {2}/{1}.hmmsearch.out --tblout {2}/{1}.hmmsearch.tblout {3} {0}".format(
        inputfile, prefix, wd, hmmModel, otherPara)
    print("RUN command: %s\n" % cmd)
    obj = Popen(cmd, shell=True)
    obj.wait()
    print("hmmsearch done!")
    return "%s/%s.hmmsearch.tblout" % (wd, prefix)


def runMmseqsEasysearch(inputfile, prefix, wd, hmmDB, otherPara="-s 7 --max-seqs 1 --alignment-mode 3 --alignment-output-mode 0 --min-aln-len 40  --cov-mode 0 --greedy-best-hits 1 --threads 30"):
    '''
    Aim: run mmseq easy search for a fasta database

    Usage: runMmseqsEasysearch(inputfile,prefix,wd,hmmDB,otherPara)
        inputfile: a protein set from a metabin or a genome 
        prefix: sample ientifier, will be used as a prefix to the output.
        wd: work path where put the result
        hmmDB: a fasta database built from mmseqs creatdb command
        otherPara: parameter for run the mmseqs easy search.
            default: "-s 7 --max-seqs 1 --alignment-mode 3 --alignment-output-mode 0 --min-aln-len 40  
                      --cov-mode 0 --greedy-best-hits 1 --threads 30"

    Return: output file path (*.m8)
    '''
    checkEnv("mmseqs")
    mkdirs(wd)

    tmpdir = "%s_tmp" % prefix
    mkdirs(tmpdir)
    output = os.path.join(wd, "%s.m8" % prefix)

    cmd = "mmseqs easy-search {0} {1} {2} {3} {4} && rm -rf {3} && echo 'mmesqs search done!'".format(
        inputfile, hmmDB, output, tmpdir, otherPara)
    print("RUN command: %s\n" % cmd)
    obj = Popen(cmd, shell=True)
    obj.wait()
    return output


def load_scoreD(score_file):
    '''
    Aim: Parse the score file 
    Return: dict. d[member] = [temperate_score, lytic_score] 
    '''
    d = {}
    with open(score_file) as f:
        for line in f:
            ref, tmperate, virulent, members = line.strip("\n").split("\t")
            for member in members.split(";"):
                d[member] = [tmperate, virulent]
    return d


# def load_mmseq2_cluster_tsv(mmseq2_tsv):
#    d = defaultdict(list)
#    with open(mmseq2_tsv) as f:
#        for line in f:
#            pre,member = line.strip("\n").split("\t")
#            d[pre].append(member)
#    return d


def load_hmmsearch_opt(hmmsearch_opt, creteria=1e-5):
    '''
    Aim: parse the hmmersearch output. this file contain multiple columns. is the output from -tblout parameter

    Return: dict.  d[refname][annoacc] = description
    '''
    print("loading mmsearch output")
    annoD = defaultdict(dict)
    with open(hmmsearch_opt) as f:
        for line in f:
            if not line.startswith("#"):
                t = re.split("\s+", line.strip("\n"))
                target_name, target_accession, query_name, accession, Evalue, score, bias, bst_Evalue, bst_score, bst_bias,\
                    exp, reg, clu, ov, env, dom, rep, inc, *description_of_target = t
                accession = accession.split(".")[0]
                # print(target_name,Evalue,bst_Evalue)
                if float(Evalue) <= float(creteria) and float(bst_Evalue) <= float(creteria):
                    annoD[target_name][accession] = query_name
        return annoD


def calcaulate_score(mmseqOpt, scoreD, creteria=1e-5):
    '''
    Aim: calculate P(temperate|GC1,GC2,...,GCN) and P(virulent||GC1,GC2,...,GCN) for a given mmseqOpt

    Usage: calcaulate_score(mmseqOpt,scoreD,creteria=1e-5)
        mmseqOpt: mmseq easy seach output with our database from function runMmseqsEasysearch() 
        scoreD: the return dict from function load_scoreD()
        creteria: just select e-value greater than creteria

    Return: p_total_temperate,p_total_lytic,label
    '''
    # load the prior probability
    p_prior_temperate, p_prior_lytic = scoreD.get("Prior_probability", "NA")

    p_temperate = []
    p_lytic = []
    p_temperate.append(p_prior_temperate)
    p_lytic.append(p_prior_lytic)

    # find the smallest e-value result for the query => return a dict
    d = {}
    with open(mmseqOpt) as f:
        for line in f:
            query, ref, iden, length, mismatch, gap, qstart, qend, sstart, send, evalue, bit_score = line.split(
                "\t")
            if query not in d:
                if float(evalue) <= creteria:
                    d[query] = [ref, evalue]
            else:
                if float(evalue) < float(d[query][-1]):
                    d[query] = [ref, evalue]

    # enumerate the dict to add the score of each lifestyle to lists
    for query, values in d.items():
        ref, evalue = values
        p_temperate_gc, p_lytic_gc = scoreD[ref]
        if p_temperate_gc != 0 and p_lytic_gc != 0:
            p_temperate.append(p_temperate_gc)
            p_lytic.append(p_lytic_gc)

    # calculate P(temperate|GC1,GC2,...,GCN) and P(virulent||GC1,GC2,...,GCN)
    label = "NA"
    p_total_temperate, p_total_lytic = 0, 0
    if len(p_temperate) == 1:
        label = "NotEnoughInfo"
    else:
        for t, l in zip(p_temperate, p_lytic):
            p_total_temperate += float(t)
            p_total_lytic += float(l)

    if p_total_temperate > p_total_lytic:
        label = "Temperate"
    elif p_total_temperate < p_total_lytic:
        label = "Virulent"
    else:
        label = "NotEnoughInfo"
    return p_total_temperate, p_total_lytic, label


def chunk_list(inputlist, chunksize=10):
    d = {}
    n = 0
    for i in range(0, len(inputlist), chunksize):
        t = inputlist[i:i+chunksize]
        d[n] = t
        n += 1
    return d


prophage_relate = ['PF00552', 'PF00589', 'PF00665', 'PF02022', 'PF02899',
                   'PF02920', 'PF09003', 'PF12482', 'PF12834', 'PF12835', 'PF13009', 'PF13102',
                   'PF13333', 'PF13495', 'PF13683', 'PF13976', 'PF14659', 'PF14882', 'PF16795',
                   'PF17921', 'PF18103', 'PF18644', 'PF18697', 'PF06806', 'PF07825', 'PF09035']


def bayes_classifier_single(inputfile, prefix, wd, pfam_creteria=1e-5, mmseqs_creteria=1e-5):
    '''
    Aim: single predict lifestyle

    Process: have 2 phase
        phase 1: align to the integrase and excisionase(16 pfam id in total) -- hmmer
        phase 2: align to our protein database -- mmseqs easy-search
        if temperate appear in either one phase --> final will be temperate 

    Usage: bayes_classifier_single(inputfile,prefix,wd,pfam_creteria=1e-5,mmseqs_creteria=1e-5)
        inputfile: protein set from a bin|a genome.
        prefix: sample ientifier, will be used as a prefix to the output.
        wd: work path where put the result
        pfam_creteria: creteria to filter pfam evalue greater than x (default: 1e-5)
        mmseqs_creteria: creteria to filter mmseqs evalue greater than x (default: 1e-5)

    Return:[prefix,inte_label,excision_label,pfam_label,p_total_temperate,p_total_lytic,bc_label,final_label] 
    '''
    mkdirs(wd)
    pfam_label, bc_label, final_label = "Virulent", "Virulent", "Virulent"
    fileDir = os.path.dirname(os.path.abspath(__file__))

    # phase 1 hmmsearch for pfam
    integrase_hmm = os.path.join(fileDir, "db/integrase_pfv34.hmm")
    excisionase_hmm = os.path.join(fileDir, "db/excisionase_pfv34.hmm")

    # run hmmsearch for integrease and excisionase
    pfam_wd = os.path.join(wd, "BC_pfam")
    inte_opt = runHmmsearch(inputfile, "%s.BC_integrase" %
                            prefix, pfam_wd, integrase_hmm)
    excision_opt = runHmmsearch(
        inputfile, "%s.BC_excisionase" % prefix, pfam_wd, excisionase_hmm)
    # print(inte_opt,excision_opt)

    # load pfam result
    inte_label, excision_label = [0, 0]
    inte_annoD = load_hmmsearch_opt(inte_opt, creteria=pfam_creteria)
    excision_annoD = load_hmmsearch_opt(excision_opt, creteria=pfam_creteria)

    if inte_annoD:
        inte_label = len(inte_annoD)
    if excision_annoD:
        excision_label = len(excision_annoD)

    if inte_label or excision_label:
        pfam_label = "Temperate"
    # print(prefix,inte_label,excision_label,pfam_label)

    # phase 2 bayes classifier
    # run mmseqs search
    bc_mmseqsDB = os.path.join(fileDir, "db/bayes_mmseqs_index/all.protein")
    mmseqs_wd = os.path.join(wd, "BC_mmseqs")
    mmseqs_prefix = "%s.BC_mmseqs" % prefix
    mmseq_opt = runMmseqsEasysearch(inputfile, mmseqs_prefix, mmseqs_wd, bc_mmseqsDB,
                                    otherPara="-s 7 --max-seqs 1 --alignment-mode 3 --alignment-output-mode 0 --min-aln-len 40  --cov-mode 0 --greedy-best-hits 1 --threads 30")
    # print(mmseq_opt)

    # load score file
    score_file = os.path.join(fileDir, "db/score.tsv")
    member2scoreD = load_scoreD(score_file)

    # read mmseqs easy search file
    p_total_temperate, p_total_lytic, bc_label = calcaulate_score(
        mmseq_opt, member2scoreD, creteria=mmseqs_creteria)
    print(prefix, inte_label, excision_label, pfam_label,
          p_total_temperate, p_total_lytic, bc_label)

    # phase 3 combine result pfam and bayes classifier
    if pfam_label == "Temperate" or bc_label == "Temperate":
        final_label = "Temperate"
    return [prefix, inte_label, excision_label, pfam_label, p_total_temperate, p_total_lytic, bc_label, final_label]


def bayes_classifier_batch(inputfile, wd, summaryfile="BC_predict.summary", pfam_creteria=1e-5, mmseqs_creteria=1e-5):
    '''
    Aim: batch predict lifestyle
    Usage: bayes_classifier_batch(inputfile,wd,summaryfile,pfam_creteria=1e-5,mmseq_creteria=1e-5)
        inputfile: a tab seperate file contain two column.
            first column: sample name;
            second column: path of the protein file;
        wd: work path where put the result
        summaryfile: file name for summary of the predict output. location will be under the wd path
        pfam_creteria: creteria to filter pfam evalue greater than x (default: 1e-5)
        mmseqs_creteria: creteria to filter mmseqs evalue greater than x (default: 1e-5)
    '''

    mkdirs(wd)
    opt = open(os.path.join(wd, summaryfile), "w")
    header = "sample_name\tintegrase_number\texcisionase_number\tpfam_label\tbc_temperate\tbc_virulent\tbc_label\tfinal_label\tpath\n"
    opt.write(header)

    with open(inputfile) as f:
        for line in f:
            sample_name, path = line.strip("\n").split("\t")
            res = bayes_classifier_single(
                path, sample_name, wd, pfam_creteria, mmseqs_creteria)
            prefix, inte_label, excision_label, pfam_label, p_total_temperate, p_total_lytic, bc_label, final_label = res
            res.append(path)
            opt.write("\t".join([str(i) for i in res])+"\n")
    opt.close()


if __name__ == "__main__":
    # for single predict
    #bayes_classifier_single("./example/simulate_art_sample1.21.faa", "simulate_art_sample1.21", "./test")
    # bayes_classifier_single("./example/simulate_art_sample1.1.faa","simulate_art_sample1.1","./test")
    # bayes_classifier_single("./example/simulate_art_sample1.5.faa","simulate_art_sample1.5","./test")
    # for batch predict
    bayes_classifier_batch("./example/example.list","./batch_test","BC_predict.summary")

#!/usr/bin/env python

import argparse, sys, csv, re
#import numpy as np
from os.path import exists
from collections import defaultdict

def main():

## parse input arguments ##

    parser = argparse.ArgumentParser(description='Calculate predictive value (PV), true +ve rate (TPR), false +ve rate (FPR) and +ve likelihood ratio (LR) of genotypes (GTs) \
        as potential markers for discrete traits. Takes flapjack formatted genotype, phenotype and map files. For a given marker site, the PV (a.k.a precision) of a GT for a particular \
        trait is: PV=n(GT|trait)/total(GT) =TP/(TP+FP) where a positive test (P) represents a known GT call. Can be thought of as the proportion of all samples with a \
        particular GT that actually have the trait. The TPR (=TP/(TP+FN)) can also be considered the flipped PV, i.e. the proportion of all samples with the trait actually positive for the GT \
        (a.k.a. recall or sensitivity). FPR=FP/(FP+TN) which is the proportion of FPs out of all samples without trait. False discovery rate (FDR), false -ve rate (FNR) and true -ve rate (TNR) \
        are reciprocals of PV, TPR and FPR (FDR=1-PV, FNR=1-TPR, TNR=1-FPR). LR=TPR/FPR and is used for assessing the value of a positive GT call in usefully changing the probability that a \
        the trait exists in a test sample. An ideal marker should have both PV and TPR close to 1, with a LR of 10 or more (assuming a sufficient number and diversity of input samples). \
        Genotype file must first be transposed using "transpose.awk". Alternate orders of alleles (i.e. A/T or T/A) are treated as distinct. Nucleotide groupings are also tested where \
        applicable, e.g. "hasC" comprises all GTs containing "C" (C, C/T, T/C) for a given site. CalcGDV Does not adjust for population stratification and/or LD which the user should factor \
        into their marker/sample selection beforehand. Writes to STDOUT.')
    parser.add_argument('-g', '--genos', required=True, type=str, help='Transposed flapjack genotype file in teb-delimited format.')
    parser.add_argument('-s', '--sites', required=False, type=str, help='Comma separated list of SNP/marker IDs to test or text file listing marker IDs in single column.')
    parser.add_argument('-r', '--range', required=False, type=str, help='Alternative to marker list, provide coordinate range(s) in which markers assessed (requires MAP file). \
        Position info is added to output. Input format: -r <MAP file>,<seq ID>:<start>-<end> OR -r <MAP file>,<BED file of multiple ranges>')
    parser.add_argument('-t', '--trait', type=str, required=True, help='Indicate phenotype file and trait column by which to partition samples for GDV calculations. \
        Any number of trait variables are allowed though must be discrete. Input format: -t <PHENO file>,<column header name>,<missing/unknown identifier string (blank if none)>')
    parser.add_argument('-f', '--filter', type=str, required=False, help='Indicate variable in phenotype file to filter samples by. \
        Will base GDVs on only samples with this variable (eg. country). Input format: -f <column header name>,<identifier string>')
    parser.add_argument('-pv', '--pvmin', type=float, required=False, help='Indicate float minimum PV which to report (PV > min per GT-trait).')
    parser.add_argument('-tpr', '--tprmin', type=float, required=False, help='Indicate float minimum TPR which to report (TPR > min per GT-trait).')
    parser.add_argument('-fpr', '--fprmax', type=float, required=False, help='Indicate float maximum FPR which to report (TPR < max per GT-trait).')
    parser.add_argument('-lr', '--lrmin', type=float, required=False, help='Indicate float minimum LR which to report (LR > min per GT-trait).')
    parser.add_argument('-n', '--numgt', type=int, required=False, help='Indicate the minimum n occurences of a GT for any given marker/site for it to be reported.')
    parser.add_argument('--strmiss', type=str, required=False, default='-', help='Indicate alternate missing data string (default="-") used in genotype file. \
        Provide empty string "" if none.')
    parser.add_argument('-v', '--verbose', action='store_true', required=False, help='Print out info lines for all markers, even if they fail minimum thresholds.')
    args = parser.parse_args()

## check required args and setup vars
    # check geno file exists
    if not exists(args.genos):
        sys.exit(f"ERROR: file '{args.genos}' not found!\n")
    # listify trait arg and strip any whitespace
    trait_args0 = args.trait.split(",", 2)
    trait_args1 = [x.lstrip() for x in trait_args0]
    # check pheno file exists
    if not exists(trait_args1[0]):
        sys.exit(f"ERROR: file '{trait_args1[0]}' not found!\n")
    # 1st comment to stdout
    print("## input cmd: " + " ".join(sys.argv))
    # initialise set of marker IDs to test
    marker_set = set()

## process sites argument ##
    # listify elements or if file-of-sites exists, open and retrieve IDs
    if args.sites:
        if "," in args.sites:
            sites_args0 = args.sites.split(",")
        elif not exists(args.sites):
            sys.stderr.write(f"WARNING: file '{args.sites}' not found. Assuming as marker ID string.\n")
            sites_args0 = [args.sites]
        else:
            with open(args.sites) as opt_file1:
                sites_args0 = opt_file1.readlines()
    # strip any whitespace/newlines and update marker set
        sites_args1 = [x.strip() for x in sites_args0]
        marker_set.update(sites_args1)

## process range argument ##
    if args.range:
        ranges_dct = {} #will have k,v as seqID:(start,end)
        coords_dct = {} #will have k,v as markerID:(seqID,pos)
        range_args0 = re.split(',|:|-', args.range)
        range_args1 = [x.lstrip() for x in range_args0]
        # check range args
        if not exists(range_args1[0]): #check map file exists
            sys.exit(f"ERROR: file '{range_args1[0]}' not found!\n")
        if ':' in args.range and len(range_args1) < 4 or len(range_args1) < 2:
            sys.exit(f"ERROR: incorrect format for option --range\n")
        elif ':' in args.range:
            ranges_dct[range_args1[1]]=(int(range_args1[2]),int(range_args1[3])) #assign seqid as k with v as (start,end)
        elif len(range_args1) == 2 and exists(range_args1[1]): #check bed file exists if only 2 args
            #parse bed file
            with open(range_args1[1]) as bed_file:
                reader_bed = csv.reader(bed_file, delimiter='\t', quoting=csv.QUOTE_NONE)
                for row in reader_bed:
                    if not row[0].lstrip().startswith('#'):
                        break #skips comment lines
                for row in reader_bed:
                    try:
                        ranges_dct[row[0]]=(int(row[1]),int(row[2])) #assign seqid as k with v as (start,end)
                    except IndexError:
                        sys.stderr.write(f"WARNING: index error parsing '{range_args1[1]}'.\n")
        else:
            sys.stderr.write(f"WARNING: file '{range_args1[1]}' not found. Assuming as sequence ID.\n")
            ranges_dct[range_args1[1]]=(0,float('inf')) #if no range given and 2nd arg not a file, set as seqid with range 0 to inf

        #parse map file
        with open(range_args1[0], 'r') as opt_file2:
            reader_map = csv.reader(opt_file2, delimiter='\t', quoting=csv.QUOTE_NONE)
            n_line = 0
            for row in reader_map:
                n_line += 1
                if not row[0].lstrip().startswith('#'):
                    break #skips comment lines
            sites_in_range = []
            for row in reader_map:
                n_line += 1
                try:
                    marker = row[0]
                    seqid = row[1]
                    pos = int(row[2])
                    if seqid in ranges_dct and pos >= ranges_dct[seqid][0] and pos <= ranges_dct[seqid][1]:
                        sites_in_range.append(marker)
                        coords_dct[marker]=(seqid,pos) #add position info to coords dict
                except (TypeError, ValueError, IndexError) as error:
                    sys.stderr.write(f"ERROR: '{range_args1[0]}' line {n_line} ({error})\n")
                    continue
        # update marker set
        marker_set.update(sites_in_range)

## handle sites and range args exception, check pvmin, tprmin, lrmin args ##
    if args.sites or args.range:
        n_sites = len(marker_set)
        if n_sites < 1:
            sys.exit(f"ERROR: 0 sites added to list! Ensure correct inputs of --sites or --range. Exiting...\n")
        else:
            print(f"## {n_sites} markers added to list to test.")
    else:
        marker_set = False
        sys.stderr.write(f"WARNING: marker set not specified - all markers in '{args.genos}' will be assessed.\n")
    if any([args.pvmin, args.tprmin, args.fprmax, args.lrmin, args.numgt]):
        print_gt_lst = True #list of GTs for each marker added to header if cutoff option or verbose is used
        if any(x < 0.0 or x > 1.0 for x in (x for x in (args.pvmin, args.tprmin, args.fprmax)) if x is not None):
            sys.exit(f"ERROR: invalid input for cutoff(s) - pvmin|tprmin|fprmax must be between 0 and 1. Exiting...\n")
    else:
        print_gt_lst = False

## process trait argument ##
    # assign and check col name and identifier strs
    try:
        colname_trait = trait_args1[1]
    except IndexError:
        sys.exit(f"ERROR: --trait option missing required arguments!\n")
    if len(trait_args1) == 3:
        miss_trait_id = trait_args1[2]
    else:
        miss_trait_id = "" #if no value for missing var given, assign as empty string
        sys.stderr.write(f"WARNING: identifier for missing/unknown trait not given - all non-empty vars considered as traits.\n")
    # open pheno file and retrieve header
    with open(trait_args1[0], 'r') as req_file1:
        reader_pheno = csv.reader(req_file1, delimiter='\t', quoting=csv.QUOTE_NONE)
        n_line = 0
        for row in reader_pheno:
            n_line += 1
            if not row[0].lstrip().startswith('#'):
                break #skips comment lines
        header_pheno = [x for x in row if x.strip()] #remove any empty or whitespace-only elements
    # retrieve index of col name
        try:
            colidx_trait = header_pheno.index(colname_trait) + 1 #adding one skips the sample column which does not have a header name
        except ValueError:
            sys.exit(f"ERROR: '{colname_trait}' was not found in the pheno file header!\n")

## process filter argument ##   
    # listify elements and strip any whitespace
        if args.filter:
            filt_args0 = args.filter.split(",", 1)
            filt_args1 = [x.lstrip() for x in filt_args0]
    # assign and check col name and identifier strs
            colname_filt = filt_args1[0]
            if len(filt_args1) == 2:
                filt_var_id = filt_args1[1]
    # retrieve index of col name
                try:
                    colidx_filt = header_pheno.index(colname_filt) + 1 #adding one skips the sample column which does not have a header name
                except ValueError:
                    sys.exit(f"ERROR: '{colname_filt}' was not found in the pheno file header!\n")
            else:
                args.filter = False
                sys.stderr.write(f"WARNING: identifier for --filter not given - all samples are considered.\n")

## parse whole pheno file ##
        n_samples = 0
        trait_samplesets = defaultdict(set)
        for row in reader_pheno:
            n_line += 1
            try:
                if args.filter and row[colidx_filt] != filt_var_id:
                    continue
                if row[colidx_trait] == miss_trait_id:
                    continue
                else:
                    sample = row[0]
                    trait = row[colidx_trait]
                    trait_samplesets[trait].add(sample) #each unique trait becomes a key for the set of corresponding samples (added iteratively according to their trait value)
                    n_samples += 1
            except (TypeError, ValueError, IndexError) as error:
                sys.stderr.write(f"ERROR: '{trait_args1[0]}' line {n_line} ({error})\n")
        if n_samples < 2:
            sys.exit(f"ERROR: number of samples added ({n_samples}) too low! Ensure arguments match target strings in '{trait_args1[0]}'.\n")
        print(f"## {n_samples} samples added from '{trait_args1[0]}' after filtering.\n## ", end="")
        for trait in trait_samplesets.keys():
            print(trait, end=", ")
        print(f"traits added from '{trait_args1[0]}' with no. of samples as ", end="")
        for trait in trait_samplesets:
            print(len(trait_samplesets[trait]), end=", ")
        print("respectively.")

## parse genotype data ##
    # open geno file and retrieve header
    with open(args.genos, 'r') as req_file2:
        strip_lines = (line.lstrip() for line in req_file2) #remove leading whitespace after transposition from original fj format
        reader_geno = csv.reader(strip_lines, delimiter='\t', quoting=csv.QUOTE_NONE)
        header_geno = next(reader_geno)
    # map samples in samplesets to column index in geno file using header list
        sampleset_indices = defaultdict(list)
        for trait in trait_samplesets.keys():
            for sample in trait_samplesets[trait]:
                try:
                    colidx_samp = header_geno.index(sample) + 1 #adding one skips the marker id column which does not have a header name
                    sampleset_indices[trait].append(colidx_samp)
                except ValueError:
                    sys.stderr.write(f"WARNING: no match for '{sample}' in '{args.genos}' header! Skipping...\n")
                    continue
        pooled_indices = sum(sampleset_indices.values(), []) # create list of all sample indices pooled
        samples_avail = len(pooled_indices)
        if samples_avail < 2:
            sys.exit(f"ERROR: <2 samples matched to '{args.genos}' header! Ensure sample names in '{trait_args1[0]}' and '{args.genos}' match and latter is correctly formatted (samples as colnames, markers as rownames, no comment lines).\n")
        print(f"## {samples_avail}/{n_samples} samples found in '{args.genos}'.")
    # sort indices so they are in ascending order (I know not necessary since idx lookup is O(1) anyway but I just like order ok)
        samplesets_avail = dict()
        for trait in sampleset_indices:
            sampleset_indices[trait].sort()
            samplesets_avail[trait] = len(sampleset_indices[trait])
        pooled_indices.sort()
        traits = sorted(list(sampleset_indices.keys())) #list of trait names
        header_out = '\n#<ID>\t<GT>\t<n>\t<trait>\t<PV>\t<TPR>\t<FPR>\t<LR>' #header str for later
        if args.range:
            header_out += '\t<seq>\t<pos>' #add position fields if range(s) given

## iterate over marker genotypes (GTs) ##
        nucleos = ['A', 'T', 'C', 'G']
        n_line = 1 # nth line starts at 1 to include header line
        if marker_set:
            num_markers = len(marker_set)
        else:
            num_markers = -1
        seen_count = 0
        out_count = 0
        for row in reader_geno:
            n_line += 1
            out_str = '' #multiline output str to add to - will only be printed at final step if at least 1 GT meets criteria (hit_count > 0)
            try:
                if marker_set == False or row[0] in marker_set:
                    marker = row[0]
                    seen_count += 1
                    # initialize master dict of trait keys whose values are inner dicts of GT:PV key:value pairs
                    row_master_dct = dict() 
                    # collate GTs of all eligible samples for marker (row)
                    GT_all_lst = [row[i] for i in pooled_indices if row[i] != args.strmiss]
                    GT_all_set = set(GT_all_lst)
                    GT_all_dct = dict.fromkeys(GT_all_set, 0) #each unique element assigned starting value 0
                    row_avail = len(GT_all_lst)
                    row_missing = samples_avail - row_avail
                    #add 1st info comment for marker to out string
                    out_str += f"#markerID={marker},"
                    if args.range:
                        try:
                            out_str += f"pos={coords_dct[marker][0]}_{coords_dct[marker][1]}," #add position info if range(s) given
                        except KeyError:
                            out_str += f"pos=n/a,"
                    out_str += f"GTavail={row_avail},GTmissing={row_missing}" 
                    # collect counts of GTs across all eligible samples (excl missing)
                    for gt in GT_all_set:
                        GT_all_dct[gt] = GT_all_lst.count(gt)
                    # sum counts of conditional GTs (e.g. has 'A' = A|N/A|A/N) across all eligible samples and add to dict
                    conGT_all = {} # initialise seperate dict for conditional GT counts (prevents double-counting cond GTs from previous iterations if were to add directly to GT_sub_dct)
                    for nt in nucleos:
                        has_nt = 0
                        non_nt = 0
                        hmatch = 0 #n matching GTs
                        nmatch = 0 #n non-matching GTs
                        for gt, count in GT_all_dct.items():
                            if nt in gt:
                                hmatch += 1
                                has_nt += count
                            else:
                                nmatch += 1
                                non_nt += count
                        if hmatch > 1 and has_nt > 0:
                            conGT_all['has' + nt] = has_nt # only report 'has' if there are >1 matching GTs and nt occurs at least once
                            if nmatch > 1:
                                conGT_all['non' + nt] = non_nt # only report 'non' if there are >1 non-match GTs
                    GT_all_dct.update(conGT_all) # merge GT and conditional GT dicts
                    # collate GTs of each sample subset (trait) for marker (row)
                    for trait in traits:
                        row_master_dct[trait] = {} # initialise inner dict for trait
                        GT_sub_lst = [row[i] for i in sampleset_indices[trait] if row[i] != args.strmiss]
                        GT_sub_set = set(GT_sub_lst)
                        GT_sub_dct = dict.fromkeys(GT_sub_set, 0) #each unique element assigned starting value 0
                        sub_avail = len(GT_sub_lst)
                        sub_missing = samplesets_avail[trait] - sub_avail
                        out_str += f"\n#trait={trait},GTavail={sub_avail},GTmissing={sub_missing}" #add per-trait info comment to out str
                    # collect counts of GTs across subset samples (excl missing)
                        for gt in GT_sub_set:
                            GT_sub_dct[gt] = GT_sub_lst.count(gt) # equal to TP
                    # sum counts of conditional GTs across subset samples and add to dict
                        conGT_sub = {} # initialise seperate dict for conditional GT counts (prevents double-counting cond GTs from previous iterations if were to add directly to GT_sub_dct)
                        for nt in nucleos:
                            has_nt = 0
                            non_nt = 0
                            for gt, count in GT_sub_dct.items():
                                if nt in gt:
                                    has_nt += count
                                else:
                                    non_nt += count
                            conGT_sub['non' + nt] = non_nt
                            conGT_sub['has' + nt] = has_nt
                        GT_sub_dct.update(conGT_sub) # merge GT and conditional GT dicts
                    # calculate stats, add GT:[PV, TPR, FPR, LR] list to inner trait dict of outer master dict
                        for gt in GT_sub_dct.keys():
                            if gt in GT_all_dct:
                                try:
                                    #TP = GT_all_dct[gt]
                                    PV = GT_sub_dct[gt] / GT_all_dct[gt] #proportion of samples with gt actually having the trait
                                    FP = GT_all_dct[gt] - GT_sub_dct[gt] #samples with gt that do not have trait
                                    FN = sub_avail - GT_sub_dct[gt] #samples with trait that do not have gt
                                    TN = row_avail - (GT_all_dct[gt] + FN) #samples with neither trait nor gt
                                    TPR = GT_sub_dct[gt] / sub_avail #proportion of TPs out of all samples with trait
                                    FPR = FP / (FP + TN) #proportion of FPs out of all samples without trait
                                    #TNR = 1 - FPR
                                    #FNR = 1 - TPR
                                    #FDR = 1 - PV
                                    if FPR > 0.0:
                                        LR = TPR / FPR
                                    else:
                                        LR = float('nan')
                                except ZeroDivisionError:
                                    continue
                                row_master_dct[trait][gt] = [round(PV, 3), round(TPR, 3), round(FPR, 3), round(LR, 3)]
                    # extract all GT keys and sort first alphanumericlly and then by str length - ensures consistency in order GT results printed
                    GT_key_lst = list(GT_all_dct.keys())
                    GT_key_lst.sort(key=str)
                    GT_key_lst.sort(key=len)
                    if print_gt_lst or args.verbose:
                        out_str += f"\n#GTs={','.join(GT_key_lst)}"
                    # print results for marker (row) to stdout if GTs pass cutoffs
                    out_str += header_out #add header to out str
                    hit_count = 0
                    for gt in GT_key_lst:
                        if args.numgt and GT_all_dct[gt] < args.numgt:
                            continue
                        for trait in traits:
                            try:
                                stats_out=row_master_dct[trait][gt] #for GT, retrieve stats for each trait
                            except KeyError:
                                continue
                            if args.pvmin and stats_out[0] < args.pvmin:
                                continue #skip if PV < user defined min
                            if args.tprmin and stats_out[1] < args.tprmin:
                                continue #skip if TPRs < user defined min
                            if args.fprmax and stats_out[2] > args.fprmax:
                                continue #skip if FPRs >= user defined max
                            if args.lrmin and stats_out[3] < args.lrmin:
                                continue #skip if LR < user defined min
                            out_str += f"\n{marker}\t{gt}\t{GT_all_dct[gt]}\t{trait}\t" #add GT result cols to out row
                            out_str += '\t'.join([str(x) for x in stats_out]) #add trait stats result cols to out row
                            if args.range:
                                try:
                                    out_str += f"\t{coords_dct[marker][0]}\t{coords_dct[marker][1]}" #add position info to out row if range(s) given
                                except KeyError:
                                    out_str += "\tn/a\tnan"
                            hit_count += 1
                    if hit_count > 0 or args.verbose:
                        out_count += hit_count
                        print(out_str) #print final output for row if conditions met or verbose==True
            except (ValueError, IndexError) as error:
                sys.stderr.write(f"ERROR: '{args.genos}' line {n_line} ({error})\n")
                continue
            # stop iterating through rows if all markers found and print final counts
            if seen_count == num_markers and num_markers > -1:
                print(f"## end: {seen_count}/{num_markers} markers found up to line {n_line} of '{args.genos}'. Results for {out_count} genotypes reported.")
                sys.stderr.write(f"DONE: {seen_count}/{num_markers} markers found up to line {n_line} of '{args.genos}'. GDVs for {out_count} genotypes reported.\n")
                break
        # when end of file reached before seen == num markers or marker not given, print final counts 
        if seen_count < num_markers:
            print(f"## end: {seen_count}/{num_markers} markers found up to last line ({n_line}) of '{args.genos}'. Results for {out_count} genotypes reported.")
            sys.stderr.write(f"DONE: {seen_count}/{num_markers} markers found up to last line ({n_line}) of '{args.genos}'. GDVs for {out_count} genotypes reported.\n")
        if num_markers == -1:
            print(f"## end: {seen_count} markers found up to last line ({n_line}) of '{args.genos}'. Results for {out_count} genotypes reported.")
            sys.stderr.write(f"DONE: {seen_count} markers found up to last line ({n_line}) of '{args.genos}'. GDVs for {out_count} genotypes reported.\n")          

if __name__ == '__main__':
    main()

#per marker results output should be like:
#markerID=<row[0]>,GTavail=<row_avail>,GTmissing=<row_missing>
#trait=<trait1>,GTavail=<sub_avail>,GTmissing=<sub_missing>
#...
#<ID>  <GT>    <n>  <trait> <PV>   <TPR> <FPR>   <LR>   <seq>   <pos>
#<row[0]> <gt>  <GT_all_dct[gt]> <trait1>    <row_master_dct[trait1][gt][0]>    <row_master_dct[trait1][gt][1]> <row_master_dct[trait1][gt][2]> <row_master_dct[trait1][gt][3]> <coords_dct[row[0]][0]> <coords_dct[row[1]][1]>
#<row[0]> <gt>  <GT_all_dct[gt]> <trait2>    <row_master_dct[trait2][gt][0]>    <row_master_dct[trait2][gt][1]> <row_master_dct[trait2][gt][2]> <row_master_dct[trait2][gt][3]> <coords_dct[row[0]][0]> <coords_dct[row[1]][1]>
#...
#i.e. info lines, header line, GTs as rows, stats for each trait as extra rows of same GT
#seq and pos fields only added if range args given - seq and pos values sourced from map file provided with range(s)


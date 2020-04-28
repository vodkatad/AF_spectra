#!/usr/bin/env python3

import argparse
import sys

def next_bed_entry(bedcn):
    line = bedcn.readline()
    if line != '':
        line.rstrip('\n')
        entry = line.split('\t')
        return([entry[0], int(entry[1]), int(entry[2]), float(entry[3])])
    return(None)

# return the first entry that overlap current and True if the bed is finished
# schiavismo di fidanzate
last_entry = None
already_given = False
def get_next_overlapping_bed(current, bedcn, verbose):
    global last_entry
    global already_given
    if last_entry == None:
        last_entry = next_bed_entry(bedcn)
        if last_entry == None:
            return(None, True)

    result = None 
    if last_entry != None and last_entry[0] == current[0] and last_entry[1] < current[2] and current[1] < last_entry[2]:
        if verbose:
            print('evaluating overlap {} {} {}'.format(last_entry[0], last_entry[1], last_entry[2]), file=sys.stderr)
        result = last_entry
        if last_entry[2] <= current[2]:
            last_entry = None
        else:
            if already_given:
                result = None
            already_given = not already_given

    return(result, False)
        
def get_next_bin(current, binlen, chrs):
     # we have space for another bin
    if current[2] + binlen < chrs[current[0]]: # check ends here TODO
        return (current[0], current[2], current[2] + binlen)
    # we need to get the next chr
    elif current[2] == chrs[current[0]]:
        k = list(chrs.keys()) # ordered?
        chri = k.index(current[0])
        if chri + 1 < len(k):
            return (k[chri+1], 0, binlen) 
        else:
            return None
    # we need the last piece of this chr
    else:
        return(current[0], current[2], chrs[current[0]])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Project bed file with cn on a fixed size binning")

    parser.add_argument('--bedcn', '-c', dest='bedcn', action='store', help='a bed file with cn in the name [4th] field. Has to be sorted!', required=True)
    parser.add_argument('--bin', '-b', dest='bin', action='store', type=int, default=10000, help='size in bp of the bin where you want to project the cn', required=False) 
    #parser.add_argument('--bin', '-b', dest='bin', action='store', type=int, default=10000, help='size in bp of the bin where you want to project the cn', required=False) 
    #parser.add_argument('--chrlen', '-l', dest='chrlen', action='store', help='file with size of chrs', required=False) 
    # TODO add possibility to be flexible with chr and chr lengths
    parser.add_argument('--verbose', '-v', action='store_true',  help='verbose execution')

    #CHRS = ['chr' + x for x in range(1,22)]

    #[egrassi@occam matched_normals_30x_downsampled]>head -n 22  ~/bit/task/annotations/dataset/gnomad/GRCh38.d1.vd1.allchr.bed | cut -f 1,2 | awk '{print "@",$1,"@",":",$2}' | tr "@" '"' | tr -d " "
    chrlen = {
        "chr1": 248956422,
        "chr2": 242193529,
        "chr3": 198295559,
        "chr4": 190214555,
        "chr5": 181538259,
        "chr6": 170805979,
        "chr7": 159345973,
        "chr8": 145138636,
        "chr9": 138394717,
        "chr10": 133797422,
        "chr11": 135086622,
        "chr12": 133275309,
        "chr13": 114364328,
        "chr14": 107043718,
        "chr15": 101991189,
        "chr16": 90338345,
        "chr17": 83257441,
        "chr18": 80373285,
        "chr19": 58617616,
        "chr20": 64444167,
        "chr21": 46709983,
        "chr22": 50818468
    }

    args = parser.parse_args()
    if (args.verbose):
        print('cn_bed\t{}'.format(args.bedcn), file=sys.stderr)
        print('bin_size\t{}'.format(args.bin), file=sys.stderr)

    
    # Tuple with current bin: chr, b, e. 0 based end excluded
    next_bin = None
    done = False
    current = None
    with open(args.bedcn, 'r') as bedcn:
        #entry = next_bed_entry(bedcn) # ???
        while not done:
            # First chr to be considered from the bed directly
            if current == None:
                current = ('chr1', 0, args.bin)
            else:
                # We should always have completely consumed the previous bin in our previous loop
                current = get_next_bin(current, args.bin, chrlen)
                # End of bins
                if current == None:
                    done = True
                    break
            overlap = []
            overlaplen = []
            # overlap check  a0 <= b1 && b0 <= a1;
            # https://fgiesen.wordpress.com/2011/10/16/checking-for-interval-overlap/
            # < and not <= for end excluded
            while True:
                entry, done = get_next_overlapping_bed(current, bedcn, args.verbose)
                if entry is None:
                    break
                overlap.append(entry)
                overlaplen.append(min(entry[2], current[2]) - max(entry[1], current[1]))

            # manage the overlapping entries. We get a weighted average of their cn,
            # where the weight is the overlap length
            if len(overlap) != 0:
                cn = 0
                ovlen = 0
                for i in range(0, len(overlap)):
                    cn += overlap[i][3] * overlaplen[i]
                    ovlen += overlaplen[i]
                cn =  cn / ovlen
                print('{}\t{}\t{}\t{}'.format(current[0],current[1],current[2], cn))
            else:
                print('{}\t{}\t{}\t{}'.format(current[0],current[1],current[2], 0))
            
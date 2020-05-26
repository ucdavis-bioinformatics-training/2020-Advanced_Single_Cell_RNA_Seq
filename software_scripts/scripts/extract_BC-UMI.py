#!/usr/bin/env python
'''
'extract' the first 'length'  basepairs
(presumably the cell barcode and UMI sequence)
from 'read' and places the sequence and quality
as a tag in the read id.

'insert' place the extracted length from reads
 back on the read.

input/output is hts_stream tab file.

WARNING: This script has not been fully tested!!!!
'''

import sys
import os
from optparse import OptionParser  # http://docs.python.org/library/optparse.html
import gzip


usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option('--insert', help="insert extracted tags into library",
                  action="store_true", dest="insert",default=False)
parser.add_option('--extract', help="extract to tags",
                  action="store_true", dest="extract", default=False)
parser.add_option('--read', help="the read the barcode and umi can be found on",
                  action="store", type="int", dest="read")
parser.add_option('--length', help="the total length of the BC and UMI",
                  action="store", type="int", dest="length")

(options,  args) = parser.parse_args()  # uncomment this line for command line support

intab = sys.stdin
outtab = sys.stdout

if options.insert and options.extract:
    parser.error("options --insert and --extract are mutually exclusive")

if options.extract and not options.length:
    parser.error("when extracting need to specify length")

read = options.read
length = options.length
extract = options.extract
insert = options.insert

def reverseComplement(s):
    """
    given a seqeucne with 'A', 'C', 'T', and 'G' return the reverse complement
    """
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    letters = list(s)
    try:
        letters = [basecomplement[base] for base in letters]
    except:
        raise
    return ''.join(letters[::-1])

def reverse(s):
    """
    given a sequence return the reverse
    """
    letters = list(s)
    return ''.join(letters[::-1])

for line in intab:
    #Comment/header lines start with @
    line = line.strip('\n').split('\t')
    if len(line) > 2 and len(line) < 9:
        #Handle SE:
        if (len(line) == 3 or len(line) == 4):
            tags = []
            tags2 = ''
            ID1 = line[0]
            Seq1 = line[1]
            Qual1 = line[2]
            if (len(line) == 4):
                tags.extend(line[3].split('~'))
                if len(tags) > 1:
                    tags2 = tags[1]
                tags = tags[0].split("|")

            if (extract):
                if len(Seq1) < length:
                    continue
                # TR:Z:SEQ
                # TY:Z:Qual
                if (read == 1):
                    # If barcode|UMI is on Read 1 then extract from first part of read
                    tags.append("TR:Z:" + Seq1[:length])
                    tags.append("TY:Z:" + Qual1[:length])
                    Seq1 = Seq1[length:]
                    Qual1 = Qual1[length:]
                elif (read == 2):
                    # If barcode|UNI is on Read 2 the extract from end of the read
                    tags.append("TR:Z:" + reverseComplement(Seq1[-length:]))
                    tags.append("TY:Z:" + reverse(Qual1[-length:]))
                    Seq1 = Seq1[:(len(Seq1)-length)]
                    Qual1 = Qual1[:(len(Qual1)-length)]
                else:
                    sys.exit("read is not a 1 or 2")

                outtab.write(ID1 + '\t')
                outtab.write(Seq1 + '\t')
                outtab.write(Qual1 + '\t')
                outtab.write('|'.join(tags))
                if tags2 != '':
                    outtab.write('~' + tags2)
                outtab.write('\n')
                continue
            if (insert):
                tr = [i for i in tags if "TR:Z:" in i]
                if len(tr) != 1:
                    sys.exit("Didn't find the TR tag")
                ty = [i for i in tags if "TY:Z:" in i]
                tags.remove(tr[0])
                tags.remove(ty[0])

                final_tags = ''
                if len(tags) > 0:
                    final_tags = '|'.join(tags) + '\t'
                if tags2 != '':
                    final_tags += tags2

                if (read == 1):
                    outtab.write(ID1 + '\t')
                    outtab.write(tr[0][5:] + '\t')
                    outtab.write(ty[0][5:] + '\t')
                    outtab.write(ID1.replace(" 1:", " 2:") + '\t')
                    outtab.write(reverseComplement(Seq1) + '\t')
                    outtab.write(reverse(Qual1))
                elif (read == 2):
                    outtab.write(ID1 + '\t')
                    outtab.write(Seq1 + '\t')
                    outtab.write(Qual1 + '\t')
                    outtab.write(ID1.replace(" 1:", " 2:") + '\t')
                    outtab.write(tr[0][5:] + '\t')
                    outtab.write(ty[0][5:])
                else:
                    sys.exit("read is not a 1 or 2")

                if final_tags != '':
                    outtab.write('\t' + final_tags)
                outtab.write('\n')
                continue
        #Handle PE:
        if (len(line) == 6 or len(line) == 8):
            tags1 = []
            tags2 = []
            ID1 = line[0]
            # TODO, if Seq1 or Seq2 are 'N' (no seq) then don't append
            Seq1 = line[1]
            Qual1 = line[2]
            ID2 = line[3]
            Seq2 = line[4]
            Qual2 = line[5]
            if (len(line) == 8):
                tags1.extend(line[6].split('|'))
                tags2.extend(line[7].split('|'))

            if (extract):
                if read == 1 and len(Seq1) < length:
                    continue
                elif read == 2 and len(Seq2) < length:
                    continue
                # TR:Z:SEQ
                # TY:Z:Qual
                outtab.write(ID1 + '\t')
                if (read == 1):
                    # If barcode|UMI is on Read 1 then extract from first part of read
                    tags1.append("TR:Z:" + Seq1[:length])
                    tags1.append("TY:Z:" + Qual1[:length])
                    Seq1 = Seq1[length:]
                    Qual1 = Qual1[length:]
                    if (len(Seq1) > 0):
                        outtab.write(Seq1 + '\t')
                        outtab.write(Qual1 + '\t')
                        outtab.write(ID2 + '\t')
                        outtab.write(Seq2 + '\t')
                        outtab.write(Qual2 + '\t')
                    else:
                        outtab.write(reverseComplement(Seq2) + '\t')
                        outtab.write(reverse(Qual2) + '\t')

                    outtab.write('|'.join(tags1))
                    outtab.write('\t' if len(Seq1) > 0 else '~')
                    outtab.write('|'.join(tags2))
                    outtab.write('\n')
                elif (read == 2):
                    # If barcode|UNI is on Read 2 the extract from end of the read
                    tags1.append("TR:Z:" + Seq2[:length])
                    tags1.append("TY:Z:" + Qual2[:length])
                    Seq2 = Seq2[length:]
                    Qual2 = Qual2[length:]
                    outtab.write(Seq1 + '\t')
                    outtab.write(Qual1 + '\t')
                    if (len(Seq2) > 0):
                        outtab.write(ID2 + '\t')
                        outtab.write(Seq2 + '\t')
                        outtab.write(Qual2 + '\t')
                    outtab.write('|'.join(tags1))
                    outtab.write('\t' if len(Seq2) > 0 else '~')
                    outtab.write('|'.join(tags2))
                    outtab.write('\n')
                else:
                    sys.exit("read is not a 1 or 2")
                continue
            if (insert):
                tr = [i for i in tags1 if "TR:Z:" in i]
                if len(tr) != 1:
                    sys.exit("Didn't find the TR tag")
                ty = [i for i in tags1 if "TY:Z:" in i]
                tags1.remove(tr[0])
                tags1.remove(ty[0])
                final_tags = ''
                if len(tags1) > 0:
                    final_tags = '|'.join(tags1) + '\t'
                if tags2 != '':
                    final_tags += '|'.join(tags2)

                if (read == 1):
                    outtab.write(ID1 + '\t')
                    outtab.write(tr[0][5:] + Seq1 + '\t')
                    outtab.write(ty[0][5:] + Qual1 + '\t')
                    outtab.write(ID2 + '\t')
                    outtab.write(Seq2 + '\t')
                    outtab.write(Qual2)
                elif (read == 2):
                    outtab.write(ID1 + '\t')
                    outtab.write(Seq1 + '\t')
                    outtab.write(Qual1 + '\t')
                    outtab.write(ID2 + '\t')
                    outtab.write(tr[0][5:] + Seq2 + '\t')
                    outtab.write(ty[0][5:] + Qual2)
                else:
                    sys.exit("read is not a 1 or 2")

                if final_tags != '':
                    outtab.write('\t' + final_tags)
                outtab.write('\n')
                continue
    else:
        sys.exit("Tab file malformed")

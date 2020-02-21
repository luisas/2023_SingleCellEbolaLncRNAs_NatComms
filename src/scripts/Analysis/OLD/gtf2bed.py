#!/usr/bin/env python
"""GTF2BED"""

import os, sys

i=0
for line in sys.stdin:
  i+=1
  if line.startswith('#'):
    continue

  try:
    contig,pred,f,s,e,dot,strand,score,comment = line.split('\t')
  except:
    sys.stderr.write( "Warning: Wrong line %s: %s\n" % ( i,str(line.split('\t')) ) )
    continue

  if f == "gene":

    #print(len(comment.split('"')))
    g = comment.split('"')[1]
    if(len(comment.split('"'))==11):
      biotype = comment.split('"')[9]
      gene_name = comment.split('"')[5]
    else:
     biotype = comment.split('"')[7]
     gene_name = "NA"
    s,e = int(s),int(e)
    # BED is 0-based, half-open
    s -= 1
    print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % ( contig,str(s),str(e),strand,f,g,gene_name,biotype,pred ))

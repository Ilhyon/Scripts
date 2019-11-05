#!/usr/bin/env python
# -*- coding: utf-8 -*-:

import os
import sys

def do_blast (query_dir, db, outdir):

    filelist = os.listdir(query_dir)
    for queryfile in filelist:
		#~ print str(queryfile)
        out = outdir + queryfile.split('.')[0] + '.blastresult.xml'
        cmd = "blastn -query " + query_dir + queryfile + " -db " + db + " -out "+ out + " -outfmt 5" + " -word_size 4"
        #~ print cmd
        os.system(cmd)

do_blast(sys.argv[1], sys.argv[2], sys.argv[3])

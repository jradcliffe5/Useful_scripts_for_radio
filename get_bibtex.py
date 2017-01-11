#!/usr/bin/python

# Program to fetch bibtex refs from ADS
# On commandline, provide a text file with one bibcode per line

import urllib, string, fileinput, sys

list=[]
outfile=open(sys.argv[1]+'.bibtex', 'w')
for x in fileinput.input(sys.argv[1]):
    y=string.replace(x, "&", "%26")
    print "Fetching "+y
    url='http://cdsads.u-strasbg.fr/cgi-bin/nph-bib_query?bibcode='+y[0:(len(y)-1)]+'&data_type=BIBTEX&db_key=AST%26nocookieset=1'
    ref=urllib.urlopen(url)
    text=ref.readlines()
    # remove CR from list
    text.pop()
    for i in text:
	print i,
	outfile.write(i)
    ref.close()

outfile.close()

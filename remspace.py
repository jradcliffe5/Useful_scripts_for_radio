#!/usr/bin/python

# Program to replace spaces in filenames by "_"
# If no arguments are specified, consider every file
# in current working directory
# List cwd when ready

import sys
import string
import os

if len(sys.argv)>1:
    files=sys.argv[1:]
else:
    files=os.listdir(os.getcwd())

for file in files:
    newname=string.replace(file, " ", "_")
    os.rename(file, newname)

os.system("ls -al")

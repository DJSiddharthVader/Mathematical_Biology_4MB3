#!/bin/bash
if [ $# -eq 0 ]; then
    echo "usage: bash wordcount.sh yourfile.tex"
else
    words="$(detex $1 | wc -w)"
    echo "Document has $words words"
fi

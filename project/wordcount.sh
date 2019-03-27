#!/bin/bash
if [ $# -eq 0 ]; then
    echo "usage: wordcount.sh yourfile.tex"
else
    words="$(detex $1 | wc -l)"
    echo "Document has $words words"
fi

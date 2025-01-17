#!/bin/bash
if [ $# -eq 0 ]; then
    echo "usage: bash wordcount.sh yourfile.tex"
else
    words="$(detex $1  | grep -v "block" | grep -v "blank" | grep -v "draw" | grep -v "node" | grep -v "d[SIRW]" | wc -w)"
    echo "Document has $words words"
fi

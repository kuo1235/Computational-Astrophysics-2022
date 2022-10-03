#!/bin/bash

echo -e "Hello, please type the path: \c" 
read path

i=0
j=0

for something in $path/*; do
    if [ -f "$something" ]
    then
    	echo "$something is a file"
	i=$((i+1))	
    elif [ -d "$something" ]
    then
	echo "$something is a directory"
	j=$((j+1))	
    fi
done

echo "There are $i files and $j dirctories"

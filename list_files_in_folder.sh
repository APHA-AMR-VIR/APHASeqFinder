#!/bin/bash
for file in $(ls reference_key.txt)
do rm $file
done
for line in $(ls references)
do echo $((counter+=1)) $line >> reference_key.txt
done
cat reference_key.txt

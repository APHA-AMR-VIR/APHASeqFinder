#!/bin/bash
number_of_references=$(wc -l < reference_key.txt)
#echo $number_of_references
echo "Pick the reference you want to run from this list:" 
sh list_files_in_folder.sh
read -t 10 reference_pick

if [[ $reference_pick > "$number_of_references" ]]; then echo "Exiting because reference not in list"; exit 1; elif [[ $reference_pick == "1" ]]; then reference="AMRDatabase_20190702_and_EnteroPLasmids_20190514.fna";elif [[ $reference_pick == "2" ]]; then reference="VirulenceFactors_MAbuOun_Aug2015List.fna";elif [[ $reference_pick == "3" ]]; then reference="BacMet_HeavyMetal_resistance_20180803.fas";
 elif [[ -n ${reference_pick//[0-9]/} ]]; then echo "Exiting because input contained non-numerical characters!"; exit 1; elif [[ $reference_pick == "" ]]; then echo "Exiting because database not chosen"; exit 1
fi
echo $reference

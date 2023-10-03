#!/bin/bash

pdb_examples=(
"./input/A1KHF2.pdb"
"./input/pdb1su4.ent"
"./input/pdb1ak2.ent"
"./input/pdb1ake.ent"
"./input/pdb1aky.ent"
"./input/pdb1dvr.ent"
"./input/pdb1e15.ent"
"./input/pdb1kko.ent"
"./input/pdb1n55.ent"
"./input/pdb2ak3.ent"
"./input/pdb2aky.ent"
"./input/pdb2eck.ent"
"./input/pdb3ch0.ent"
"./input/pdb3cwn.ent"
"./input/pdb4ake.ent"
)

idx=1
for entry in "${pdb_examples[@]}";
do
    entry_id=$(basename "$entry" | cut -d'.' --fields 1)
    timestamp=$(date)
    echo "[$timestamp] $idx/${#pdb_examples[@]}:  $entry"

    out_dir="./output/""$entry_id"
    mkdir -p "$out_dir"

    webnma sa -p "$out_dir" "$entry"

    fluc_plot="./output/""$entry_id""/WEBnma_analyses_results/displacement/fluctuations.png"
    mv "$fluc_plot" "./output/""$entry_id""_fluctuations_updated.png"

    ((++idx))
done

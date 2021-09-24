## Helper script for testing webnma3 performance

test_report='./webnma3_performance_report.txt'

######### Testcase 0: 4 PDBs #####
test_dir="../tests/data_profile_alignment/input"
## Run 10 times
count=1
while [[ "$count" -le 10 ]]; do
    echo "=============== Testcase0: comparative analyses: "$count"=============="  >> $test_report
    (time webnma_api ca ${test_dir}/*.pdb -s -p $test_dir) >>$test_report 2>&1
    count=$((count + 1))
done
rm -r ${test_dir}/WEB*


###### Test case 1: 8 pdbs ####
test_dir='./dataset_PDB8'

rm -rf $test_dir
mkdir $test_dir

list='2ak3 4ake 1ak2 1dvr 1ake 2eck 1aky 2aky'
pdbs=""
for pdb in $list
do
    # download pdb and select the chain 'A' / 'B'
    if [ $pdb = "1dvr" ]; then
	webnma_api dl $pdb -c 'B' -p $test_dir   
    else
	webnma_api dl $pdb -c 'A' -p $test_dir
    fi

    # rename the file
    mv ${test_dir}/pdb${pdb}_selected.ent ${test_dir}/${pdb}.ent

    # concat pdb names
    pdbs=${pdbs}' '${test_dir}/${pdb}.ent
done

echo "Data downloading...DONE "

### Run 10 times 
count=1
while [[ "$count" -le 10 ]]; do
    echo "=============== Testcase1: comparative analyses: "$count"=============="  >> $test_report
    (time webnma_api ca ${pdbs} -s -p $test_dir) >>$test_report 2>&1
    count=$((count + 1))
done



# ######## Test case 2: 20 pdbs #####

test_dir="./dataset_PDB20"

### Run 2 times 
count=1
while [[ "$count" -le 2 ]]; do
    echo "=============== Testcase2: comparative analyses: "$count"=============="  >> $test_report
    (time webnma_api ca ${test_dir}/*.pdb -s -p $test_dir) >>$test_report 2>&1
    count=$((count + 1))
done


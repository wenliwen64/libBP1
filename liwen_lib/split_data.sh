#this script is used to split a large directory into small diretories
#!/bin/bash
nfiles=$(ls ./test_100/|wc -l)
echo $nfiles
nparts=$[$nfiles/40]
echo $nparts

II=0
while [ $II -le $nparts ]; do
    mkdir "part_$II"
    find ./test_100/* | head -40 | xargs -i  mv {} "part_$II" 
    let II=$II+1
done

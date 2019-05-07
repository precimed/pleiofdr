#!/bin/bash

#reading trait1
trait1f=./traitFolder/trait1/
trait1n=$(ls $trait1f)

for t1 in $trait1n
do
trait1nb=$(basename "$t1")
IFS='_' read -r -a trait1_array <<< "$trait1nb"
trait1name=${trait1_array[1]}


# reading trait2
trait2f=./traitFolder/trait2/
trait2n=$(ls $trait2f)

# loop starts here
for t2 in $trait2n

do
trait2nb=$(basename "$t2")
IFS='_' read -r -a trait2_array <<< "$trait2nb"
trait2name=${trait2_array[1]}

#Replace the variables
sed -i -e "/traitfolder=/c\traitfolder=./traitFolder"\
    -i -e "/traitfile1=/c\traitfile1="/trait1/"$trait1nb" \
    -i -e "/traitname1=/c\traitname1=$trait1name"\
    -i -e "/traitfiles=/c\traitfiles={\'/trait2/$trait2nb\'}"\
    -i -e "/traitnames=/c\traitnames={\'$trait2name\'}" config.txt

#For cond_FDR
# Output directory
sed -i -e "/outputdir=/c\outputdir="results_"$trait1name"_"$trait2name"_"condfdr/" \
    -i -e "/stattype=/c\stattype=condfdr" \
    -i -e "/fdrthresh=/c\fdrthresh=0.01" config.txt

# Run Matlab now for cond_FDR
matlab -nodesktop -r "run runme.m; exit"

#For conj_FDR
# Output directory
sed -i -e "/outputdir=/c\outputdir="results_"$trait1name"_"$trait2name"_"conjfdr/"\
    -i -e "/stattype=/c\stattype=conjfdr"\
    -i -e "/fdrthresh=/c\fdrthresh=0.05" config.txt

# Run Matlab now for conj_FDR

matlab -nodesktop -r "run runme.m; exit"
done
done

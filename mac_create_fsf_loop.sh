#!/bin/bash

# Subject list for loop
subs=`cat ${1}`          # make sure txt file is full fib file names

for sub in ${subs}
do
echo $sub	
./mac_single_sub_create_fsf.sh ${sub}

echo ‘done for subject: ‘ ${sub}
done


wd="/home/abert/media1/projects/drosophila-evolution/data"

mkdir -p $wd/diyabc_inputs_11
mkdir -p $wd/diyabc_inputs_12

cat $wd/sample_header_AE_admVSnoadm.txt > $wd/diyabc_inputs_11/header.txt
cp $wd/diyabc_inputs_10/POOL_PopData_2.snp $wd/diyabc_inputs_11/

cd $wd/diyabc_executable

./diyabc-RF-linux-v1.1.54 -p $wd/diyabc_inputs_11/ -n "t:100;c:1;s:2011;f:f"
./diyabc-RF-linux-v1.1.54 -p $wd/diyabc_inputs_11/ -R "ALL" -t 100 -r 70000 -g 700 -m

cat $wd/sample_header_AE_whatadm.txt > $wd/diyabc_inputs_12/header.txt
cp $wd/diyabc_inputs_10/POOL_PopData_2.snp $wd/diyabc_inputs_12/

cd $wd/diyabc_executable

./diyabc-RF-linux-v1.1.54 -p $wd/diyabc_inputs_12/ -n "t:100;c:1;s:2011;f:f"
./diyabc-RF-linux-v1.1.54 -p $wd/diyabc_inputs_12/ -R "ALL" -t 100 -r 80000 -g 800 -m

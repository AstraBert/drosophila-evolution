wd="/home/abert/media1/projects/drosophila-evolution/data"

cd $wd/diyabc_executable

./diyabc-RF-linux-v1.1.54 -p $wd/diyabc_noinv_1a/ -n "t:100;c:1;s:2011;f:f"
./diyabc-RF-linux-v1.1.54 -p $wd/diyabc_noinv_1a/ -R "ALL" -t 100 -r 140000 -g 1400 -m

./diyabc-RF-linux-v1.1.54 -p $wd/diyabc_noinv_1b/ -n "t:100;c:1;s:2011;f:f"
./diyabc-RF-linux-v1.1.54 -p $wd/diyabc_noinv_1b/ -R "ALL" -t 100 -r 80000 -g 800 -m

cd /gatk_modified/drosophila-project/

./diyabc-RF-linux-v1.1.54 -p /gatk_modified/drosophila-project/data/pop_data/test_diyabc/ -n "t:50;c:1;s:2002;f:f"
./diyabc-RF-linux-v1.1.54 -p /gatk_modified/drosophila-project/data/pop_data/test_diyabc/ -R "ALL" -m -t 20 -g 100 -r 20000
./diyabc-RF-linux-v1.1.54 -p /gatk_modified/drosophila-project/data/pop_data/test_diyabc/ -d a:pl

for scen in 1 2 3
do
    for param in N1 N2 N3 N4 t34 t23 t12 t32 t24 t21
    do
        ./abcranger-linux-v1.16.69 -t 60 -j 40 --noob 50 --parameter N1 --chosenscen 1
    done
done

./abcranger-linux-v1.16.69 -t 60 -j 40


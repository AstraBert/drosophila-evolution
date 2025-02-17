import random
from math import ceil
import subprocess as sp

f = open("POOL_PopData_1.txt")
lines = f.readlines()
random.shuffle(lines)
if len(lines) < 10000:
	limit = ceil(len(lines)/2)
	first_block = lines[:limit]
	second_block = lines[-limit:]
else:
	first_block = lines[:10000]
	random.shuffle(lines)
	second_block = lines[-10000:]	
print("First block contains:", len(first_block), " SNPs")
print("Second block contains:", len(second_block), " SNPs")
g = open("POOL_PopData_1a.txt","w")
g.writelines(first_block)
g.close()
g = open("POOL_PopData_1b.txt","w")
g.writelines(second_block)
g.close()
sp.run("head -n 2 POOL_PopData.snp > dts1/POOL_PopData_2.snp", shell=True)
sp.run("cat POOL_PopData_1a.txt >> dts1/POOL_PopData_2.snp", shell=True)
sp.run("head -n 2 POOL_PopData.snp > dts2/POOL_PopData_3.snp", shell=True)
sp.run("cat POOL_PopData_1b.txt >> dts2/POOL_PopData_3.snp", shell=True)

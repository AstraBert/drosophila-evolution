# A = Pool5 - EERU_1
# B = Pool19 - ISR
# C = Pool20 - DGN
# O = Pool21 - DrosSim
# X = var

f = open("F4_to_compute.tsv", "w")
f.write("NUM\tDEN\n")

for i in range(1,22):
    if i != 5 and i != 19 and i != 20 and i != 21:
        f.write(f"Pool5,Pool21;Pool20,Pool{i}\tPool5,Pool21;Pool20,Pool19\n")

for i in range(1,22):
    if i != 5 and i != 16 and i != 20 and i != 21:
        f.write(f"Pool5,Pool21;Pool20,Pool{i}\tPool5,Pool21;Pool20,Pool16\n")

for i in range(1,22):
    if i != 5 and i != 17 and i != 20 and i != 21:
        f.write(f"Pool5,Pool21;Pool20,Pool{i}\tPool5,Pool21;Pool20,Pool17\n")

for i in range(1,22):
    if i != 5 and i != 18 and i != 20 and i != 21:
        f.write(f"Pool5,Pool21;Pool20,Pool{i}\tPool5,Pool21;Pool20,Pool18\n")

for i in range(1,22):
    if i != 5 and i != 9 and i != 20 and i != 21:
        f.write(f"Pool5,Pool21;Pool20,Pool{i}\tPool5,Pool21;Pool20,Pool9\n")
f.close()


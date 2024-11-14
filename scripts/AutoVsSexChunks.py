
file_path = "data/freebayes_inputs/all.chunks"
f = open(file_path, "r")
lines = f.readlines()
a = open("data/freebayes_inputs/all.auto.chunks", "w")
x = open("data/freebayes_inputs/all.X.chunks", "w")
region2status =  {"2R": a, "2L": a, "3R": a, "3L": a, "X": x}

for line in lines:
    if line.split(":")[0] in region2status:
        region2status[line.split(":")[0]].write(line)
    else:
        continue 

a.close()
x.close()

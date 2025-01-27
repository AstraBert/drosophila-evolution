input_file = "statobsRF.txt"
output_file = "statobsRF.csv"

f = open(input_file, "r", encoding="utf-8")
lines = f.readlines()
cleaned_lines = [line.strip() for line in lines if line.strip()]
header = ",".join(cleaned_lines[0].split())
values = ",".join(cleaned_lines[1].split())

g = open(output_file, "w", encoding="utf-8")
if header.endswith("\n"):
	g.write(header)
else:
	g.write(header+"\n")
if values.endswith("\n"):
	g.write(values)
else:
	g.write(values+"\n")
g.close()

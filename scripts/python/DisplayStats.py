import pandas as pd
import sys

df = pd.read_csv("statobsRF.csv")
available_stats = list(df.keys())
stats_to_compare = [stat.lower() for stat in available_stats]

print("Hi user! With me, you can display the statitstics in statobsRF.txt :)")
print("1. You can use the name of the statistic to see its value")
print("2. You can use the '/available' command to see the list of available statistics")
print("3. You can just use the '/quit' command to quit the execution of the program")

while True:
	user = input("-->")
	if user.lower() == "/available":
		print(', '.join(available_stats))		
	elif user.lower() == "/quit":
		print("Goodbye! :)")
		sys.exit(0)
	else:
		if user.lower() in stats_to_compare:
			idx = stats_to_compare.index(user.lower())
			val = df[available_stats[idx]].to_list()[0]
			print(f"The value of {available_stats[idx]} is: {val}")
		else:
			print("I cannot recognize your input: please try again!")

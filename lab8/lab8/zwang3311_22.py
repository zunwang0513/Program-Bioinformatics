#!/usr/bin/python
#Exercise 22 Check if braces are paired

#This method count { as 1 and } as -1, if paired the sum should be 0
#If anywhere in the summation process the sum is less than 0, it is not paired.
inputstr = input("Input your braces: ")
lst = [brace for brace in inputstr]
num = 0
for char in lst:
	if char == "{":
		num += 1
	else:
		num -= 1

	if num < 0:
		print("Not Paired")

if num == 0:
	print("Paired")
elif num > 0:
	print("Not Paired")

#!/usr/bin/python
#Exercise 4 2000 to 3000 divisible by 7 not by 5
a = 2000

while a <= 3000:
	if (a % 7 == 0) & (a % 5 != 0):
		print (a)
	a += 1

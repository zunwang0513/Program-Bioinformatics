#!/usr/bin/python
#Exercise 2 Print 10 Fibonacci

count = 0
n1 = 0
n2 = 1
while count < 10:
	n = n1 + n2
	print (n1)
	n1 = n2
	n2 = n
	count += 1

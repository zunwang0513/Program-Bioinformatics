#!/usr/bin/python
#Exercise 14 Print List of Fibonacci

count = int(input("Input the number of Finonacci number you want: "))
lst = []
n1 = 0
n2 = 1
while count > 0:
	lst.append(n1)
	nth = n1 + n2
	n1 = n2
	n2 = nth
	count -= 1
print(lst)

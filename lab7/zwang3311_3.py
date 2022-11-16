#!/usr/bin/python
#Exercise 3 factorial

num = int(input("Enter your value:"))
product = 1
while (num > 0):
	product *= num
	num -= 1
print (product)

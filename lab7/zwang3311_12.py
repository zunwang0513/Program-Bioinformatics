#!/usr/bin/python
#Exercise 12 Cut string in Half

string = input("Input your string: ")

if len(string) % 2:
	halfidx = (len(string) + 1) / 2
else:
	halfidx = len(string) / 2

print (string[0:int(halfidx)])

#!/usr/bin/python
#Exercise 23 Detect palindrom

word = input("Please input a word: ")

if word == word[::-1]:
	print("It is a palindrom")
else:
	print("It is not a palindrom")

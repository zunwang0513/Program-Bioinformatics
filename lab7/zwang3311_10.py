#!/usr/bin/python
#Exercise 10 Ticketing when birthday

print ("Speed limit is 50 mph")
speedlimit = 50
speed = int(input("Current speed is: "))
birthday = input("Is this birthday? y/n: ")
fine = "none"

if birthday == "y":
	speed -= 5

if speed - speedlimit > 15:
	fine = "large"
elif (speed - speedlimit < 15) & (speed - speedlimit > 0):
	fine = "small"

print (fine) 


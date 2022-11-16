#!/usr/bin/python
#Excercise 8 Sum average max and min of all numbers input

num = int(input("input your number: "))
sumofnum = num
count = 1
maxofnum = num
minofnum = num
while num != 0:
	num = int(input("input your number: "))
	sumofnum += num
	average = sumofnum / count
	count += 1
	if maxofnum < num:
		maxofnum = num
	if (minofnum > num) & (num != 0):
		minofnum = num

print ("sum is " + str(sumofnum))
print ("average is " + str(average))
print ("maximum is " + str(maxofnum))
print ("minimum is " + str(minofnum))

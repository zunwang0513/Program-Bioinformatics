#!/usr/bin/python
#Exercise 23 Detect DNA palindromes

sequencestr = input("Please input a DNA sequence: ")

#Turn it into a list
sequence = [char for char in sequencestr]

#temp make sure this sequence only contains ATCG
temp = 0
for base in sequence:
    if (base != "A") and (base != "G") and (base != "C") and (base != "T"):
        temp -= 1
if temp < 0:
    print("Not DNA sequence")
else:
    for baseidx in range(len(sequence) // 2):
        if sequence[baseidx] == "A":
            sequence[baseidx] = "T"
        elif sequence[baseidx] == "T":
            sequence[baseidx] = "A"
        elif sequence[baseidx] == "C":
            sequence[baseidx] = "G"
        elif sequence[baseidx] == "G":
            sequence[baseidx] = "C"
    if sequence == sequence [::-1]:
        print("This is a biological palindrom")
    else:
        print("This is not a biological palindrom")

#Change the first half of the list to its pairing DNA and see if the string turn out to be a palindrome

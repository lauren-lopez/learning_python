#!/usr/bin/env python3

# You are probably well aware of the 'birthday paradox'
# https://en.wikipedia.org/wiki/Birthday_problem
# Let's try simulating it
# We will have a variable number of bins (can be months or days)
# And some number of trials for the simulation
# And some number of people whose have random birthdays
# Use assert() to check parameters
# On the command line:
#	python3 birthday.py <bins> <trials> <people>

import sys
import random

assert(len(sys.argv) == 4)

bins = int(sys.argv[1])
trials = int(sys.argv[2])
people = int(sys.argv[3])

assert(bins > 0)
assert(trials > 0)
assert(people > 1)

collisions = 0                        #when two or more birthdays occur on the same day

# create an empty calendar
for t in range(trials):               # outermost loop needs to be number of trials
    #calendar = []                    # this is one way of setting the calendar. Simpelr alternative is below
    #for i in range(bins): calendar.append(0)
    calendar = [0] * bins             # this sets a 0 value for the number of bins (e.g. number of days, months etc)
    
    # insert people into calendar
    for p in range(people):           #populate the calendar with the defined number of people
        r = random.randint(0, bins-1) #defines a random position in the calendar to assign a person to
        calendar[r] += 1              #adds person to that date in the calendar (currently 0)
    
    # find shared birthdays
    same_day = False                  # another way to do this would be to set = 0, and add one each time
    for day in calendar:
        if day > 1: 
            same_day = True           #if more than one birthday on the same day, condition becomes true
            break                     # breaks the loop after True is met - may save time
    if same_day: collisions += 1      #adds one to collision count for each shared birthday
    
print(collisions/trials)              #to increase the precision of this number, increase the number of trials in the command line

"""
python3 birthday.py 365 1000 23
0.520
"""


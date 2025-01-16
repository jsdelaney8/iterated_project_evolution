### Code for paper - "A Simplified Model of Iterative Compound Optimization"
### Date : 14 October 2024
### Author : John S Delaney
### ORCHID iD: https://orcid.org/0000-0003-2218-1167
###
### Requires Python3.6 or higher (uses F-strings and the 'choices' method from module random)

import random as r
import numpy as np

# Number of set bits required in the fingerprint, number of activebits at the start of the fingerprint and
# the number of "kill" bits at the end
censor = 0.1 #0.05
nactivebits = 15
num_set = 40 # 56 roughly the mean number set in EPS H projects
nkillbits = 5

bitvectorlength = 256
killcutoff = bitvectorlength - nkillbits + 1

# Set the number of projects to generate and the size of each project
nprojects = 100 #100 usual
projectsize = 1000 #1000 usual

def main():
    # Seed the random number generator
    r.seed()

    # Make a header line for the csv file
    header = "prj,activity,activity_no_kill,seq," + ",".join([f'f{i}' for i in range(1, bitvectorlength+1)]) + ",null"
    print(header)

    # Loop over projects (nprojects projects)
    for prj in range(1, nprojects+1):

        # Initialise scores
        scores = [0] * (bitvectorlength+1)
        act_list = []
        seq = 0

        # Setup an weakly active compound as the lead (was always the same, now randomised to an extent)
        fp = initial_weak_active_setup()

        # Loop through a project, project size set to projectsize compounds
        while seq < projectsize:
            # Measure the activity of the compound - with and without kill bits
            active, active_no_kill = measure_activity(fp)

            # Output the keys as set bits and the whole string as textual fingerprint
            # Random censor - increases the spacing between compounds
            if seq == 0:
                seq += 1
                act_list.append(active)
                output_fingerprint(fp,prj,seq,active,active_no_kill)
            else:
                if r.random() < censor:
                    seq += 1
                    act_list.append(active)
                    output_fingerprint(fp,prj,seq,active,active_no_kill)

            # Increment score array
            for i in fp:
                scores[i] += active

            # Pick a bit position change candidate, ensuring it is not already set
            selected = bit_position_change_candidate(fp)

            # Find max score and freeze that bit position (protect variable)
            protect,max_score = freeze_top_score_position(scores)

            # Choose one of the ten set to change, ensuring it is not protected
            random_index = choose_bit_to_change(fp,protect)

            fp[random_index] = selected
        

def initial_weak_active_setup() :
    fp=[]
    fp.append(r.randint(1,nactivebits))
    choice = r.choices(list(range(nactivebits+1,bitvectorlength-nkillbits)) , k=num_set-1)
    for i in choice:
        fp.append(i)
    return fp

def choose_bit_to_change(fp,protect):
    random_index = None
    while random_index is None or fp[random_index] == protect:
        random_index = r.randint(0, num_set - 1)
    return random_index

def freeze_top_score_position(scores):
    max_score = -1
    for i in range(1, bitvectorlength+1):
        if scores[i] > max_score:
            protect = i
            max_score = scores[protect]
    return protect,max_score

def bit_position_change_candidate(fp):
    selected = None
    while selected is None or selected in fp:
        selected = r.randint(1, bitvectorlength)
    return selected

def measure_activity(fp):
    # Assess activity
    active = 0
    for i in fp:
        if i <= nactivebits:
           active = active+1
    active_no_kill = active
    # Apply kill bit if present
    for i in fp:
        if i >= killcutoff:
           active = 0
    return active,active_no_kill

def output_fingerprint(fp,prj,seq,active,active_no_kill):
    # Initialise text fingerprint
    blank = "0" * bitvectorlength
    # Output the keys as set bits and the whole string as textual fingerprint
    txtfp = [char for char in blank]
    for bit in fp:
        txtfp[bit-1] = "1"
    print(f"p{prj},{active},{active_no_kill},{seq},{','.join(txtfp)}")

if __name__ == "__main__":
    main()
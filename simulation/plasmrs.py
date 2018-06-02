#!/usr/bin/env python

# Plasmid Mutation Simulator

import random
import time
import argparse
import os.path

# Cell population stored as dictionary of lists:
#   Keys are tuples as number of (FP, MP) in that category of cell
#   Values are lists with (Fitness, Selection Probability, Cell Number)

parser = argparse.ArgumentParser(description='Simulate plasmid mutations')

parser.add_argument('--initial-cells', '-i', dest='initialPopulationSize', type=int,
                   default=1, help='initial number of cells')

parser.add_argument('--final-cells', '-f', dest='finalPopulationSize', type=int,
                   default=10000, help='final number of cells')

parser.add_argument('--plasmid-copy-number', '-p', dest='plasmidsPerCell', type=float,
                   default=100, help='mean number of plasmids per cell')

parser.add_argument('--plasmid-mu', '-u', dest='plasmidMutationRate', type=float,
                   default=1E-6, help='mutation rate per plasmid replication')

parser.add_argument('--plasmid-segregation-model', '-m', dest='plasmidReplicationSegregationModel', type=float,
                   default=1, help='segregation model for plasmids (1=replicate plasmids and eqully divide into daughter cells)')

parser.add_argument('--replicates', '-r', dest='replicates', type=int,
                   default=100, help='replicates')

parser.add_argument('--random-seed', '-s', dest='seednum', type=int,
                   default=-1, help='random number seed (Default = set based on time)')

parser.add_argument('--mutant-count-output', '-mo', dest='mutant_count_file_name', type=str,
                   default="mutant_counts.csv", help='Output file')
                   
parser.add_argument('--population-output', '-po', dest='population_file_name', type=str,
                   default="population.csv", help='Output file')

args = parser.parse_args()


ancestralPlasmidWithinCellFitness = 1
mutantPlasmidWithinCellFitness = 1

# Future extension to vary fitness depending 
# on how many plasmids are present
def computeCellFitness(cellType):
  return 1

def assignCellSelectionProb(states):
  popFit = 0
  #Add up all fitness values
  for key in states:
    states[key][1] = states[key][0] * states[key][2]
    popFit += states[key][1]
  #Assign probabilities
  for key in states:
    states[key][1] = float(states[key][1]) /  popFit
  return states

def pickCellToDivide(states):
  randNum = random.uniform(0, 1)
  s = 0
  for key in states:
    s += float(states[key][1])
    if s >= randNum:
      return key
  print(randNum)
  
# cell is a list with two values
def replicatePlasmidsInCell(cell,numNewPlasmids, plasmidMutationRate):

  pFitTotal = 0 # Sum total fitness of plasmids in cells
  pFitTotal += cell[0] * ancestralPlasmidWithinCellFitness 
  pFitTotal += cell[1] * mutantPlasmidWithinCellFitness

  if (pFitTotal == 0):
    return

  #print str(cell) + " " + str(pFitTotal) + "\n"

  APfit = (float(ancestralPlasmidWithinCellFitness) * cell[0]) / pFitTotal # Determine relative fitnesses of plasmids in cell
  MPfit = (float(mutantPlasmidWithinCellFitness) * cell[1]) / pFitTotal
  
  #print str(FPfit) + " " + str(MPfit) + " " + str(DPfit) + "\n"

  pTot = cell[0] + cell[1]
  for plasmid in range(0,numNewPlasmids):
    randNum = random.uniform(0, 1)
    if MPfit > randNum:
      cell[1] += 1
    else:
      randNum = random.uniform(0, 1)
      if plasmidMutationRate > randNum:
        cell[1] += 1
      else:
        cell[0] += 1  

def divide(cell, states):

  states[cell][2] -= 1 # Subtract the cell about to reproduce from the population
  if (states[cell][2] == 0):
    del states[cell]

  cell = list(cell) # Convert the tuple to a list

  #Replicate the plasmids
  
  ##For methods that replicate before segregation
  if (args.plasmidReplicationSegregationModel == 1) or (args.plasmidReplicationSegregationModel == 2):
    replicatePlasmidsInCell(cell,args.plasmidsPerCell, args.plasmidMutationRate)

  daughter1 = [0, 0] # Make two empty daughter cells
  daughter2 = [0, 0]

### Divide plasmids
  plasmidList = []
  plasmidList += cell[0] * [0]
  plasmidList += cell[1] * [1]
    
  # Equal
  if (args.plasmidReplicationSegregationModel == 1):
    random.shuffle(plasmidList)
    for p in range(0,args.plasmidsPerCell):
      daughter1[plasmidList[p]] += 1
    for p in range(args.plasmidsPerCell,2*args.plasmidsPerCell):
      daughter2[plasmidList[p]] += 1  
  # Coin flip
  elif (args.plasmidReplicationSegregationModel == 2) or (args.plasmidReplicationSegregationModel == 4):
    for p in plasmidList:
      randNum = random.uniform(0, 1)
      if randNum < 0.5:
        daughter1[p] += 1
      else:
        daughter2[p] += 1
  else:
    exit()
  
    #For methods that replicate AFTER segregate
  if (args.plasmidReplicationSegregationModel == 3) or (args.plasmidReplicationSegregationModel == 4):
    replicatePlasmidsInCell(daughter1,args.plasmidsPerCell-daughter1[0]-daughter1[1], args.plasmidMutationRate)
    replicatePlasmidsInCell(daughter2,args.plasmidsPerCell-daughter2[0]-daughter2[1], args.plasmidMutationRate)
  
  daughter1 = tuple(daughter1)
  daughter2 = tuple(daughter2)

  if daughter1 in states:
    states[daughter1][2] += 1 # If that cell type is already in states, add one to the count
  else:
    states[daughter1] = [computeCellFitness(daughter1), None, 1] # Otherwise, new key and compute fitness

  if daughter2 in states:
    states[daughter2][2] += 1 # If that cell type is already in states, add one to the count
  else:
    states[daughter2] = [computeCellFitness(daughter1), None, 1] # Otherwise, new key and compute fitness


def main():
  
  if (args.seednum == -1):
    t = int( time.time() * 1000.0 )
    seed = ( ((t & 0xff000000) >> 24) +
     ((t & 0x00ff0000) >>  8) +
     ((t & 0x0000ff00) <<  8) +
     ((t & 0x000000ff) << 24)   )
    random.seed(seed)
    print("Random seed based on time: " + int(seed))
  else:
    random.seed(args.seednum)
    print("Used defined random seed: " + int(args.seednum))

  
  #Open in append mode
  mutant_count_file_exists = os.path.isfile(args.mutant_count_file_name) 
  mutant_count_file = open(args.mutant_count_file_name, "w+")
  if not mutant_count_file_exists:
    mutant_count_file.write("N_0,N,N_p,u_p,repl,m\n")
  
  population_file_exists = os.path.isfile(args.population_file_name) 
  population_file = open(args.population_file_name, "w+")
  if  not population_file_exists:
    population_file.write("N_0,N,N_p,u_p,repl,p_wt,p_mut,cells\n")

  
  for i in range(1, args.replicates+1):
    print("replicate: " + str(i))
    states = {} # Create a dictionary of "cell states" that keeps track of cells in population


    
    # In state dictionary:
    # Each key will be plasmid count: (ancestral plasmids, mutated plasmids)
    # Each value will be a list: [the computed fitness of cells with that plasmid count, fitness relative to population, and the number of cells of that type].
    populationSize = args.initialPopulationSize
    states[(args.plasmidsPerCell, 0)] = [computeCellFitness((args.plasmidsPerCell, 0,)), None, int(args.initialPopulationSize)]

    while(populationSize < args.finalPopulationSize):
      populationSize = populationSize + 1
      states = assignCellSelectionProb(states)
      cellToDiv = pickCellToDivide(states)
      divide(cellToDiv, states)

    m = 0
    for key in states:
      print(str(key) + "\n" + str(states[key]) + "\n")
      
      population_file.write(str(args.initialPopulationSize) + "," 
        + str(args.finalPopulationSize) + "," + str(args.plasmidsPerCell) + "," 
        + str(args.plasmidMutationRate) + "," + str(i) + "," 
        + str(key[0]) + "," + str(key[1]) + "," + str(states[key][2]) + "\n")

      if (key[1] > 0):
        m = m + states[key][2]
    mutant_count_file.write(str(args.initialPopulationSize) + "," 
      + str(args.finalPopulationSize) + "," + str(args.plasmidsPerCell) + "," 
      + str(args.plasmidMutationRate) + "," + str(i) + "," + str(m) + "\n")
      
    population_file.flush()
    mutant_count_file.flush()
    

main()

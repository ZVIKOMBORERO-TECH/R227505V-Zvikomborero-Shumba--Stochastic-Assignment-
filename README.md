# HASTS416/RM Markov Chain Questions

## Overview
This repository contains solutions for the Markov Chain questions with three parts:
- Q1/A1: 5-State Markov Chain Analysis
- Q2/A2: 7-State Markov Chain Analysis  
- Q3/A3: Non-Homogeneous Traffic Model

## Files
- `assignment_markov_analysis.R` - Main analysis code with complete implementations
- `run_assignment.R` - Script to execute the analysis
- `A1_*.png` - Generated visualizations for Question 1 (Q1/A1)
- `A2_*.png` - Generated visualizations for Question 2 (Q2/A2)

## Usage
To run the complete analysis:

```r
# Set working directory to the question folder
setwd("path/to/R227505V-Zvikomborero-Shumba--Stochastic-Assignment-")

# Execute the analysis
source("run_assignment.R")
```

**Note**: Make sure to replace the path with the actual location where you downloaded the repository.

## Requirements
- R (version 4.0+)
- Required packages: base R packages only

## Question Structure
The analysis includes:
1. **Shared utility functions** for matrix operations and Markov chain simulation
2. **Question 1 (Q1/A1)**: 5-state Markov chain with convergence analysis and trajectory visualization
3. **Question 2 (Q2/A2)**: 7-state Markov chain with similar analysis
4. **Question 3 (Q3/A3)**: Non-homogeneous traffic flow model

## Output Files
- Convergence plots showing marginal probabilities over time
- State transition diagrams
- Multiple trajectory simulations

## Author
Student Name: Zvikomborero Shumba
Student ID: R227505V
Course: HASTS416/RM

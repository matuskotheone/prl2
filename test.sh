#!/bin/bash

# Check if the number of arguments is correct
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <argument 1> <argument 2>"
    exit 1
fi

# Assign command-line arguments to variables
numbers=4
arg1=$1
arg2=$2

# Compile MPI program
mpic++ --prefix /usr/local/share/OpenMPI -o life life.cpp

# Run MPI program with specified number of processors and arguments
mpirun --prefix /usr/local/share/OpenMPI -np "$numbers" life "$arg1" "$arg2"

# Clean up
rm -f life

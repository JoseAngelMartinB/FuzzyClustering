# FuzzyClustering

## Description
Development of a Java program to calculate the **fuzzy communities** of each vertex of a given graph, using an approach that allows each vertex of the graph to belong to multiple communities at the same time.

## Author
Developed by **José Ángel Martín Baos** for **SciCom (Scientific Computing Group)** of the University of Castilla-La Mancha in July 2016.

## Input data
The input data should be the adjacency matrix of the desired graph, and should be stored in a separated file.

## Some extra information
* The file `FuzzyClustering_2.java` is a improved version of the algoritm that allows parrallel computing.
* The method `findCommunities_GC` implements the same functionality of as the method `findCommunities` but using CSR (Compressed Sparse Row) format to represent the adjacency matrix, with the consecuent improvement in the memory, while executing the program.

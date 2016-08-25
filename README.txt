*********************************************************************************************

  FlexiCoClustering - A flexible model-based clustering of mixed binary and continuous data
               Fatin Zainul Abidin & David Westhead: www.bioinformatics...
               
*********************************************************************************************
|                                                                                           |
|     FlexiCoClustering is made available under a GPL license, for full details of the      |
|     license see the included LICENSE text file.                                           |
|                                                                                           |
|-------------------------------------------------------------------------------------------|

Contents of release:

    FlexiCoClustering 
    | FlexiCoClustering.jar - the main executable file. See usage below.
    | src - directory containing all the java source.
    | javadoc - directory containing the javadoc for this release.
    | README - this file.
    | LICENSE - details of the GPL license this software is covered by.
    | example - directory containing example file required to run this program.
    | ext - directory containing all external libraries reguired for this program.
    | Documentation - the FlexiCoClustering manual file

Usage:

    To run the demo (using JDK 1.8.0), use the following command:

    On Windows command promt
    java -jar <path to the FlexiCoClustering.jar>\FlexiCoClustering.jar

    On any Linux/Unix distributions
    java -jar <path to the FlexiCoClustering.jar>/FlexiCoClustering.jar

    On the GUI options:

    1.  "Name a clustfile" : Name the modfile as i.e. Modfile.txt
        "Name an outfile"  : Name an outfile as i.e. Output.txt
        "Select inputFile" : Select <path to the FlexiCoClustering.jar>\exampla\input.txt
    2.  Press the "Execute" button. On default, this program will run for 100 iterations
        (MaxTemps) and produces and update heatmap image files (.png) in every 10 iterations 
        interval for binary and continuous variables. The Output.txt and Modfile.txt 
        will be updated as the program progress and once the program terminated respectively.
    3.  If more then the initially specified number of iterations is required after step 2 had
        finished, re run the program by first inserting current run step number in the NRun
        and then press the "Execute" button again.
        
Input file (i.e. Input.txt):

    A space " " separated file.

    +--------------------------------------------------------------------------+
    |BinaryInputs: FRYL OR13H1 PRPF8 CACNA2D3 SEMA4A TUBA3C PPP1R3A RYR3 FOXP1 |
    |ContinuousInputs: ELANE ZFY JPH1 LIFR ASS1 DUSP27 LHX6 NAPSB SLC28A3      |
    |TCGA1 0 0 0 1 0 0 0 0 1 3.41 2.81 8.91 9.91 13.44 2.33 9.09 9.82 0.01     |
    |TCGA2 1 0 0 0 0 0 0 0 1 3.41 2.81 8.91 9.91 13.44 2.33 9.09 9.82 0.01     |
    |TCGA3 0 0 0 0 0 1 0 0 1 3.41 2.81 8.91 9.91 13.44 2.33 9.09 9.82 0.01     |
    |TCGA4 0 0 0 1 0 0 0 0 0 3.41 2.81 8.91 9.91 13.44 2.33 9.09 9.82 0.01     |
    |TCGA5 0 0 0 0 1 0 0 0 1 3.41 2.81 8.91 9.91 13.44 2.33 9.09 9.82 0.01     |
    +--------------------------------------------------------------------------+
    Third row onwards: 
    1. First column: The data points (i.e. gene names, sample names) 
    2. Second column onwards are the binary inputs ("1" and "0") followed by continulous numbers

*********************************************************************************************

 FlexiCoClustering v.0.1        Fatin Zainul Abidin & David Westhead: www.bioinformatics...
                                                                           
*********************************************************************************************
|                                                                                           |
|     FlexiCoClustering is made available under a GPL license, for full details of the      |
|     license see the included LICENSE text file.                                           |
|                                                                                           |
|-------------------------------------------------------------------------------------------|


For the latest version of FlexiCoClustering and any tutorials and examples that exist, see 
the FlexiCoClustering project page at http://www.bioinformatics.... 


Contents of release:

    FlexiCoClustering v.0.1
    | FlexiCoClustering-v.0.1.jar - the library. Put this anywhere on your class path.
    | src - directory containing all the java source.
    | javadoc - directory containing the javadoc for this release.
    | README - this file.
    | LICENSE - details of the GPL license this software is covered by.
    | example - directory containing example files required to run this program.
    | ext - directory containing all external libraries reguired for this program.
    | Documentation - the FlexiCoClustering v.0.1 manual file


Usage:

    To run the demo (using JDK 1.6.0 or later), use the following command:

    java -jar <path to the FlexiCoClustering-v.0.1.jar>\FlexiCoClustering-v.0.1.jar
    
    On the GUI options:

    1.  "Name a clustfile" : Name the modfile as i.e. Modfile.txt
        "Name an outfile"  : Name an outfile as i.e. Output.txt
        "Select inputFile" : Select <path to the FlexiCoClustering-v.0.1.jar>\input.txt
    2.  Press the "Execute button". This program will run for 100 iterations and produces heatmap image file in every 10 iterations interval for binary inputs and continuous inputs. The Output.txt and Modfile.txt will be updated as the program progress and once the program terminated respectively.

   
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

   Third line onwards:
       1.  First column: The data points names (i.e. gene names, sample names)
       2.  Second column onwards are the binary inputs ("1" and "0", continulous numbers)



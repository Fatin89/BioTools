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

    To run the demo , use the following command:

    java  <path to the FlexiCoClustering>.ClustProg <Input.txt> <Output.txt>
    
    Press enter
    
    To re-submit the code after program terminated:
    1.  Change Nrun: '0' to Nrun: '1'
    2.  Increase the MaxTemps to a higher value
    
   
Input file (i.e. Input.txt):

    A space " " separated file.

    +--------------------------------------------------------------------------+
    |     |
    |TCGA1 0 0 0 1 0 0 0 0 1 3.41 2.81 8.91 9.91 13.44 2.33 9.09 9.82 0.01     |
    |TCGA2 1 0 0 0 0 0 0 0 1 3.41 2.81 8.91 9.91 13.44 2.33 9.09 9.82 0.01     |
    |TCGA3 0 0 0 0 0 1 0 0 1 3.41 2.81 8.91 9.91 13.44 2.33 9.09 9.82 0.01     |
    |TCGA4 0 0 0 1 0 0 0 0 0 3.41 2.81 8.91 9.91 13.44 2.33 9.09 9.82 0.01     |
    |TCGA5 0 0 0 0 1 0 0 0 1 3.41 2.81 8.91 9.91 13.44 2.33 9.09 9.82 0.01     |
    +--------------------------------------------------------------------------+

   Third line onwards:
       1.  First column: The data points names (i.e. gene names, sample names)
       2.  Second column onwards are the binary inputs ("1" and "0", continulous numbers)


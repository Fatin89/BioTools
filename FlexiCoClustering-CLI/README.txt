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

    FlexiCoClustering-CLI
    | FlexiCoClustering- the package folder
    | README - this file.
    | LICENSE - details of the GPL license this software is covered by.
    | Documentation - the FlexiCoClustering manual file


Usage:

    To run the demo , use the following command:

    java  <path to the FlexiCoClustering.jar>/FlexiCoClustering.jar <input.txt> <Output.txt>
    
    Press enter
    
    To re-submit the code after program terminated:
    1.  Change Nrun: '0' to Nrun: '1'
    2.  Increase the MaxTemps to a higher value
    
   
Input file (i.e. Input.txt):

    A space " " separated file.

    +--------------------------------------------------------------------------+
    | NItems: 170                                                              |
    | NBinary: 9                                                               |
    | NContinuous: 9                                                           |
    | Agglomerative                                                            |
    | IC: AIC                                                                  |
    | MergeSplitProbability: 0.25                                              |
    | MaximumIterations: 100                                                   |
    | StartTemp: 500                                                           |
    | TempFactor: 0.999                                                        |
    | MaxTemps: 1000                                                           |
    | MaxRepIters: 2000                                                        |
    | Seed: 1                                                                  |
    | EMIterations: 100                                                        |
    | OutInterval: 1                                                           |
    | ClustFile: clustfile.txt                                                 |
    | SinglePoint: 0                                                           |
    | NClusters: 0                                                             |
    | NormExp: 1                                                               |
    | MixCoefficient: 1                                                        |
    | Nrun: 0                                                                  |
    |TCGA1 0 0 0 1 0 0 0 0 1 3.41 2.81 8.91 9.91 13.44 2.33 9.09 9.82 0.01     |
    |TCGA2 1 0 0 0 0 0 0 0 1 3.41 2.81 8.91 9.91 13.44 2.33 9.09 9.82 0.01     |
    |TCGA3 0 0 0 0 0 1 0 0 1 3.41 2.81 8.91 9.91 13.44 2.33 9.09 9.82 0.01     |
    |TCGA4 0 0 0 1 0 0 0 0 0 3.41 2.81 8.91 9.91 13.44 2.33 9.09 9.82 0.01     |
    |TCGA5 0 0 0 0 1 0 0 0 1 3.41 2.81 8.91 9.91 13.44 2.33 9.09 9.82 0.01     |
    +--------------------------------------------------------------------------+

   Third line onwards:
       1.  First column: The data points names (i.e. gene names, sample names)
       2.  Second column onwards are the binary inputs ("1" and "0", continulous numbers)


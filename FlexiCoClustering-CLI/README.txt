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
    | input.txt - a demo input text file


Usage:

    To run the demo , use the following command:

    java  <path to the FlexiCoClustering.jar>/FlexiCoClustering.jar <input.txt> <Output.txt>
    
    Press enter
    
    To re-submit the code after program terminated:
    1.  Change Nrun: '0' to Nrun: '1'
    2.  Increase the MaxTemps to a higher value than previous run (if the MaxTemps iteration was completed) or else, just use default MaxTemps.
    3.  The program can be terminated at any time by pressing Ctrl+C
    
   
Input file (i.e. input.txt):

A space " " separated file.
-----------------------------------------------------------------------
NItems: 200
NBinary: 20
NContinuous: 20
Agglomerative
IC: AIC
MergeSplitProbability: 0.25
MaximumIterations: 100
StartTemp: 500
TempFactor: 0.99
MaxTemps: 100000
MaxRepIters: 100
Seed: 1
EMIterations: 100
OutInterval: 1
ClustFile: clustfile.txt
NormExp: 0
Nrun: 0  
GENE0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 1 1 0 1 1 0 -0.4765386533 -1.20066852 -0.422161301 -0.4769262762 -0.3915312885 1.6614625942 -0.5264303615 0.5148318734 -0.7380265044 0.7983998534 1.6788213258 0.3525189371 -0.5343971289 0.8074409234 -1.1545894616 1.1516130154 -1.313722178 1.6782261766 -1.0162295721 -0.3512150164
GENE1 0 0 0 1 1 1 1 0 0 0 0 0 0 0 1 1 0 1 1 0 -0.4808406431 -1.1876763307 -0.4423410512 -0.4751493839 -0.396965996 1.6709779194 -0.526275095 0.5121057944 -0.7399474163 0.8027740397 1.6742977045 0.3562186725 -0.5220498526 0.7931213765 -1.1539646663 1.1424829456 -1.3068357316 1.6858205787 -1.0318587232 -0.3533332355
GENE2 0 0 0 1 1 1 1 0 0 0 0 0 0 0 1 1 0 1 1 0 -0.4877185071 -1.1811041112 -0.4370126771 -0.4743384742 -0.4249390029 1.6630460356 -0.5143488862 0.5137408259 -0.7354947281 0.8006116622 1.6781314018 0.3520833077 -0.5260607663 0.7879168079 -1.154671976 1.1459990165 -1.3188525326 1.6714563372 -1.044445159 -0.3566300553
GENE3 0 0 0 1 1 1 1 0 0 0 0 0 0 0 1 1 0 1 1 0 -0.4805147615 -1.1938070599 -0.4334079898 -0.4817564761 -0.3792389503 1.6778144603 -0.5152447322 0.5473686209 -0.7566644854 0.8010619739 1.6901398941 0.3370148097 -0.5116772569 0.8111656649 -1.1637947577 1.131828782 -1.3101856866 1.6750770174 -1.0320545106 -0.3419100239
GENE4 0 0 0 1 1 1 1 0 0 0 0 0 0 0 1 1 0 1 1 0 -0.4954518968 -1.2091381875 -0.4345381315 -0.4727608557 -0.3922457347 1.6801477645 -0.5250288841 0.4939493614 -0.7337396017 0.8093641675 1.6790758671 0.3440635754 -0.5116291603 0.8148480595 -1.14248843 1.1433251814 -1.3061449941 1.663806288 -1.0192842641 -0.3672015256
GENE5 0 0 0 1 1 1 1 0 0 0 0 0 0 0 1 1 0 1 1 0 -0.4911375599 -1.1879389496 -0.4386060351 -0.4864677453 -0.4001632414 1.6651027605 -0.5332935566 0.5170789559 -0.7514515822 0.8180966415 1.6853029648 0.3478484738 -0.5372745165 0.7983662732 -1.1589839099 1.1460077692 -1.3206090553 1.6758920674 -1.0344307014 -0.3657283278
GENE6 0 0 0 1 1 1 1 0 0 0 0 0 0 0 1 1 0 1 1 0 -0.4838485309 -1.2002114897 -0.4423914419 -0.4924432623 -0.3876738113 1.6971197285 -0.5438864896 0.5147098215 -0.7462507198 0.8138986462 1.676762964 0.3441529717 -0.52999529 0.7997874468 -1.1431666742 1.1527270638 -1.331168246 1.6877146448 -1.030679756 -0.3600469108
GENE7 0 0 0 1 1 1 1 0 0 0 0 0 0 0 1 1 0 1 1 0 -0.4897891989 -1.1889511486 -0.4456117067 -0.4717567455 -0.3914853631 1.6945979498 -0.5326287141 0.4935343707 -0.7427462324 0.8139988046 1.6744598616 0.3427291599 -0.5209096308 0.8071768653 -1.1632877862 1.140879741 -1.3144533559 1.6775802363 -1.0152447552 -0.3544501325
                                        .
                                        .
                                        .
-----------------------------------------------------------------------------
      18th line onwards
       1.  First column: The data points names (i.e. gene names, sample names)
       2.  Second column onwards are the binary inputs ("1" and "0") followed by continuous values    
     


/*
 *  Copyright 2016 Fatin Zainul Abidin & David Westhead (www.bioinformatics.leeds.ac.uk)
 *  Licensed under GNU Lesser General Public License
 * 
 *  This file is part of FlexiCoClustering - The model-based flexible co-clustering api for Java.
 *
 *  FlexiCoClustering is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published 
 *  by the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  FlexiCoClustering is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with FlexiCoClustering.  If not, see <http://www.gnu.org/licenses/>.
 */
 


package FlexiCoClustering;



import java.util.*;

import java.io.*;

import java.util.logging.Level;

import java.util.logging.Logger;



public class EMRefiner {

    private int [][] reg_inputs; /* array holding regulatory input information, dimension is number of genes by number of inputs */

    private double [][] exp_patterns; /* array holding expression patterns dimension number of genes by number of 'time' points */

    private int ngen,ninps,npnts; /* number of genes, number of reg inputs, number of expression 'time' points */

    private int nmods; /* number of modules */



/* instance methods and variables for refinement using the expectation maximization algorithm (EM) */

    private double [][] margin_density;

    private double [][][] normpars,old_normpars,new_normpars;

    private double [][] bern_ps,old_bern_ps,new_bern_ps;

    private double [] mix_pars,old_mix_pars,new_mix_pars;

    private double [][] prob_density;

    private int [] mpm_list;

    private int em_iters;

    private static final double SMALL_VAR=0.00001;

    private static final double SMALL=0.00001;





    

    

    

    public EMRefiner(double [][] exp_patrns, int [][] reginps, int ngenes, int ninputs, int nexp, int nmodules, double [] mix,

                                double [][] bernp, double [][][] npars, int iterations) {

        ngen = ngenes;

        ninps = ninputs;

        npnts = nexp;

        /* no need to copy these arrays - we will not change them */

        exp_patterns = exp_patrns;

        reg_inputs = reginps;

        nmods = nmodules;

        em_iters = iterations;

        /* set up arrays holding parameters and marginal density */

        bern_ps = new double [nmods][ninps];

        normpars= new double [nmods][npnts][2];

        mix_pars = new double [nmods];

        old_bern_ps = new double [nmods][ninps];

        old_normpars= new double [nmods][npnts][2];

        old_mix_pars = new double [nmods];

        new_bern_ps = new double [nmods][ninps];

        new_normpars= new double [nmods][npnts][2];

        new_mix_pars = new double [nmods];

        margin_density = new double [ngen][nmods];

        prob_density= new double [ngen][nmods];

        mpm_list=new int [ngen];

        /* copy arrays */

        copy_double(mix,mix_pars);

        copy_double(bernp,bern_ps);

        copy_double(npars,normpars);

       







    }

    public void refineEM(PrintStream ps) {



        ps.println("Initial parameter values for EMRefiner");

        printPars(ps,mix_pars,bern_ps,normpars);

        calc_margin();

        ps.println("Initial Marginal density");



        for ( int i=0; i<ngen ; i++ ) {

            int mpm = -1;

            double mp = -1.0;

            ps.print(" " + i + " ");

            for ( int j=0; j<nmods ; j++  ) {

                if (margin_density[i][j] > mp) {mpm=j;mp=margin_density[i][j];}

                ps.printf("%8.7f ", margin_density[i][j] );

            }

            ps.println(" :  " + mpm);

        }

                /*ps.println("Initial Probability densities (Pdens)");

                for ( int i=0; i<ngen ; i++ ) {

                        ps.print(" " + i + " ");

                        for ( int j=0; j<nmods ; j++  ) {

                                ps.printf("%9.7f ", prob_density[i][j] );

                        }

                        ps.print("\n");

                }*/





/*changes start here*/



        int niters = em_iters;

        int iter = 0;

        double abs_diff = 100;

        while ( abs_diff > 0.0 ) {

                calc_margin();

                copy_double(mix_pars,old_mix_pars);

                copy_double(bern_ps,old_bern_ps);

                copy_double(normpars,old_normpars);

                update_parameters();

                copy_double(mix_pars,new_mix_pars);

                copy_double(bern_ps,new_bern_ps);

                copy_double(normpars,new_normpars);

                ps.println("Parameter values for EMRefiner for iteration: " + iter);

                printPars(ps,mix_pars,bern_ps,normpars);



                ps.println("Marginal density for iteration: " + iter);

                for ( int i=0; i<ngen ; i++ ) {

                    int mpm = -1;

                    double mp = -1.0;

                    ps.print(" " + i + " ");

                    for ( int j=0; j<nmods ; j++  ) {

                        if (margin_density[i][j] > mp) {mpm=j;mp=margin_density[i][j];}

                            ps.printf("%8.7f ", margin_density[i][j] );

                    }

                    ps.println(" :  " + mpm);

                    mpm_list[i]=mpm;

                }

                for ( int i=0; i<ngen ; i++ ) {

                    ps.print(mpm_list[i] + " ");

                }

                ps.print("\n");



                double final_diff = 0.0;

                for ( int j=0; j<nmods ; j++ ) {

                    double temp_mix_pars = 0.0;

                    double temp_bern = 0.0;

                    double temp_mean = 0.0;

                    double temp_vars = 0.0;

                    int count = 0;

                    temp_mix_pars = temp_mix_pars + DiffDens ( old_mix_pars[j],new_mix_pars[j] );

                    for (int i=0; i<ninps ; i++ ) {

                        temp_bern=temp_bern + DiffDens ( old_bern_ps[j][i],new_bern_ps[j][i] );

                        count = count + 1;

                    }

                    for (int i=0; i<npnts ; i++ ) {

                        temp_mean=temp_mean + DiffDens ( old_normpars[j][i][0],new_normpars[j][i][0] );

                        count = count + 1;

                    }

                    for (int i=0; i<npnts ; i++ ) {

                        temp_vars=temp_vars + DiffDens ( old_normpars[j][i][1],new_normpars[j][i][1] );

                            count = count + 1;

                    }

                    final_diff= (final_diff + temp_mix_pars + temp_bern + temp_mean + temp_vars) / count;

                }

                final_diff = final_diff/nmods;

                abs_diff = final_diff;

                ps.println("Summation of absolute differences for iteration " + iter + " is " + final_diff);

                ps.print("\n");

                iter = iter +1;

                if ( iter == 5000 ) {

                    break;

                }

        }

    }



    private void update_parameters() {

        updateBernoulli();

        updateNormal();

        updateMix();

    }



    private void updateBernoulli() {



        for ( int k=0; k<nmods ; k++ ) {

            for ( int j=0 ; j<ninps ; j++ ) {

                double sum0=0.0;

                double sum1=0.0;

                for ( int i=0 ; i<ngen ; i++ ) {

                    if ( reg_inputs[i][j] == 0 ) sum0 = sum0 + margin_density[i][k];

                        else sum1 = sum1 + margin_density[i][k];

                }

                bern_ps[k][j] = sum1/(sum0+sum1);

            }

        }

    }



     private void updateNormal() {



        for ( int k=0; k<nmods ; k++ ) {

            for ( int j=0 ; j<npnts ; j++ ) {

                double mean=0.0;

                double sum=0.0;

                double var=0.0;

                for ( int i=0 ; i<ngen ; i++ ) {

                    sum = sum + margin_density[i][k];

                        mean = mean + exp_patterns[i][j]*margin_density[i][k];

                }

                mean = mean/sum;

                for ( int i=0 ; i<ngen ; i++ ) {

                    var = var + (exp_patterns[i][j]-mean)*(exp_patterns[i][j]-mean)*margin_density[i][k];

                    /*if ( var < SMALL_VAR ) var=SMALL_VAR;*/

                }

                var = var/sum;

                /*if ( var < SMALL_VAR ) var=SMALL_VAR;*/

                normpars[k][j][0] = mean;

                normpars[k][j][1] = var;

            }

        }



    }

     

    private void updateMix() {

        for (int k=0 ; k<nmods ; k++ ) {

            mix_pars[k]=0.0;

                for ( int i=0; i<ngen ; i++ ) {

                    mix_pars[k] = mix_pars[k] +margin_density[i][k];

                }

            mix_pars[k] = mix_pars[k]/ngen;

        }

    }



    /* calculate the marginal density from starting parameters */

    /* this is prob gene i generated from mix  component j calcd by Bayes rule */

    private void calc_margin( ) {

        for ( int i=0 ; i<ngen ; i++ ) {

            for ( int j=0 ; j<nmods ; j++ ) {

                double normconst = 0.0;



            for ( int k=0 ; k<nmods ; k++ ) {

                normconst = normconst + mix_pars[k]*pdens(i,k);



            }



            margin_density[i][j] = mix_pars[j]*pdens(i,j)/normconst;

            prob_density[i][j] = pdens(i,j);

            }

        }



    }





    /* calculate the probability density function */

    /* this calculates the prob density for gene gen given it is in module mdule */

    private double pdens ( int gen , int mdule )    {



        double prob=1.0;

            for (int i=0; i<ninps ; i++ ) {

                prob=prob*BernProb( reg_inputs[gen][i], bern_ps[mdule][i] );

            }

            for ( int i=0; i<npnts ; i++ ) {

                prob=prob*NormDens( exp_patterns[gen][i],normpars[mdule][i][0],normpars[mdule][i][1] );

            }

        return prob;

    }

    

    

    /*  a utility method to copy one integer array into another */

    private void copy_int(int [] from, int [] to ) {

        for (int i=0; i<from.length ; i++) to[i]=from[i];

    }

/*  a utility method to copy one double array into another */

    private void copy_double(double [] from, double [] to ) {

        for (int i=0; i<from.length ; i++) to[i]=from[i];

    }



    private void copy_double ( double [][] from, double [][] to ) {

        for (int i=0; i<from.length ; i++ ) copy_double ( from[i], to[i] );

    }



    private void copy_double ( double [][][] from, double [][][] to ) {

        for (int i=0; i<from.length ; i++ ) copy_double ( from[i], to[i] );

    }



    /* normal density */

    private double NormDens( double x, double mean, double var ) {

        if ( var < SMALL_VAR ) var=SMALL_VAR;

        double dens = Math.exp(-1.0*(x-mean)*(x-mean)/(2*var));

        dens = dens/Math.sqrt(2*Math.PI*var);

         /*dens = dens/Math.sqrt(2*Math.PI*var);*/

        return dens;

    }



    private double BernProb ( int x, double par ) {

        return par*x + (1.0-par)*(1-x);

    }



    private double DiffDens ( double curr,  double  old ) {  /********************convergence criterion*********************/

        double diff = 0.0;

        diff = diff + Math.abs(curr - old);

        return diff;

    }



    private void printPars( PrintStream ps, double [] mpars, double [][] bpars, double [][][] npars ) {

        ps.println("EM refinement with " + nmods + " modules");

        ps.println("Mixing parameters");

        for (int i=0; i<mpars.length ; i++) {

            ps.printf("%8.7f  ",mpars[i]);

        }

        ps.println(" ");

        ps.println("Bernoulli parameters ");

        for ( int i=0; i<nmods ; i++ ) {

            ps.print("Module " + i + ":");

            for (int j=0 ; j<ninps ; j++ ) {

                ps.printf("%8.7f ",bpars[i][j]);

            }

            ps.println();

        }

        ps.println("Normal parameters ");

        for ( int i=0; i<nmods ; i++ ) {

            ps.print("Module " + i + ":");

            for (int j=0 ; j<npnts ; j++ ) {

                ps.printf("(%8.7f,%8.7f) ",npars[i][j][0],npars[i][j][1] );

            }

            ps.print("\n");

        }

        ps.println();

    }

                                                                   

    







}





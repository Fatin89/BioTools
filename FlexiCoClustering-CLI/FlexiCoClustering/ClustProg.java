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
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.nio.*;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.FutureTask;
import java.util.concurrent.TimeUnit;
import java.util.Collections;
import java.util.Comparator;
import java.util.TreeSet;



public class ClustProg {

	private int [][] reg_inputs;
	private int [] modules,old_modules,best_modules;
	private double [][] exp_patterns;
	private String [] genenames;
	private int ngen,ninps,npnts;
        private int nmodules;
        private int mixcoeff;
        private int norm;
        private int EMrefine;
        private int em_iters;
	private boolean agglom=false;
	private String infoc;
	private int singlepoint=0;
	private int max_iters;
	private int out_int=0;
	private static int seed=0;
	private double msProb=0.25; /**/
	private double StartTemp=100.0;/**/
	private double tempFact=0.9; /**/
	private int MaxTemps=100; /**/
	private String modfile="Modfile.txt";
        private PrintStream outps=System.out;
	private static final double LARGE=10.0;
        private static final double SMALL=0.00001;
        private static final double SMALL_VAR=0.00001;
	private static final int repeat = 1;
	private static int nrun;
	private int rep;	
	private String output = "temp";
	 private static PrintStream ps0;
	public ClustProg(String infile, String outfile) throws IOException {
		readInput(infile);
		output = outfile;
		//System.out.println(output);
		if (outfile != null) outps = new PrintStream(new FileOutputStream(outfile,true));
	}

	private class ModScore1 { 
		int nmods;
		double score;
		ModScore1 (int n, double s ) {
			nmods=n;
			score=s;
		}
	}

	private class ModScore2 { 
		int nmods;
		double score;
		ModScore2 (int n, double s ) {
			nmods=n;
			score=s;
		}
	}

	
	
	public static void main(String [] argv) {
	
		if ( argv.length < 1) { /* if there's no argument given*/
			System.out.println("Usage java regulation.RegProg input.txt [output.txt]");
		} else {	
			String infile = argv[0];
			try {
				String outf = null;
				if (argv.length>1) {
					outf = argv[1];
					

					
				}
				ClustProg rp = new ClustProg(infile,outf);
				
				if (nrun !=  0){
					System.out.println("Running FlexiCoClustering re-submission...");
					System.out.println("Press ctrl+C to terminate" );
					System.out.println("Please INCREASE the MaxTemps to a higher value" );
					System.out.println("=================================================" );
					rp.run();
				}
				else{
					System.out.println("Running FlexiCoClustering...");
					System.out.println("Press ctrl+C to terminate" );
					System.out.println("Please change Nrun to '1' to re-submit and INCREASE the MaxTemps to a higher value" );
					System.out.println("=================================================" );
					rp.doIt();
				}
			
			} catch (IOException ioe) {
				System.out.println("Error reading input file: " + ioe.getMessage());
			}
		}

	}

	private void readInput(String file) throws IOException {

		File inf = new File(file);
		BufferedReader br = new BufferedReader(new FileReader(inf));
		String line = br.readLine();
		ngen = parseintline(line,"NItems");
		line = br.readLine();
		ninps = parseintline(line,"NBinary");
		line = br.readLine();
		npnts = parseintline(line,"NContinuous");
		line = br.readLine();
		agglom=false;
		if ( line.startsWith("Agglomerative") ) agglom=true;
		line = br.readLine();
		infoc = parseStringline(line,"IC");
		line = br.readLine();
		msProb = parsedoubleline(line,"MergeSplitProbability");
		line = br.readLine();
		max_iters = parseintline(line,"MaximumIterations");
		line = br.readLine();
		StartTemp = parsedoubleline(line,"StartTemp");
		line = br.readLine();
		tempFact = parsedoubleline(line,"TempFactor");
		line = br.readLine();
		MaxTemps= parseintline(line,"MaxTemps");
		line = br.readLine();
		rep = parseintline(line,"MaxRepIters");
		line = br.readLine();
		seed = parseintline(line,"Seed");
                line = br.readLine();
                em_iters = parseintline(line,"EMIterations");
		line = br.readLine();
		out_int = parseintline(line,"OutInterval");
		line = br.readLine();
		modfile = parseStringline(line,"ClustFile");
		line = br.readLine();
		singlepoint = parseintline(line,"SinglePoint");
                line = br.readLine();
                nmodules = parseintline(line,"NClusters");
                line = br.readLine();
                norm = parseintline(line,"NormExp");
                line = br.readLine();
                mixcoeff = parseintline(line,"MixCoefficient");
                line = br.readLine();
                nrun = parseintline(line,"Nrun");

		reg_inputs = new int [ngen][ninps];
		exp_patterns = new double [ngen][npnts];
		modules = new int [ngen];
		best_modules = new int [ngen];
		old_modules = new int [ngen];
		genenames = new String [ngen];

		StringTokenizer st=null;
		String tok=null;
		
		for (int i=0; i<ngen ; i++ ) {
			line = br.readLine();
			st = new StringTokenizer(line," ");
			genenames[i] = new String(st.nextToken());
			for ( int j=0; j<ninps ; j++ ) {
				tok = st.nextToken();
				reg_inputs[i][j] = Integer.parseInt(tok);
			}
			for ( int j=0; j<npnts ; j++) {
				tok = st.nextToken();
				exp_patterns[i][j] = Double.parseDouble(tok);
			}
		}

		br.close();

	}

	private void printRunPars( PrintStream ps ) {

		ps.println("====================================================================");
		ps.println("========================== RUNTIME PARAMETERS ======================");
		ps.println("\nNumber of Items                     : " + ngen);
		ps.println("Number of binary     inputs           : " + ninps);
		ps.println("Number of continuous values           : " + npnts);
		if (agglom) ps.println("Operation mode                        : Agglomerative");
		ps.println("Score function                        : " + infoc);
		ps.println("Merge/Split move probability          : " + msProb);
		ps.println("Maximum iterations per temperature    : " + max_iters);
		ps.println("Starting temperature                  : " + StartTemp);
		ps.println("Temperature reduction factor          : " + tempFact);
		ps.println("Maximum temperatures                  : " + MaxTemps);
		ps.println("Maximum score repeats                 : " + rep);
		ps.println("Random seed                           : " + seed);
                ps.println("EMIterations                          : " + em_iters);
		ps.println("Output interval                       : " + out_int);
		ps.println("Module output file                    : " + modfile);
		if (singlepoint>0) ps.println("Single point energy evaluation modules: " + singlepoint);
                ps.println("Number of modules                     : " + nmodules);
                ps.println("NormExp                               : " + norm);
                ps.println("MixCoefficient                        : " + mixcoeff);
                //ps.println("EMrefine                              : " + EMrefine);
		ps.println("====================================================================");
		ps.println("===================================================");

	}

	private int parseintline( String line, String id ) throws IOException {
		StringTokenizer st = new StringTokenizer(line,": ");
		String tok = st.nextToken();
		if ( !tok.equals(id) ) throw new IOException("Unexpected input: " + line);
		tok = st.nextToken();
		return Integer.parseInt(tok);
	}

	private double parsedoubleline( String line, String id ) throws IOException {
		StringTokenizer st = new StringTokenizer(line,": ");
		String tok = st.nextToken();
		if ( !tok.equals(id) ) throw new IOException("Unexpected input: " + line);
		tok = st.nextToken();
		return Double.parseDouble(tok);
	}

	private String parseStringline( String line, String id ) throws IOException {
		StringTokenizer st = new StringTokenizer(line,": ");
		String tok = st.nextToken();
		if ( !tok.equals(id) ) throw new IOException("Unexpected input: " + line);
		tok = st.nextToken();
		return tok;
	}
	
	private  ArrayList parallel2( Random rgn ){
		int count = 0;
		int [] temp_modules = new int [ngen];
		int []old_parr_modules = new int [ngen];
		copy_int(modules,temp_modules);
		/*System.out.println("~~~~~~~~~~~~~~~~~~~~");
		System.out.println(score1(temp_modules).score);
		System.out.println("~~~~~~~~~~~~~~~~~~~~");*/
		double SCORE =  0.0;
		SCORE = SCORE + score1(temp_modules).score;
		ArrayList Final_modules = new ArrayList ();
                while (count < 1) { 
			ArrayList<Integer> Best_modules = new ArrayList<Integer>(); 
			changeModule(rgn,temp_modules);
			
			
			SCORE = score1(temp_modules).score;
			for ( int j=0 ; j<ngen ; j++){
					Best_modules.add(temp_modules[j]);
			}
			Final_modules = Best_modules;
				
					
			
			/*System.out.println("####################");	
			System.out.println(score1(temp_modules).score);
			System.out.println("####################");*/	
		
			count ++;	
                }        
		/*System.out.println((Integer) Final_modules.get(0));*/
		
		
		 

		
		return Final_modules;
	}        








	private  ArrayList parallel( Random rgn , double newtemp , int mc){
					
		int count = 0;
		int [] temp_modules = new int [ngen];
		int []old_parr_modules = new int [ngen];
		copy_int(modules,temp_modules);
		/*System.out.println("~~~~~~~~~~~~~~~~~~~~");
		System.out.println(score1(temp_modules).score);
		System.out.println("~~~~~~~~~~~~~~~~~~~~");*/
		double SCORE =  0.0;
		SCORE = SCORE + score1(temp_modules).score;
		ArrayList Final_modules = new ArrayList ();
		double beta=1.0/StartTemp;
                while (count<max_iters) {  
			ArrayList<Integer> Best_modules = new ArrayList<Integer>();
			copy_int(temp_modules,old_parr_modules);
			for (int i=0 ; i< mc; i++){
				if ( rgn.nextDouble() > msProb ){
					changeModule(rgn,temp_modules);
				}
				else { 
					if ( rgn.nextDouble()>0.5 ) {
						mergeModules(rgn,temp_modules);
					}
					else {
						splitModule(rgn,temp_modules);
					}
				}
			}
			double new_SCORE = score1(temp_modules).score;
      			double diff = new_SCORE-SCORE;
			
			if ( (new_SCORE < SCORE ) || (rgn.nextDouble() < Math.exp(-beta*diff)) ){
				SCORE = new_SCORE;
				for ( int j=0 ; j<ngen ; j++){
					Best_modules.add(temp_modules[j]);
				}
				Final_modules = Best_modules;
			}
			else { 
				copy_int(old_parr_modules,temp_modules);
				SCORE = score1(temp_modules).score;
				for ( int j=0 ; j<ngen ; j++){
					Best_modules.add(temp_modules[j]);
				}
				Final_modules = Best_modules;
				
					
			}
			
			/*System.out.println("####################");	
			System.out.println(score1(temp_modules).score);
			System.out.println("####################");*/	
		
			count ++;	
                }        
		/*System.out.println((Integer) Final_modules.get(0));*/
		
		
		 

		
		return Final_modules;
	}

			
	private int [] sort( ArrayList<ArrayList> result) {
		/*System.out.println("Start sorting");*/
		int [] lowest= convertArray(result.get(0)); 
		for (int i=1; i< result.size() ;i++) {
			if ((Double) score1(convertArray(result.get(i))).score <  (Double) score1(lowest).score){
				lowest = convertArray(result.get(i));
			}
		}
		/*System.out.println((Double) score1(lowest).score);
		System.out.println("End sorting");*/
		return lowest;
	}
	


	private void doIt () {
                double startTime = System.currentTimeMillis();
		printRunPars(outps);
		int count_temp=0;
		double temp=StartTemp;
		int nExp = norm;
                if (nExp != 0) normExp();
                init_modules();
		if ( singlepoint>0 ) {
			int gpm = ngen/singlepoint;
			for ( int i=0 ; i<ngen ; i++ ) modules[i] = i/gpm;
			ModScore2 ms2 = score2(); 
			outps.println("===================================================");
			outps.println("Single point score");
			printModules(outps,ms2,modules,false);
			outps.println("===================================================");
			return;
		}

		
		Random rg1 = new Random(seed);
		ModScore2 ms2 = score2();
		outps.println("Starting modules " + ms2.nmods + " Total score " + ms2.score);
		int accept = 0;
		double ar = 0.0;
		double scr = 0.0;
		int nmods = 0;
		double oldscr = 0.0, diff = 0.0;
		int count_score = 0;
		ModScore2 oldscore = ms2;
		ModScore2 bestscore = ms2;
		double final_score = 0.0;
		int final_nmods = 0;
		int count_iter=0;
    		/*buat parallelization utk every thread has its own SA for n iters*/
		for ( count_temp=0 ; count_temp<MaxTemps ; count_temp++ ) {
			count_iter=count_iter+1;
			final double newtemp = StartTemp*tempFact;
			/*if (count_iter/out_int ==1) {
			outps.println("===================================================");
			outps.println("Starting temperature: " + temp);
			outps.println("===================================================");
			}*/
			int processors = 5;//Runtime.getRuntime().availableProcessors();
			ExecutorService executor = Executors.newFixedThreadPool(processors);
			ArrayList<Callable<ArrayList>> callList = new ArrayList<Callable<ArrayList>> ();
			List<Future<ArrayList>> future = new ArrayList<Future<ArrayList>> ();
			for (int i=1; i<processors+1;i ++){
				
				final int p = 2;
				if (i < processors/2){
					final int k = 1;
					callList.add(new Callable<ArrayList> (){
					
						public  ArrayList call() throws InterruptedException, ExecutionException {
							Random rgn = new Random(seed);
							int c = k;
							/*System.out.println(parallel(rgn));*/
        						return parallel(rgn,newtemp,c);
						}
					});
				}
				else {
					if (i> (processors/2)+2){
						final int k = 2;
						callList.add(new Callable<ArrayList> (){
						
							public  ArrayList call() throws InterruptedException, ExecutionException {
								Random rgn = new Random(seed);
								int c = k;
								/*System.out.println(parallel(rgn));*/
        							return parallel(rgn,newtemp,k);
							}
						});
					}    		
					else {
						final int k = 3;
						callList.add(new Callable<ArrayList> (){
						
							public  ArrayList call() throws InterruptedException, ExecutionException {
								Random rgn = new Random(seed);
								int c = k;
								/*System.out.println(parallel(rgn));*/
        							return parallel(rgn,newtemp,k);
							}
						});
					}    	
				}
									
				
			
			}
			double new_score = bestscore.score;
			ArrayList<ArrayList> results = new ArrayList<ArrayList>();
			try {
    				future=executor.invokeAll(callList);
    				executor.shutdown();
    				/*executor.awaitTermination(20, TimeUnit.SECONDS);*/
				for ( int i=0; i < future.size(); i++){
					results.add(future.get(i).get());
					
				}
								
				
					
			}
			catch (InterruptedException e){
				e.printStackTrace();
				System.err.println("tasks interrupted");
			}
			catch (ExecutionException e){
  				System.err.println("tasks interrupted");
			}

			
			int [] SortedResults = sort(results);
			for (int i=0;i<ngen;i++){
				modules[i]= SortedResults[i];
			}
	
			
			
			oldscore = ms2;
			ms2 = score2();
			scr = ms2.score;
			if (scr < oldscore.score){
				bestscore = new ModScore2(ms2.nmods,scr);
				final_score = bestscore.score;
				final_nmods = bestscore.nmods;
				copy_int(modules,best_modules);
			}
			else {
				ms2 = oldscore;
				copy_int(old_modules,modules);
			}
			StartTemp=newtemp;
		
			
			if (bestscore.score == new_score){
				count_score = count_score + 1;
			}
			else{
				count_score = 1;
			}
  


			ms2 = bestscore;
			scr = ms2.score;
			copy_int(best_modules, modules);



                       if (count_iter/out_int ==1) {
                                outps.println("===================================================");
                                outps.println("Iteration: " + count_temp);
                                printModules(outps,ms2,modules,false);
                                outps.println("===================================================");
                       }
  



                       if (count_iter/out_int ==1) {
                                outps.println("Current best score: " + bestscore.score);
                                outps.println("Current best modules: " + bestscore.nmods);
                                outps.println("Current best score number of hits found: " + count_score);
                                double endTime   = System.currentTimeMillis();
                                double totalTime = endTime - startTime;
                                outps.println("Current time taken: " + (totalTime/1000.0)/60.0 + "minutes");
                                outps.println("Current temperature step: " + count_temp);
                                outps.println("===================================================");

                                outps.println("===================================================");
                                outps.println("Finishing temperature: " + newtemp);
                                ar = (double) accept/ (double) max_iters;
                                outps.println("Acceptance ratio: : " + ar);
                                outps.println("===================================================");
                                count_iter=0;
                        }
                        




			if (count_score == rep) {
				break;
			}
		}
		outps.println("\n\nBest solution found");
		printModules(outps,bestscore,best_modules,false);
		try {
			PrintStream mps = new PrintStream(new FileOutputStream(modfile));
			mps.println("\n\nBest solution found");
			printModules(mps,bestscore,best_modules,true);
			mps.close();
		}
		catch (FileNotFoundException fnf ) {
			outps.println("Error: Could not open module output file");
		}
                /* now refine the answer with Expectation Maximization */
                /* variables holding parameters and marginal density */
              	
               
               	int emr = em_iters;
		if ( emr!=0 ) {
                	
			int nmodules=countModules(best_modules);
                	double [][] bernpars = new double [nmodules][ninps];
                	double [][][] normpars= new double [nmodules][npnts][2];
                	double [] mixpars = new double [nmodules];

                	/* initialise the distribution parameters */

			init_parameters (best_modules,mixpars,bernpars,normpars);
			try{ps0 = new PrintStream(new BufferedOutputStream(new FileOutputStream("EMRefinement.txt")),true);
                	EMRefiner emrf = new EMRefiner(exp_patterns,reg_inputs,ngen,ninps,npnts,nmodules,mixpars,bernpars,normpars,em_iters);
                	
			emrf.refineEM(ps0);
			}
			catch ( FileNotFoundException fnf ){
				System.out.println("EMrefinement cannot proceed");
			}
			
		}
                double endTime   = System.currentTimeMillis();
                double totalTime = endTime - startTime;
		outps.println("Total time: " + (totalTime/1000.0)/60.0 + "minutes");
		return;
	}

	private void run() throws IOException,FileNotFoundException {
		int nExp = norm;
                if (nExp != 0) normExp();
		double startTime = System.currentTimeMillis();
		FileReader fr = new FileReader(output);
		BufferedReader br = new BufferedReader(fr);
		StringTokenizer st=null;
		String tok=null;
		String line =null;
		List<String> tmp = new ArrayList<String>();
		int count =0;
		do {
			line = br.readLine();
			tmp.add(line);
		}	
		while (line!=null);
			        
		List<String> temp_modules = new ArrayList<String>();
		
		for (int i=tmp.size()-1; i>=0 ; i--) {
			if (tmp.get(i) != null){
				if (tmp.get(i).length() == 7){
					count ++;
				}
			}
			if (count == 3){
				temp_modules.add(tmp.get(i));
			}
		}		
		String[] temp1 = temp_modules.get(16).split(" ");
		int [] last_modules = new int[temp1.length];
		for(int i = 0; i < temp1.length; i++) {
   			last_modules[i] = Integer.parseInt(temp1[i]);
		}
		
		copy_int(last_modules,modules);
		copy_int(modules,best_modules);
		copy_int(modules,old_modules);

		/*System.out.println(temp_modules.get(0));
		System.out.println(temp_modules.get(1));
		System.out.println(temp_modules.get(2));	
		System.out.println(temp_modules.get(3));
		System.out.println(temp_modules.get(4));
		System.out.println(temp_modules.get(5));
		System.out.println(temp_modules.get(6));
		System.out.println(temp_modules.get(7));
		System.out.println(temp_modules.get(8));
		System.out.println(temp_modules.get(9));
		System.out.println(temp_modules.get(10));	
		System.out.println(temp_modules.get(11));
		System.out.println(temp_modules.get(12));
		System.out.println(temp_modules.get(13));
		System.out.println(temp_modules.get(14));
		System.out.println(temp_modules.get(15));
		System.out.println(temp_modules.get(16));	*/
		
		double time_difference=Double.parseDouble(temp_modules.get(11).split(":")[1].split("m")[0]);
		double StartTemp=Double.parseDouble(temp_modules.get(7).split(":")[1].split(" ")[1]);
		int count_score=Integer.parseInt(temp_modules.get(12).split(":")[1].split(" ")[1]);
		/*System.out.println(temp);
		System.out.println(count_temp);*/

		Random rg1 = new Random(seed);
		ModScore2 ms2 = score2();
		outps.append("Starting modules " + ms2.nmods + " Total score " + ms2.score);



		int accept = 0;
		double ar = 0.0;
		double scr = 0.0;
		int nmods = 0;
		double oldscr = 0.0, diff = 0.0;
		
		ModScore2 oldscore = ms2;
		ModScore2 bestscore = ms2;
		double final_score = 0.0;
		int final_nmods = 0;
		int count_iter=0;
    		/*buat parallelization utk every thread has its own SA for n iters*/
		for ( int count_temp=Integer.parseInt(temp_modules.get(10).split(":")[1].split(" ")[1]);count_temp<MaxTemps ; count_temp++ ) {
			count_iter=count_iter+1;
			final double newtemp = StartTemp*tempFact;
			/*if (count_iter/out_int ==1) {
			outps.println("===================================================");
			outps.println("Starting temperature: " + temp);
			outps.println("===================================================");
			}*/
			int processors = 4; //Runtime.getRuntime().availableProcessors();
			ExecutorService executor = Executors.newFixedThreadPool(processors);
			ArrayList<Callable<ArrayList>> callList = new ArrayList<Callable<ArrayList>> ();
			List<Future<ArrayList>> future = new ArrayList<Future<ArrayList>> ();
			for (int i=1; i<processors+1;i ++){
				
				final int p = 2;
				if (i < processors/2){
					final int k = 1;
					callList.add(new Callable<ArrayList> (){
					
						public  ArrayList call() throws InterruptedException, ExecutionException {
							Random rgn = new Random(seed);
							int c = k;
							/*System.out.println(parallel(rgn));*/
        						return parallel(rgn,newtemp,c);
						}
					});
				}
				else {
					if (i> (processors/2)+2){
						final int k = 2;
						callList.add(new Callable<ArrayList> (){
						
							public  ArrayList call() throws InterruptedException, ExecutionException {
								Random rgn = new Random(seed);
								int c = k;
								/*System.out.println(parallel(rgn));*/
        							return parallel(rgn,newtemp,k);
							}
						});
					}    		
					else {
						final int k = 3;
						callList.add(new Callable<ArrayList> (){
						
							public  ArrayList call() throws InterruptedException, ExecutionException {
								Random rgn = new Random(seed);
								int c = k;
								/*System.out.println(parallel(rgn));*/
        							return parallel(rgn,newtemp,k);
							}
						});
					}    	
				}
									
				
			
			}
			double new_score = bestscore.score;
			ArrayList<ArrayList> results = new ArrayList<ArrayList>();
			try {
    				future=executor.invokeAll(callList);
    				executor.shutdown();
    				/*executor.awaitTermination(20, TimeUnit.SECONDS);*/
				for ( int i=0; i < future.size(); i++){
					results.add(future.get(i).get());
					
				}
								
				
					
			}
			catch (InterruptedException e){
				e.printStackTrace();
				System.err.println("tasks interrupted");
			}
			catch (ExecutionException e){
  				System.err.println("tasks interrupted");
			}

			
			int [] SortedResults = sort(results);
			for (int i=0;i<ngen;i++){
				modules[i]= SortedResults[i];
			}
	
			
			
			oldscore = ms2;
			ms2 = score2();
			scr = ms2.score;
			if (scr < oldscore.score){
				bestscore = new ModScore2(ms2.nmods,scr);
				final_score = bestscore.score;
				final_nmods = bestscore.nmods;
				copy_int(modules,best_modules);
			}
			else {
				ms2 = oldscore;
				copy_int(old_modules,modules);
			}
			StartTemp=newtemp;
		
			
			if (bestscore.score == new_score){
				count_score = count_score + 1;
			}
			else{
				count_score = 1;
			}
  


			ms2 = bestscore;
			scr = ms2.score;
			copy_int(best_modules, modules);



                       if (count_iter/out_int ==1) {
                                outps.append('\n'+"===================================================");
                                outps.append('\n'+"Iteration: " + count_temp+'\n');
                                printModules(outps,ms2,modules,false);
                                outps.append("===================================================");
                       }
  



                       if (count_iter/out_int ==1) {
                                outps.append('\n'+"Current best score: " + bestscore.score);
                                outps.append('\n'+"Current best modules: " + bestscore.nmods);
                                outps.append('\n'+"Current best score number of hits found: " + count_score);
                                double endTime   = System.currentTimeMillis();
                                double totalTime = endTime - startTime  ;
				double [] TotalTime = new double[2];
				TotalTime[0]=(totalTime/1000.0)/60.0;
				TotalTime[1]=time_difference;
				double sum = 0.0;
				for (double i : TotalTime)
    					sum += i;
                                outps.append('\n'+"Current time taken: " + sum + "minutes");
				//System.out.println(time_difference);
				//System.out.println(time_difference+1);
                                outps.append('\n'+"Current temperature step: " + count_temp);
                                outps.append('\n'+"===================================================");

                                outps.append('\n'+"===================================================");
                                outps.append('\n'+"Finishing temperature: " + newtemp);
                                ar = (double) accept/ (double) max_iters;
                                outps.append('\n'+"Acceptance ratio: : " + ar);
                                outps.append('\n'+"===================================================");
                                count_iter=0;
                        }
                        




			if (count_score == rep) {
				break;
			}
		}
		outps.append("\n\nBest solution found");
		printModules(outps,bestscore,best_modules,false);
		try {
			PrintStream mps = new PrintStream(new FileOutputStream(modfile));
			mps.println("\n\nBest solution found");
			printModules(mps,bestscore,best_modules,true);
			mps.close();
		}
		catch (FileNotFoundException fnf ) {
			outps.append("Error: Could not open module output file");
		}
                /* now refine the answer with Expectation Maximization */
                /* variables holding parameters and marginal density */
              	
               
               	int emr = em_iters;
		if ( emr!=0 ) {
                	
			int nmodules=countModules(best_modules);
                	double [][] bernpars = new double [nmodules][ninps];
                	double [][][] normpars= new double [nmodules][npnts][2];
                	double [] mixpars = new double [nmodules];

                	/* initialise the distribution parameters */
			init_parameters (best_modules,mixpars,bernpars,normpars);
                        try{ps0 = new PrintStream(new BufferedOutputStream(new FileOutputStream("EMRefinement.txt")),true);
                        EMRefiner emrf = new EMRefiner(exp_patterns,reg_inputs,ngen,ninps,npnts,nmodules,mixpars,bernpars,normpars,em_iters);

                        emrf.refineEM(ps0);
                        }
                        catch ( FileNotFoundException fnf ){
                                System.out.println("EMrefinement cannot proceed");
                        }
         

		}
                double endTime   = System.currentTimeMillis();
                double totalTime = endTime - startTime;
		double [] TotalTime = new double[2];
		TotalTime[0]=(totalTime/1000.0)/60.0;
		TotalTime[1]=time_difference;
		double sum = 0.0;
		for (double i : TotalTime)
    			sum += i;
		outps.println("Total time: " + sum + "minutes");
		
		return;


		
		
	}
		
		
		
	private int [] convertArray ( ArrayList vals ) {
		int [] listVal = new int[vals.size()];
		for ( int i=0 ; i < vals.size() ; i++ ) {
			listVal[i] = (Integer) vals.get(i);
		}
		return listVal;
	}
				

	private void changeModule (Random rgns, int [] temp_module) {
		int rndgene = rgns.nextInt(ngen);
		int rndmod =  rgns.nextInt(ngen);
		temp_module[rndgene] = rndmod;
	}

	private void mergeModules (Random rgns, int [] temp_module) {
		int rndgene1 = rgns.nextInt(ngen);
		int rndgene2 = rgns.nextInt(ngen);
		int mod1 = temp_module[rndgene1];
		int mod2 = temp_module[rndgene2];
		for (int i=0; i<ngen ; i++ ) {
			if ( temp_module[i] == mod2 ) temp_module[i]=mod1;
		}
	} 

	private void splitModule (Random rgns, int [] temp_module) {
		int rndgene = rgns.nextInt(ngen);
		int mod = temp_module[rndgene];
		int newmod = rgns.nextInt(ngen);
		for (int i=0; i<ngen ; i++ ) {
			if ( (temp_module[i] == mod) && rgns.nextDouble()>0.5 ) temp_module[i]=newmod;
		}
	}
 
	private void printModules(PrintStream ps, ModScore2 ms2, int [] mods, boolean full) {

		
		if (!full){
			ps.append("Solution with " + ms2.nmods + " modules");
			ps.append('\n'+"Total score " + ms2.score);
			ps.append("\nModules"+'\n');			
			for ( int i=0; i<mods.length ; i++) ps.append(mods[i] + " ");
			ps.append("\n");
		}
		if (full) {
			ps.println("Solution with " + ms2.nmods + " modules");
			ps.println("Total score " + ms2.score);
			ps.println("\nModules");
			for (int i=0, count=0; i<ngen ; i++ ) {
				int [] mod = indices(mods,i);
				if (mod[0] == -1 ) continue;
				count++;
				double [] modr = calcModReg(mod);
				double [][] mode = calcModExp(mod);
				ps.println("Module " + count);
				ps.print("Genes: ");
				for ( int j=0; mod[j]!=-1 ; j++ ) ps.print(genenames[mod[j]] + "\t");
				ps.println("\nRegulation");
				for ( int j=0; j<ninps ; j++ ) ps.printf("%8.7f\t",modr[j]);
				ps.print("\nExpression\n");
				for ( int j=0; j<npnts ; j++ ) ps.printf("%12.9f\t",mode[j][0]);
				ps.print("\n");
				for ( int j=0; j<npnts ; j++ ) ps.printf("%12.9f\t",mode[j][1]);
				ps.print("\nGene patterns\n");
				for ( int j=0; mod[j]!=-1 ; j++ ) {
					ps.printf("%-10s : ",genenames[mod[j]]);
					for (int k=0; k<ninps ; k++) ps.print(reg_inputs[mod[j]][k] + " ");
					ps.print(": ");
					for (int k=0; k<npnts ; k++) ps.printf("%12.9f ",exp_patterns[mod[j]][k]);
					ps.print("\n");
				}
				ps.print("\n\n");
			}
		}
	}



	private void printModules1(PrintStream ps, ModScore1 ms1, int [] mods, boolean full) {

		
		if (!full){
			ps.append('\n'+"Solution with " + ms1.nmods + " modules");
			ps.append('\n'+"Total score " + ms1.score);
			ps.append("\nModules"+'\n');			
			for ( int i=0; i<mods.length ; i++) ps.append(mods[i] + " ");
			ps.append("\n");
		}
		if (full) {
			ps.println("Solution with " + ms1.nmods + " modules");
			ps.println("Total score " + ms1.score);
			ps.println("\nModules");
			for (int i=0, count=0; i<ngen ; i++ ) {
				int [] mod = indices(mods,i);
				if (mod[0] == -1 ) continue;
				count++;
				double [] modr = calcModReg(mod);
				double [][] mode = calcModExp(mod);
				ps.println("Module " + count);
				ps.print("Genes: ");
				for ( int j=0; mod[j]!=-1 ; j++ ) ps.print(genenames[mod[j]] + "\t");
				ps.println("\nRegulation");
				for ( int j=0; j<ninps ; j++ ) ps.printf("%8.7f\t",modr[j]);
				ps.print("\nExpression\n");
				for ( int j=0; j<npnts ; j++ ) ps.printf("%12.9f\t",mode[j][0]);
				ps.print("\n");
				for ( int j=0; j<npnts ; j++ ) ps.printf("%12.9f\t",mode[j][1]);
				ps.print("\nGene patterns\n");
				for ( int j=0; mod[j]!=-1 ; j++ ) {
					ps.printf("%-10s : ",genenames[mod[j]]);
					for (int k=0; k<ninps ; k++) ps.print(reg_inputs[mod[j]][k] + " ");
					ps.print(": ");
					for (int k=0; k<npnts ; k++) ps.printf("%12.9f ",exp_patterns[mod[j]][k]);
					ps.print("\n");
				}
				ps.print("\n\n");
			}
		}
	}





	private void copy_int(int [] from, int [] to ) {
		for (int i=0; i<from.length ; i++) to[i]=from[i];
	}
         
        private void copy_double(double [] from, double [] to ) {
		for (int i=0; i<from.length ; i++) to[i]=from[i];
	}

	private void normExp() {
		double mn, sd;
		for ( int i=0; i<ngen ; i++) {
			mn = mean(exp_patterns[i], npnts);	
			sd = Math.sqrt(var(exp_patterns[i], npnts));	
			for ( int j=0; j<npnts ; j++) {
				exp_patterns[i][j] = (exp_patterns[i][j]-mn)/sd;
			}
		}
	}

	private void init_modules() {
              
                int nmods = nmodules;
                if (nmods != 0){
		        if (agglom) for ( int i=0, mod=-1 ; i<ngen ; i++ ){ 
                                int gpm = ngen/nmods;
				if ( (i % gpm) == 0 ) mod++;
				modules[i]=mod;
                        }
                        else{ 
                                for (int i=0, mod=-1 ; i<ngen; i++){
                                int gpm = ngen/nmods;
				if ( (i % gpm) == 0 ) mod++;
				modules[i]=mod;
                                }
                        }
                }
		
                else{
                        if (agglom) for ( int i=0 ; i<ngen ; i++ ) modules[i] = i;
			else for (int i=0; i<ngen ; i++)modules[i] = 0;
                        
                }	
	        
                    
		
                
                copy_int(modules,best_modules);
		copy_int(modules,old_modules);
        }

	private int modSize ( int [] module ) {
		int count = 0;
		for ( count = 0; module[count] != -1 ; count++ );
		return count;
	} 

        private int countModules( int [] mods ) {
		int count=0;
		for ( int i=0; i<ngen ; i++ ) {
			int [] mod = indices(mods,i);
			if ( mod[0] != -1 ) count ++;
		}	
		return count;
	}
      
 

	/* this is the AIC/BIC score */
	private ModScore1 score1 (int [] module) {
		int count=0;
		double s=0.0,sm=0.0,alpha=0.0;
                int mf = mixcoeff;
    		
                if (mf!=0) {
                	for (int i=0; i<ngen ; i++) {
				int [] mod = indices(module,i);
				if ( mod[0] != -1 ) {
					/* score the module, increment the total score and the number of modules */
					sm = score_module(mod);
					alpha = (double) modSize(mod)/ (double) ngen;
					sm = sm + Math.log(alpha);
					s = s + sm; 
					count++;
				}
			}
                	int npars = count*(ninps+2*npnts) + count - 1; 
			if (infoc.equals("AIC"))s = (-2*s) + (2.0*npars);
        		if (infoc.equals("AIC2.5"))s = (-2*s) + (2.0*npars);
        		if (infoc.equals("AIC3"))s = (-2*s) + (3.0*npars);
        		if (infoc.equals("AIC4"))s = (-2*s) + (4.0*npars);
       			if (infoc.equals("AIC5"))s = (-2*s) + (5.0*npars);
       			if (infoc.equals("HQC"))s = (-2*s) + (2.0*npars*Math.log(Math.log(ngen)));
        		if (infoc.equals("BIC"))s = npars*Math.log(ngen)-2*s;
       			if (infoc.equals("CAIC"))s = (-2*s) + (npars*(Math.log(ngen)+1));
		}
		
		return new ModScore1(count,s);
	}


	private ModScore2 score2 () {
		int count=0;
		double s=0.0,sm=0.0,alpha=0.0;
                int mf = mixcoeff;
    		
                if (mf!=0) {
                	for (int i=0; i<ngen ; i++) {
				int [] mod = indices(modules,i);
				if ( mod[0] != -1 ) {
					/* score the module, increment the total score and the number of modules */
					sm = score_module(mod);
					alpha = (double) modSize(mod)/ (double) ngen;
					sm = sm + Math.log(alpha);
					s = s + sm; 
					count++;
				}
			}
                	int npars = count*(ninps+2*npnts) + count - 1; 
			if (infoc.equals("AIC"))s = (-2*s) + (2.0*npars);
        		if (infoc.equals("AIC2.5"))s = (-2*s) + (2.0*npars);
        		if (infoc.equals("AIC3"))s = (-2*s) + (3.0*npars);
        		if (infoc.equals("AIC4"))s = (-2*s) + (4.0*npars);
       			if (infoc.equals("AIC5"))s = (-2*s) + (5.0*npars);
       			if (infoc.equals("HQC"))s = (-2*s) + (2.0*npars*Math.log(Math.log(ngen)));
        		if (infoc.equals("BIC"))s = npars*Math.log(ngen)-2*s;
       			if (infoc.equals("CAIC"))s = (-2*s) + (npars*(Math.log(ngen)+1));
		}
		
		return new ModScore2(count,s);
	}







	private double score_module(int [] module) {
		double s=0.0;
		s = score_mod_exp(module);
		s =  s + score_mod_reg(module);
		
		return s;
	}

	private double score_mod_exp(int [] module) {
		double [][] consexp = calcModExp(module);
		int gene;
		double score=0.0,s=0.0;
		/* assume normal distribution about consensus */ 
		for (int i=0; module[i] != -1 ; i++) {
			gene=module[i];
			for ( int j=0; j<npnts ; j++) {
				s = log_normdensity(exp_patterns[gene][j],consexp[j][0],consexp[j][1]);
				score += s;
			}
		}
		return score;
	}




	private double score_mod_reg(int [] module) {
		int gene;
		double [] modreg = calcModReg(module);
		double s;
		/* use match and mismatch probabilities */
		double score=0.0;
		for (int i=0; module[i] != -1 ; i++) {
			gene=module[i];
                       	for ( int j=0; j<ninps ; j++) {
				if ( reg_inputs[gene][j] == 1 ) s = Math.log(modreg[j]);
				else s = Math.log(1-modreg[j]);
				score += s;
				
			}
			
		}
		return score;
	}


	double [][] calcModExp(int [] module) {
		/* just calculate average at each time point */
		int gene;
		double [] exp = new double [ngen];
		double [][] modexp = new double [npnts][2];
		for ( int j=0 ; j<npnts ; j++ ) {
			int count=0;
			for ( int i=0; module[i] != -1 ; i++ ) {
				gene=module[i];
				exp[i] = exp_patterns[gene][j]; 
				count++;
			}
			modexp[j][0]=mean(exp,count);
			modexp[j][1]=var(exp,count);
		}
		return modexp;
	}

	double [] calcModReg(int [] module) {
		double [] regp = new double [ninps];
		int gene;
		for (int i=0; i<regp.length ; i++) regp[i]=0.0;
		for ( int j=0 ; j<ninps ; j++ ) {
			int i=0,count=0;
			for ( i=0; module[i] != -1 ; i++ ) {
				gene=module[i];
				if ( reg_inputs[gene][j] == 1 ) count++;
			}
			regp[j] = (double) count/(double) i;
			if (regp[j] < SMALL) regp[j] = SMALL;
			if (regp[j] > (1.0-SMALL) ) regp[j] = 1.0-SMALL;
		}
		return regp;
	}
        
        private void init_parameters ( int [] mods, double [] mixpars, double [][] bernpars, double [][][] normpars ) {
                int count=0;
                double [] work;
                double [][] work1;
                for (int i=0; i<ngen ; i++ ) {
                        int [] mod = indices(mods,i);
                        if ( mod[0] == -1 ) continue;
                        count++;
                        int size = modSize(mod);
                        mixpars[count-1] = (double) size / (double) (ngen);
                        work = calcModReg(mod);
                        copy_double(work,bernpars[count-1]);
                        work1 = calcModExp(mod);
                        for ( int j=0 ; j<npnts ; j++ )
                                for ( int k=0; k<2 ; k++ ) normpars[count-1][j][k] = work1[j][k];
                }
        }



	private int [] indices ( int [] arr, int i) {
		int k=0,j=0;
		int len = arr.length;
		int [] ind = new int [len+1];
		for ( j=0 ; j<=len ; j++) ind[j]=-1;
		for ( j=0 ; j<len ; j++ ) {
			if ( arr[j] == i ) {
				k++;
				ind[k-1] = j;
			}
		}
		return ind;
	}
	

	private double mean( double [] x, int size) {
		double mean=0.0;
		for ( int i=0; i<size ; i++) mean += x[i];
		mean = mean/size;
		return mean;
	}

	private double var( double [] x, int size ) {
		if (size<2) return LARGE;
		double m=mean(x,size);
		double s=0;
		for (int i=0; i<size ; i++ ) s += (x[i]-m)*(x[i]-m);
		s = s/(double) (size-1);
		if ( s < SMALL_VAR ) s=SMALL_VAR;
		return s;
	}

	private double log_normdensity (double x, double mean, double var ) {
		double logdens;
		logdens = (-0.5*Math.log(2*Math.PI*var)) - ((x-mean)*(x-mean)/(2*var));
		return logdens; 
	}


}


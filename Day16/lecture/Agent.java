/**
 * @author 
 * @see Copyright 2007
 *
 */
package Belief2;


/**
 * Import statements: import libraries that aid in the development of this class
 * java. libraries are provided by java.
 * uchicago. libraries are part of the RePast simulation toolkit
 */
import java.awt.Color;
import java.util.Vector;
import java.util.ArrayList;

import uchicago.src.sim.gui.Drawable;
import uchicago.src.sim.gui.SimGraphics;
import uchicago.src.sim.space.Object2DGrid;
import uchicago.src.sim.util.Random;

public class Agent implements Drawable {

  	private Color agentShade;							// unused now

	// some class variables

	private static int nextID = 0;					   	// for generating unique ID's for agents.
  	private static Model model;                       	// model object
	private static Object2DGrid world;  			 	// agent space (grid)
  	private static Resource resource; 					// forest resource object
	private static Institution institution;				// the institution
  	private static Vector<Agent> hhAgentList;         	// list of all other household agents in the model

	private static double extractionLevelMax = 3.0;     // multiplier used in determining Cmax (constraint on agent's extraction level)

	/**
	 * @param scaledCurrentExtractionLevel the scaledCurrentExtractionLevel to set
	 */
	public void setScaledCurrentExtractionLevel(double scaledCurrentExtractionLevel) {
		this.scaledCurrentExtractionLevel = scaledCurrentExtractionLevel;
	}

	private int ID;										// an id number for the agent
  	private int x, y;									// grid location 0,0 is top left		

  	private int searchType;							  
  	private int numRandomOthers;					    // number random others it picks to get a norm (searchType=0)

	private Vector<Agent> neighborList;					// all of the agent's 'neighbors'

  	private double hhSize;								// the size of the household, the number of occupants
	private double leisureEndowment;                    // the number of hours of leisure each household is endowed with. This serves as a source of labor for gathering fuelwood

  	private double Cmax;								// the highest extraction level permitted to the household (absolute terms)
  	private double subsistenceRequirement;				// The extraction level beyond which marginal benefits of consumption become lower
  	private double scaledSubsistenceRequirement;		// A scaled version of the subsistenceRequirement
  	private double decreasingReturnsFactor;				// < 1.0 means diminishing returns to consumption beyond subsistenceRequirement
  														// =1.0 means constant returns even beyond subsistenceRequirement
	private double absAmountTryToExtract;				// amount the agent wants to extract (absolute terms)
	private double absAmountExtracted;					// amount agent does extract (absolute terms)
  	private double scaledCurrentExtractionLevel;		// fraction of Cmax consumed this step
  	private double scaledPreviousExtractionLevel;		// fraction of Cmax consumed last step
  	private double scaledSocialNormExtractionLevel;    	// fraction of "their" Cmax that neighbors extracted
  	private double scaledSustainabilityLevel;			//  fraction of Cmax that institution tells agent to extract
  	private double absSustainabilityLevel;				//  amount that institutution tells agent to extract
  	
  	private double socEnvRatio;
  	private double alphaConsumption;
  	private double alphaLeisure;
	private double alphaNormSust;
	private int optimizer;								// when this variable = 1 the agent optimizes, when it is 0 the agent optimizes over a set of extraction levels	private 
	int cognitiveAbility;							// no of extraction level values that make up the bounded set which the agent maximizes over if the "optimizer=0"
	private int normalSearch;							// =0 implies sample extractino levels chosen from uniform dist, =1 implies normal dist N(scaledPreviousExtractionLevel, normalSearchStd)
	private double normalSearchStd;						// Std Dev of normal distribution if normalSearch = 1	
	private int useBurnInPeriod;						// =1 implies uses uniform dist for sampling extractino levels for first 20 periods, then uses normal dist
														// but only if normalSearch = 1 
	
	private int timeStepsAsLowExtractor;                    // Number of time-steps spent as a low extractor  
	private int timeStepsAsMedExtractor;                    // Number of timesteps spent as a Medium extractor
	private int timeStepsAsHighExtractor;                   // Number of timesteps spent as a HIgh extractor
	
	private double scaledTotalAmountExtracted;					// Total scaled amount extracted ove an entire run
	
	private int clusterID;								 // =1 or 2 
	
	private int timeClock;
	

	// effectively used as alias' for terms in optimization equation	
	private	double Ans;                                     // weight for the norm vs institutional sustainability (alphaNormSust)
	private double G;                                       // the fuelwood gathering speed for the household
	private double L;                                       // the total leisure time for the household
	private double Al;                                      // the weight placed on leisure by the household
	private double Aser;                                    // the ratio of of 
	private double N;                                       // socialNormExtractionLevel
	private double S;                                       // sustainability level
	private double Ac;                                      // weight the household places on consumption (alphaConsumption)

		
	
	/**
	 * Constructor: Initializes the agent location and space it will be located in
	 * @param xValue x axis location
	 * @param yValue y axis location
	 * @param spc	Space (grid) the agent is located in 
	 */
  	public Agent (int xValue, int yValue, int opt, int cogAbil, int normSearch, double normSearchStd, double sEnvRatio, double aNormSust, 
  			double alphaCons, double alphaLeis, int sType, int householdSize, double prevScaledExtLevel, int cID, double  decRetFact, int useBurnIn){

		ID = nextID++;

    	x = xValue;
    	y = yValue;
    	agentShade = new Color(1,0,0);
    	clusterID=cID;
    	
    	timeClock = 0;

    	numRandomOthers = model.getNumRandomOthers();
    	searchType = sType;

    	hhSize = householdSize;
//		leisureEndowment = 8*30 + (hhSize-2) * 4 * 30;      // 8 hrs per day labor available from adult women, half that number for children (assuming 50% relative efficicency of work). Adult males assumed not to collect fuelwood
		leisureEndowment = 5*30;   
		optimizer = opt;
    	cognitiveAbility = cogAbil;
    	normalSearch = normSearch;
    	normalSearchStd = normSearchStd;
    	useBurnInPeriod = useBurnIn;
    	
    	socEnvRatio = sEnvRatio;
    	alphaNormSust = aNormSust;
    	alphaConsumption = alphaCons;
    	alphaLeisure = alphaLeis;	
    	
    	decreasingReturnsFactor = decRetFact;
    	
    	subsistenceRequirement = (hhSize*240.0/16.0)/600.0; // this will be in m3 per month
   		// =   (hhSize persons) * (240 MJ/person-month/16 MJ/kg)/(600 kg/m3) 
    	    	
    	Cmax = extractionLevelMax*subsistenceRequirement;		// The upper bound for consumption
    	scaledSubsistenceRequirement = 1.0/extractionLevelMax;  // A scaled consumption level above scaledSubsistenceRequirement
    															// uses more than subsistenceRequirement m3 / month of wood
   		
       	scaledPreviousExtractionLevel = prevScaledExtLevel;
    	scaledCurrentExtractionLevel = prevScaledExtLevel;
    	scaledSustainabilityLevel = 0;
    	

		scaledTotalAmountExtracted = 0.0;

		timeStepsAsLowExtractor = 0;
		timeStepsAsMedExtractor = 0;
		timeStepsAsHighExtractor = 0;
	}

	/////////////////////////////////////////////////////////////////////////////////////
	// setUpNeighborList
	// set up agent's list of neighbors it listens to for advice, based on searchType.
	// This list is fixed for the entire run.
	// searchType=1:  
    //    setup Moore neighbors, then each has probability of
    //    being replaced by random other of probPickRandomNeighbor (pPRN).
	//    Thus pPRN=0 is pure Moore, pPRN=1 is pure random, but fixed!
	//    
	// NB: this assumed the world is full of agents.
	//     assume not a torus, so agents a edges/corners will have < 4.
	//
	public void setUpNeighborList() {

		neighborList = (Vector<Agent>) world.getMooreNeighbors( x, y, false );
		double probRewire = model.getProbPickRandomNeighbor();

		if ( model.getRDebug() > 1 ) 
			System.out.printf( "-> setupNeighborList for %d (%d,%d): %d neighbors, probRewire=%.3f.\n", 
							   ID, x, y,  neighborList.size(), probRewire );

		// now rewire each neighbor as needed
		for ( int i = 0; i < neighborList.size(); ++ i ) {
			if ( Model.getUniformDoubleFromTo( 0.0, 1.0 ) < probRewire ) {
				// rewiring this one, so pick a random other (not me!)
				Agent other = pickRandomOtherAgent();
				if (  model.getRDebug() > 2 ) {
					Agent oldNbor = neighborList.get(i);
					System.out.printf( "    replace %d (%d,%d) with %d (%d,%d).\n",
							oldNbor.getID(), oldNbor.getX(), oldNbor.getY(),
							other.getID(), other.getX(), other.getY()	);
				}
				neighborList.set( i, other );
			}
		}

		if ( model.rDebug > 2 ) {
			System.out.printf( "    neighbors: " );
			printNeighbors();
			System.out.printf( "\n" );
		}

		if ( model.getRDebug() > 1 ) 
			System.out.printf( "<- setupNeighborList for agent %d done.\n", ID );
	}

	///////////////////////////////////////////////////////////////////////////////////
	// pickRandomOtherAgent
	// pick a random agent from hhAgentList, but be sure its not the calling agent!
	//
	public Agent pickRandomOtherAgent () {
		Agent other;
		do {
		   	other = hhAgentList.get( Model.getUniformIntFromTo( 0, hhAgentList.size()-1  ) );
		} while ( other == this );
		return other;
	}

	public void printNeighbors () {
	   	for ( Agent o : neighborList ) 
	   		System.out.printf( " %d@(%d,%d)", o.getID(), o.getX(), o.getY() );
	}

	//////////////////////////////////////////////////////////////////////////////////
	// isNeighbor other
	// return true if other is on my neighborList, else false.
	//
	public boolean isNeighbor ( Agent other ) {
		if ( neighborList.contains( other ) )
			return true;
		return false;
	}


  	/**
  	 * What is written to the screen if the programmer asks that the agent be printed 
  	 * to the console. Here it is just the agents x,y coordinates.
  	 */
   	public String toString() {
		return x + " " + y;
  	}

	////////////////////////////////////////////////////////////////////////////////////
	/**
	 * The step procedure is called from the main menu for each agent.
	 * - first set scaledSocialNormExtractionLevel
	 * - then figure out an extraction level to maximize utility
	 *   if optimizer
	 *      calculate the extraction level
	 *   else not optimizer
	 *      pick best extraction level from a sample
	 * - then extract
	 *      
	 */
  	public void step() {
  		
  		timeClock++;
  		
		if ( model.getRDebug() > 2 ) 
			System.out.printf( "-> Agent %d step (searchType=%d).\n", ID, searchType );

  		// set the extraction level recommended by social norms
  		switch ( searchType ) {
			case Model.searchTypeRandom:
  				scaledSocialNormExtractionLevel = getRandomScaledSocialNormExtractionLevel();
				break;

			case Model.searchTypeVonNeumann:
				scaledSocialNormExtractionLevel = getNeighborsScaledSocialNormExtractionLevel();
				break;

			default:
				System.err.printf("\n\nERROR (Agent-step): illegal searchtype %d.\n\n", searchType );
				System.exit ( -1 );
		}
  		

		//  Now choose extraction level for present timestep (optimally or with limited cognition, 
		// depending on "optimizer=1" or otherwise). 
		// If non-optimizer, choose best from "cognitiveAbility" number of random picks  

  		if (optimizer==1){
			// Optimizing Agent, this agent has full information and knows the equation 
			// to determine the optimum level to extract from the forest
  			// DOES NOT YET WORK WITH THE NEW GATHERING TIME ALGORITHM

 			Ans = alphaNormSust;
// 			G = model.WALKING_SPEED;
 			L =  leisureEndowment;
 			Al = alphaLeisure;
 			Aser = socEnvRatio;
 			N = scaledSocialNormExtractionLevel;
 			S = scaledSustainabilityLevel;              
 			Ac = alphaConsumption;
  	  	  		
			// Solutions to Quadratic Eqn of form (B +- D)/2A
			// the solution to this equation should determine the optimum level of forest extraction (fuelwood)
 			double B = Ans*G*L + Cmax*Al*(1 + Aser*N + S - Aser*S) + Ac*(G*L + Cmax*(1 + Aser* N + S - Aser*S));
 			double A = (Ac + Al + Ans) * Cmax;
 			double C =  Ac * G*L * (1 + Aser*(N - S) + S);
  			double D = Math.sqrt(B*B - 4*A*C );
  			
			// the equation yields two solutions and we are only interested in largest of the two
  			double x1 = Math.min((B + D)/(2*A), 1.0);
  			double x2 = Math.min((B - D)/(2*A), 1.0);
  			
  			double bestUtility;
  			if ( getUtility(x1) > getUtility(x2)){                                             // if the utility from solution1 is greater than utility from solution2
  				scaledCurrentExtractionLevel = x1;                                                   // then set the current extraction level to solution1
  				bestUtility = getUtility(x1);                                                  // then set the hh utility value to that obtained using solution 1
  			}
  			
  			else {
  				scaledCurrentExtractionLevel = x2;                                                   // otherwise solution2 is optimal and set the current extraction level to it
  				bestUtility = getUtility(x2);                                                  // set the best utility to the utility value obtained using solution 2
  			}
  			
  			if (getUtility(1.0) > bestUtility){                                                // if the utility from full extraction is > then the best utility then
				scaledCurrentExtractionLevel = 1.0;                                                  // set the current extraction level to the maximum value 1.0
  			}
  			
			// System.out.print("  x1 = " + x1 + "  x2 = " + x2 + "  Level = " +currentExtractionLevel +"\n");
  		}
  

		// Limited cognition,
		// in this case the household agent evaluates a number 
		// of different random extraction level choices
		// and selects the extraction level that maximizes its utility function.
  		else {
  			double bestUtil = -1;
		  	double randExtractionLevel;
			scaledCurrentExtractionLevel = -1;    // LEAVE THIS OUT FOR AGENTS WHO GET BETTER OVER TIME???
			

			
			if (normalSearch==0 || (normalSearch==1 && useBurnInPeriod == 1 && timeClock <=20))
	  			for (int i=0; i < cognitiveAbility; i++ ) {
	  				randExtractionLevel = Random.uniform.nextDoubleFromTo(0, 1);
	  				if ( getUtility( randExtractionLevel ) > bestUtil ) {
	  					scaledCurrentExtractionLevel = randExtractionLevel;
	  					bestUtil = getUtility( randExtractionLevel );  
	  				}
	
	  			}
			else
	  			for (int i=0; i < cognitiveAbility; i++ ) {
	  				randExtractionLevel = 99999;
	  				while (randExtractionLevel>1.0 || randExtractionLevel < 0.0) {
		  				randExtractionLevel = Random.normal.nextDouble(scaledPreviousExtractionLevel, normalSearchStd);
	  				}
	  				if ( getUtility( randExtractionLevel ) > bestUtil ) {
	  					scaledCurrentExtractionLevel = randExtractionLevel;
	  					bestUtil = getUtility( randExtractionLevel );  
	  				}
	
	  			}

  		}	  		
  		
		// Now extract resource
		// If no extractable resource left, extract nothing
        // Otherwise limit the household to extracting between 0 & Cmax
		// (gr asked: should we substitute some other value for extractionLevelMax, like a standard deviation?
		// NOTE: convert level (which is in 0,1) into absolute amount

		absAmountTryToExtract =  Cmax * scaledCurrentExtractionLevel;

  		if ( resource.getAbsoluteQuantity() <= model.getResMinimumValue() ) {
			scaledCurrentExtractionLevel = 0.0;
			absAmountExtracted = 0.0;
  			}
		else {  // ask resource to give us what we ask for, but it might not have that much,
			// so revise currentExtractionLevel if needed.
			absAmountExtracted = resource.decrementResource( absAmountTryToExtract );
			if ( absAmountExtracted != absAmountTryToExtract )
				scaledCurrentExtractionLevel = ( absAmountExtracted / Cmax );
		}
		
		scaledTotalAmountExtracted += scaledCurrentExtractionLevel;  // increment the running total

		if ( scaledCurrentExtractionLevel <= model.lowExtrThres )  
			timeStepsAsLowExtractor = timeStepsAsLowExtractor + 1;
		else if ( scaledCurrentExtractionLevel >= model.highExtrThres ) 
			timeStepsAsHighExtractor = timeStepsAsHighExtractor + 1;
		else 
			timeStepsAsMedExtractor = timeStepsAsMedExtractor + 1;
		
		if ( model.getRDebug() > 2 )
			System.out.printf( "<- Agent-step done.\n" );

  	}

	//////////////////////////////////////////////////////////////////////////////////////
	// getUtility
	//
	// ARGUMENT sampleLevel: extraction level value between range 0-1, where 1 is full extraction
	//
	// PURPOSE: This method calculates the individual components of the utility
	// for transparency and then combines each of those components along with the 
	// extractionlevel value to determine the households utility given that extractionLevel.
	//
	// OUTPUT: the households utility at the given extraction level is returned
	// within the range of 0-1.

  	public double getUtility( double sampleLevel ) {
  		double  consumptionComponent;
  		if (sampleLevel>=scaledSubsistenceRequirement){
	  	   consumptionComponent = Math.pow (scaledSubsistenceRequirement + decreasingReturnsFactor*(sampleLevel-scaledSubsistenceRequirement), alphaConsumption);   
	       // if decreasingReturnsFactor = 1.0, then constant even beyond subsistenceRequirement
	  	   //  decreasingReturnsFactor < 1.0 then decreasing returns beyond subsistenceRequirement
  		}
  		else {
  			 consumptionComponent = Math.pow (sampleLevel, alphaConsumption);
  		}

  		double penaltyCost = 0.0;
  		if (Agent.model.getHasPenalty() == 1) {
			penaltyCost = model.getPenaltyCost(sampleLevel);
		} else {
			penaltyCost = 0;
		}
  		
       double sampleGatheringTime;
       sampleGatheringTime = Math.min(leisureEndowment, (sampleLevel*Cmax*600/30)*resource.getHeadloadGatheringTime());
       
       // Fit the penalty cost to the tail of the leisure function
       double leisureUsed = Math.min(1, (sampleGatheringTime+penaltyCost)/leisureEndowment);
       double leisureComponent = Math.pow( 1 - leisureUsed, alphaLeisure);
       // leisureComponent explanation -- sampleLevel*Cmax*600/30 = number of trips per month

		double socialNormComponent = (socEnvRatio/2)*(1-Math.abs(sampleLevel - scaledSocialNormExtractionLevel));
  		
  		double sustComponent = ((1-socEnvRatio)/2)*(1-Math.max(0, sampleLevel - scaledSustainabilityLevel));

		// System.out.println("\n socEnvRatio = " + socEnvRatio + "\n alphaConsumption = " + alphaConsumption + 
		// "\n alphaLeisure = " + alphaLeisure + "\n alphaNormSust = " + alphaNormSust);
		// System.out.println("\n Consumption Component = " + consumptionComponent + "\n Leisure Component = " + leisureComponent + 
		// "\n socialNormComponent = " + socialNormComponent + "\n sustComponent = " + sustComponent);  		
		// System.out.println("Total Utility = " + consumptionComponent * leisureComponent * Math.pow((socialNormComponent + sustComponent),alphaNormSust));
  		
		double utility = consumptionComponent * leisureComponent * 
				  Math.pow((socialNormComponent + sustComponent),alphaNormSust); 

  		return utility;
  	}
  	
  	public double getUtilityNew( double sampleLevel ) {
  		double nonTradedComponent;

		if (sampleLevel > scaledSubsistenceRequirement) {
			nonTradedComponent = scaledSubsistenceRequirement
					+ decreasingReturnsFactor
					* (sampleLevel - scaledSubsistenceRequirement);
		} else {
			nonTradedComponent = sampleLevel;
		}

		nonTradedComponent = Math.pow(nonTradedComponent, alphaConsumption);

		double tradedComponent;

		double gatherCost = (sampleLevel * Cmax * 600 / 30)
				* resource.getHeadloadGatheringTime() / leisureEndowment;

		double penaltyCost = 0.0;
		if (Agent.model.getHasPenalty() == 1) {
			penaltyCost = model.getPenaltyCost(sampleLevel);
		} else {
			penaltyCost = 0;
		}
		tradedComponent = Math.pow(Math.max(1 - (gatherCost + penaltyCost)
				/ leisureEndowment, 0), alphaLeisure);
		// System.out.printf("%f %f %f\n", scaledSustainabilityLevel,
		// penaltyCost, tradedComponent);

		double informalAdherenceComponent, formalAdherenceComponent;

		informalAdherenceComponent = (socEnvRatio / 2)
				* (1 - Math.abs(sampleLevel - scaledSocialNormExtractionLevel));
		formalAdherenceComponent = ((1 - socEnvRatio) / 2)
				* (1 - (sampleLevel - scaledSustainabilityLevel));
		double adherenceComponent = Math.pow(informalAdherenceComponent
				+ formalAdherenceComponent, alphaNormSust);

		// System.out.printf("%f %f %f\n", nonTradedComponent, tradedComponent,
		// adherenceComponent);

		return nonTradedComponent * tradedComponent * adherenceComponent;
  	}


	///////////////////////////////////////////////////////////////////////////////////////
	// postStep
	// set previousExtractionLevel to be value from current.
	//
  	public void postStep(){
  		scaledPreviousExtractionLevel = scaledCurrentExtractionLevel;
  	}
  	
	//////////////////////////////////////////////////////////////////////////////////////
  	// getRandomScaledSocialNormExtractionLevel
	// look at a random set of other agents, and from them 
	// calculate a scaledSocialNormExtractionLevel based on their scaledPreviousExtractionLevel's.
	//
  	public double getRandomScaledSocialNormExtractionLevel () {
  		Vector<Agent> neighbours = new Vector<Agent>( numRandomOthers );
  		for (int i = 0; i < numRandomOthers; i ++) {
  			neighbours.add( pickRandomOtherAgent() );
		}
  		return getScaledSocialNormExtractionLevel( neighbours );
	}
  	
	//////////////////////////////////////////////////////////////////////////////////////
	// getNeighborsScaledSocialNormExtractionLevel
	// look at all the agents neighbors (however they are set up),
	//  calculate a scaledSocialNormExtractionLevel based on their scaledPreviousExtractionLevel's
	//
  	public double getNeighborsScaledSocialNormExtractionLevel () {
		return getScaledSocialNormExtractionLevel( neighborList );
  	}
  	  	
	/////////////////////////////////////////////////////////////////////////////////////
	// get*scaledSocialNormExtractionLevel
	// get scaledPreviousExtractionLevel from neighbors, based on normType
	//    0    Mean of neighbors
	//    1    Max of neighbors
	//    2    Min of neighbors
	// for each type, just calculate it over list of provided neighbours
	//
  	public double getScaledSocialNormExtractionLevel( Vector<Agent> neighbours ) {
		if ( model.normType == 0 ) 
			return getMeanScaledSocialNormExtractionLevel( neighbours );
		else if ( model.normType == 1 ) 
			return getMaxScaledSocialNormExtractionLevel( neighbours ); 
		else 
			return getMinScaledSocialNormExtractionLevel( neighbours ); 
  	}
  	
  	public double getMeanScaledSocialNormExtractionLevel( Vector<Agent> neighbours ){
  		double total = 0;
  		for ( Agent a : neighbours )
  			total += a.getScaledPreviousExtractionLevel();
  		if ( neighbours.size() > 1 )
			total /= neighbours.size();
  		return total;
  	}

  	public double getMaxScaledSocialNormExtractionLevel( Vector<Agent> neighbours ){
  		double maxNorm = 0.0;
  		for ( Agent a : neighbours ) {
  			if ( a.getScaledPreviousExtractionLevel() > maxNorm )
  				maxNorm = a.getScaledPreviousExtractionLevel();
  		}
  		return maxNorm;
  	}
  	
  	public double getMinScaledSocialNormExtractionLevel( Vector<Agent> neighbours ){
  		double minNorm = 1.0;
  		for ( Agent a : neighbours ) {
  			if ( a.getScaledPreviousExtractionLevel() < minNorm )
  				minNorm = a.getScaledPreviousExtractionLevel();
  		}
  		return minNorm;
  	}

  	////////////////////////////////////////////////////////////////////////////////////////
	// Getters and setters

	// for class variables

	public static void setModel( Model m ) { model = m; }
	public static void setResource( Resource r ) { resource = r; }
	public static void setInstitution( Institution i ) { institution = i; }
	public static void setWorld( Object2DGrid w ) { world = w; }
	public static void setHHAgentList ( Vector<Agent> list ) { hhAgentList = list; }

	public static void setExtractionLevelMax ( double d ) { extractionLevelMax = d; }

	public int getID() { return ID; }

	// for instance variables
  	
  	public void setScaledSustainabilityLevel(double val){
  		scaledSustainabilityLevel = val;
  		absSustainabilityLevel = val * Cmax;
  	}

  	public void setAbsSustainabilityLevel(double val){
  		absSustainabilityLevel = val;
  		scaledSustainabilityLevel = val / Cmax;
  	}
  	
  	public void setSearchType(int val){
  		searchType = val;
  	}
  	
  	public double getScaledSocialNormExtractionLevel(){
  		return scaledSocialNormExtractionLevel;
  	}
  	
  	public double getScaledPreviousExtractionLevel(){
  		return scaledPreviousExtractionLevel;
  	}
  	
  	public void setScaledPreviousExtractionLevel(double val){
  		 scaledPreviousExtractionLevel = val;
  	}
  	
  	public double getScaledCurrentExtractionLevel(){
  		return scaledCurrentExtractionLevel;
  	}
  	
  	public double getScaledCurrentSustainabilityLevel() {
  		return scaledSustainabilityLevel;
  	}
  	
  	/**
	 * @return the scaledSubsistenceRequirement
	 */
	public double getScaledSubsistenceRequirement() {
		return scaledSubsistenceRequirement;
	}

	/**
	 * @param scaledSubsistenceRequirement the scaledSubsistenceRequirement to set
	 */
	public void setScaledSubsistenceRequirement(double scaledSubsistenceRequirement) {
		this.scaledSubsistenceRequirement = scaledSubsistenceRequirement;
	}

	/**
	 * @return the scaledSustainabilityLevel
	 */
	public double getScaledSustainabilityLevel() {
		return scaledSustainabilityLevel;
	}

	/**
	 * @param scaledSocialNormExtractionLevel the scaledSocialNormExtractionLevel to set
	 */
	public void setScaledSocialNormExtractionLevel(
			double scaledSocialNormExtractionLevel) {
		this.scaledSocialNormExtractionLevel = scaledSocialNormExtractionLevel;
	}

	public void sethhAgentList(Vector<Agent> val) {
  		hhAgentList = val;
  	}

	public Vector<Agent> getNeighborList () {
		return neighborList;
	}

  	public double getHHSize(){
  		return hhSize;
  	}
  	
  	public int getX() {
    	return x;
  	}
      	
  	public int getY() {
  		return y;
  	}
  	
  	public double getAlphaConsumption(){
  		return alphaConsumption;
  	}

  	public double getAlphaLeisure(){
  		return alphaLeisure;
  	}

  	public double getSocEnvRatio(){
  		return socEnvRatio;
  	}

  	public double getAlphaNormSust(){
  		return alphaNormSust;
  	}
  	

  	public int getTimeStepsAsLowExtractor(){
  		return timeStepsAsLowExtractor;
  	}  

	public int getTimeStepsAsMedExtractor(){
  		return timeStepsAsMedExtractor;
  	}  	

	public int getTimeStepsAsHighExtractor(){
  		return timeStepsAsHighExtractor;
  	}
  	
	public int getClusterID(){
  		return clusterID;
  	}
	
	public double getScaledTotalAmountExtracted () {
		return scaledTotalAmountExtracted;
	}

	public double getAbsTotalAmountExtracted () {
		return scaledTotalAmountExtracted*Cmax;
	}
	
  	public void draw(SimGraphics g) {
  		//int shade = (int)(belief*255);
  		//agentShade = new Color(shade,0,0);
  		//g.drawFastCircle ( agentShade);
  		if (scaledCurrentExtractionLevel >= model.getHighExtrThres() ) {
  			agentShade = new Color(255,0,0);
  		}
  		else if (scaledCurrentExtractionLevel <= model.getLowExtrThres()  ) {
  			agentShade = new Color(0,255,0);
  		}
  		else {
  			agentShade = new Color(0,0,255);
  		}

  		g.drawFastCircle(agentShade);
  	}
  	

  	
}

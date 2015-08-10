package Belief2;


import java.util.Vector;

import uchicago.src.sim.engine.Schedule;
import uchicago.src.sim.space.Object2DGrid;
import uchicago.src.sim.util.Random;
import cern.jet.random.Normal;
import cern.jet.random.Uniform;
import cern.jet.random.engine.MersenneTwister;

import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

public class Model extends ModelParameters {
	
	// Some public class variables to be symbolic constants
	public static final int searchTypeRandom = 0;
	public static final int searchTypeVonNeumann = 1;

	public Schedule schedule;

	public Object2DGrid world;  		   	// agent space (grid)
	public int sizeX;					   	// X-size of world
	public int sizeY;					   	// Y-Size of world
	private Vector<Agent> agentList;		// list of the agents in the grid

	private int population;                 // total number of *people* ( approx = hhSizeMean * (sizeX*sizeY) )

	private int searchType;			 		// type of neighbourhood  
											// 0 = random each step; 1 = fixed vonNeumann (perhaps with rewiring)
	private int numRandomOthers;			// number of random others for each agent for searchType=0
    private double probPickRandomNeighbor;  // prob. replace initial (eg von Neumann) neighbor with random other agent (searchType=1)

	private int clusterType;                      //  =1 (2 bands, upper band has 'clusterSize*SizeX*' number of rows - when clusterSize=1, effectively no clusters)
	private double clusterSize;		              //  = 2 (box within a box, inner box has rows=columns=clusterSize*sizeX)
	private int clusterID;
	

    int optimizer;                              // This variable sets the household decision making to optimize (1) or tells the household to optimize over a bounded set of extraction levels (0)
	int cognitiveAbility;                       // when optimizer = 0 this variable defines the number of different amounts of extraction level that make up the bounded set from which the household optimizes its utility
	public int normType;						// When an agent looks at the extraction level of the social norms then it may retrive it using 0 = mean; 1 = max; 2 = min
	private int normalSearch;				// =0 implies boundedly rational agent evaluates sample extraction levels drawn from a uniform distribution. =1 implies normal distribution centered on previous extraction level
	private double normalSearchStd;			// Standard Dev of normal distribution used to draw sample extraction levels, if normalSearch = 1
	private int useBurnInPeriod;

	/**
	 * Penalty Parameters
	 */
	int hasPenalty;		//0: no penalty levied
						//1: penalty levied
	double penaltyH;	// "threshold" - higher values make the descent sharper for smaller inputs 
	double penaltyS;	// "shape" - higher values make the curvature higher

	
//    public double WALKING_SPEED = 4.3;      // in km /hr

	/**
	 * @return the hasPenalty
	 */
	public int getHasPenalty() {
		return hasPenalty;
	}


	/**
	 * @param hasPenalty the hasPenalty to set
	 */
	public void setHasPenalty(int hasPenalty) {
		this.hasPenalty = hasPenalty;
	}


	/**
	 * @return the penaltyH
	 */
	public double getPenaltyH() {
		return penaltyH;
	}


	/**
	 * @param penaltyH the penaltyH to set
	 */
	public void setPenaltyH(double penaltyH) {
		this.penaltyH = penaltyH;
	}


	/**
	 * @return the penaltyS
	 */
	public double getPenaltyS() {
		return penaltyS;
	}


	/**
	 * @param penaltyS the penaltyS to set
	 */
	public void setPenaltyS(double penaltyS) {
		this.penaltyS = penaltyS;
	}
	
	public double getPenaltyCost(double level) {
		return (1 + Math.pow(penaltyH,penaltyS) * Math.pow(1 - level, penaltyS) / (Math.pow(1-level, penaltyS) + Math.pow(penaltyH, penaltyS)));
	}

	public double lowExtrThres;					// defines an extraction level threshold whereby agents with extraction levels lower then the threshhold are "Low" extractors
	public double highExtrThres;				// defines an extraction level threshold whereby agents with extraction levels higher then the threshold are "High" extractors	

    // HOUSEHOLD DISTRIBUTION VARIABLES
	private double hhSizeMean;                  // the mean size of agent households
	private double hhSizeStd;                   // the standard deviation of the agent household sizes


	private int useStringForAlphaMeans;
	private int useIdenticalStds;
	private int useIdenticalAserms;
	
	private double decreasingReturnsFactor;    // < 1.0 implies decreasing returns beyond subsistenceRequirement
	
	// Cluster 1
	double alphaConsumptionMean1;                // each household has a preference weight for consumption, leisure, weighing social norms versus instutional
	double alphaLeisureMean1;                    // sustainability measurements [alphaNormSustMean], and a weight applied to the overall norms and environment signals [SocEnvRatio]
	double alphaNormSustMean1;
	double alphaConsumptionStd1;
	double alphaLeisureStd1;
	double alphaNormSustStd1;
	double alphaSocEnvRatioMean1;
	double alphaSocEnvRatioStd1;	
	String alphaMeans1;

	// Cluster 2
	double alphaConsumptionMean2;                // each household has a preference weight for consumption, leisure, weighing social norms versus instutional
	double alphaLeisureMean2;                    // sustainability measurements [alphaNormSustMean], and a weight applied to the overall norms and environment signals [SocEnvRatio]
	double alphaNormSustMean2;
	double alphaConsumptionStd2;
	double alphaLeisureStd2;
	double alphaNormSustStd2;
	double alphaSocEnvRatioMean2;
	double alphaSocEnvRatioStd2;
	String alphaMeans2;
	
	double MASTER_STD; 						// overrides all the std's (if not commented out in buildmodel)


	public double iniScaledExtrLevMean;               // the mean of the distribution of extraction levels that are assigned to households at the initialization of the model
	public double iniScaledExtrLevStd;                // the standard deviation of the extraction levels that are assigned to households at the initialization of the model
	public int iniExtrType;					    // initializing extraction levels for the population may take the following forms
                                                // 0 => random uniform distribution, 1 => normal distribution, 2 => specify the number of agents of each high,med,low
	                                            // such that within those groups agents are informed uniformly between the threshold ranges of those categories
	public int noHigh;						    // No of initial high extraction level category agents, when iniExtrType=2
	public int noMed;					   	    // No of initial medium extraction level category agents, when iniExtrType=2


	// RESOURCE VARIABLES
	private Resource forestResource;            // this resource object holds all of the forest properties
  	private double resGrowthMean;               // variable to hold the mean value of the growth rate of the resource
  	private double resGrowthStd;                // variable to hold the std dev in the resource growth rate
  	private double initialResQuantity;                 // used to set the initial resource quantity
	private int resourceExtinct;                // Flag to 1 when resource goes down to zero
	private double resMinimumValue;				// reserve level of resource
	private double initialPerCapitaForestArea;           // Initial area in ha
	private double gatheringTimeFactor;			// gathering time for one headload = 2 + (resInitial/resRemaining)^gatheringTimeFactor;

		
	// INSTITUTION VARIABLES
	private Institution institution;            // we only have one institutional agent in this model and this is it
	public int instPerCapitaAllocation;		    // 0=> institution makes rules per capita, 1=> per household (equally across households regardless of size)
	public int instFrequency;				    // time separating institutional announcements (no of months)


	// ***  Measurement variables  ***

   	// summary variables that is used to hold the total number of high/med/low extracting households 
	// that will be output each timestep to an output file and/or chart
	public int numberHighExtractors;
	public int numberMedExtractors;
	public int numberLowExtractors;

  	private double totalExtractionInPeriod;     // Total resource extraction in the present time step
	private double totalExtractionSinceBeginning; // Total resource extraction since model was initalized
	private int stepFirstHitMinResource;		// time step resource level first hits resMinimumValue
	private int	 numberTimesAtMinimumRes;		// count number of steps resource <= minimum	

	
	private double meanAgentAbsTotalExtractionAmount; // mean across agents of abs total they extracted over run
	private double sdAgentAbsTotalExtractionAmount;

	private double meanAgentAbsTotalExtractionAmountByASERM;   // mean total extraction in cluster 1 (temp variable) 
	private double sdAgentAbsTotalExtractionAmountByASERM;
//	private double meanAgentScaledTotalExtractionAmountByASERM2;
//	private double sdAgentScaledTotalExtractionAmountByASERM2;
	
	private int hhMonthsLow;                    // Total household-months of Low extraction observed
	private int hhMonthsMed;                    // Total household-months of Medium extraction observed
	private int hhMonthsHigh;                   // Total household-months of High extraction observed





	/////////////////////////////////////////////////////////////////////////////
	// addModelSpecificParameters
	// add alias and long name for Model parameters you want to set at run time
	// the long name should be same as instance variable
	//
	// Note: the generic parameters from ModelParameters are already available.
	
	public void addModelSpecificParameters () {
		parametersMap.put( "X", "sizeX" );
		parametersMap.put( "Y", "sizeY" );
		parametersMap.put( "resGrM", "resGrowthMean" );
		parametersMap.put( "resGrS", "resGrowthStd" );
//		parametersMap.put( "resQ", "resQuantity" );	
		parametersMap.put( "resMinV", "resMinimumValue" );
		parametersMap.put( "iniPerCapArea", "initialPerCapitaForestArea" );
		parametersMap.put( "gTF", "gatheringTimeFactor" );

		parametersMap.put( "sT", "searchType" );	
		parametersMap.put( "nmRO", "numRandomOthers" );	
		parametersMap.put( "prPRN", "probPickRandomNeighbor" );
		
		parametersMap.put( "hhSzM", "HHSizeMean" );
		parametersMap.put( "hhSzV", "HHSizeStd" );

		parametersMap.put( "loExTh", "lowExtrThres" );	
		parametersMap.put( "hiExTh", "highExtrThres" );
		
		parametersMap.put( "useString", "useStringForAlphaMeans" );
		parametersMap.put( "useIdenticalStds", "useIdenticalStds" );
		parametersMap.put( "decRetFact", "decreasingReturnsFactor" );
		parametersMap.put( "aConM1", "alphaConsumptionMean1" );
		parametersMap.put( "aLM1", "alphaLeisureMean1" );
		parametersMap.put( "aNSM1", "alphaNormSustMean1" );		
		parametersMap.put( "aSERM1", "alphaSocEnvRatioMean1" );
		parametersMap.put( "aConS1", "alphaConsumptionStd1" );
		parametersMap.put( "aLS1", "alphaLeisureStd1" );
		parametersMap.put( "aNSS1", "alphaNormSustStd1" );
		parametersMap.put( "aSERS1", "alphaSocEnvRatioStd1" );
		parametersMap.put( "aM1", "alphaMeans1" );
		parametersMap.put( "useIdenticalAserms", "useIdenticalAserms" );
		parametersMap.put( "aConM2", "alphaConsumptionMean2" );
		parametersMap.put( "aLM2", "alphaLeisureMean2" );
		parametersMap.put( "aNSM2", "alphaNormSustMean2" );		
		parametersMap.put( "aSERM2", "alphaSocEnvRatioMean2" );
		parametersMap.put( "aConS2", "alphaConsumptionStd2" );
		parametersMap.put( "aLS2", "alphaLeisureStd2" );
		parametersMap.put( "aNSS2", "alphaNormSustStd2" );
		parametersMap.put( "aSERS2", "alphaSocEnvRatioStd2" );
		parametersMap.put( "aM2", "alphaMeans2" );
		
		parametersMap.put( "MASTER_STD", "MASTER_STD" );

		parametersMap.put( "iELevM", "iniScaledExtrLevMean" );
		parametersMap.put( "iELevS", "iniScaledExtrLevStd" );
		parametersMap.put( "iEType", "iniExtrType" );
		parametersMap.put( "noHi", "noHigh" );
		parametersMap.put( "noMed", "noMed" );
		parametersMap.put( "opt", "optimizer" );
		parametersMap.put( "cogA", "cognitiveAbility" );
		parametersMap.put( "normalSearch", "normalSearch" );
		parametersMap.put( "normalSearchStd", "normalSearchStd" );
		parametersMap.put( "insPCA", "instPerCapitaAllocation" );
		parametersMap.put( "insFreq", "instFrequency" );
		parametersMap.put( "normType", "normType" );
		parametersMap.put( "clusterType", "clusterType" );
		parametersMap.put( "cSiz", "clusterSize" );
		parametersMap.put("hasPenalty", "hasPenalty");
		parametersMap.put("penaltyS", "penaltyS");
		parametersMap.put("penaltyH", "penaltyH");
		
		parametersMap.put( "burnIn", "useBurnInPeriod" );
	}
	
	
	// control what appears in the repast parameter panel
	public String[] getInitParam () {
		String[] params = { "sizeX", "sizeY",  "resGrowthMean", "resGrowthStd", "resMinimumValue", "initialPerCapitaForestArea", "gatheringTimeFactor",
					"socEnvRatio", "searchType", "numRandomOthers", "probPickRandomNeighbor", 
                    "hhSizeMean", "hhSizeStd", "lowExtrThres", "highExtrThres", "useStringForAlphaMeans", "useIdenticalStds", "useIdenticalAserms", "decreasingReturnsFactor",
					"alphaConsumptionMean1", "alphaLeisureMean1", "alphaNormSustMean1", "alphaMeans1", "alphaConsumptionStd1",
					"alphaLeisureStd1", "alphaNormSustStd1", "alphaSocEnvRatioMean1", "alphaSocEnvRatioStd1",  
					"alphaConsumptionMean2", "alphaLeisureMean2", "alphaNormSustMean2", "alphaMeans2", "alphaConsumptionStd2",
					"alphaLeisureStd2", "alphaNormSustStd2", "alphaSocEnvRatioMean2", "alphaSocEnvRatioStd2",  "MASTER_STD",
					"iniScaledExtrLevMean", "iniScaledExtrLevStd", "iniExtrType", "noHigh", "noMed", "optimizer", "cognitiveAbility", "normalSearch", "normalSearchStd",
							"instPerCapitaAllocation", "instFrequency", "normType", "clusterType", "clusterSize", "useBurnInPeriod",
							"hasPenalty", "penaltyH", "penaltyS",
				// these are from the super class:
					"rDebug", "seed" };
		return params;
	}

	////////////////////////////////////////////////////////////////////////////
	// constructor, if need to do anything special.
	public Model () {
	}

	///////////////////////////////////////////////////////////////////////////
	// setup
	// set defaults after a run start or restart
	// see the variable declarations at the top of this .java class for variable definitions.

	public void setup () {
		super.setup();
		System.gc();                                                                
		agentList = new Vector<Agent>();	

		setParameterDefaultValues();

		if ( rDebug > 0 )
			System.out.printf( "==> setup...\n" );
		schedule = null;

		System.gc ();                                                      // garabage collection of discarded objects
		super.setup();                                                     // THIS SHOULD BE CALLED after setting defaults in setup().
		schedule = new Schedule (1);                                       // create AFTER calling super.setup()

		if ( rDebug > 0 )
			System.out.printf( "\n<=== Model-setup() done.\n" );
	}

	///////////////////////////////////////////////////////////////////////////
	// setParameterDefaultValues
	// 
	public void setParameterDefaultValues () {
		// Village 19 has the most households from the changar.sta survey.
		// Village 19 has 100 unique household id values, so lets make the grid size 10x10
		sizeX = 25;																	
		sizeY = 25;																	

		hhSizeMean = 4.75;
		hhSizeStd = 2.0;
										
		searchType = 1;
		numRandomOthers = 4;                // 4 of them
		probPickRandomNeighbor = 0.0;

		normType = 0;
		clusterType = 1;
		clusterSize = 0.25;

		resGrowthMean = 0.002;  // Equivalient to 2.5% per year from USFS Birdsey 1992
		resGrowthStd = 0.001;
		resMinimumValue = 250;
		initialPerCapitaForestArea = 0.75;           // in ha / person
		initialResQuantity = 0.5*4.1*(initialPerCapitaForestArea*sizeX*sizeY*hhSizeMean)*10000/600; 
						// = 0.5 proportion in branches * 4.1 kg / m2 * initialForestArea ha * 10^4 m2/ha /  600 kg/m3
		gatheringTimeFactor = 1.0;

		useStringForAlphaMeans=1;
		useIdenticalStds=1;
		decreasingReturnsFactor=0.2;
		
		// Cluster 1 
		alphaMeans1 = "0.33_0.33_0.33";  // {consumption, leisure, norm-sust};
		alphaSocEnvRatioMean1 = 1.0;

		
		useIdenticalAserms=0;

		// Cluster 2 		
		alphaMeans2 = "0.33_0.33_0.33";  // {consumption, leisure, norm-sust};
		alphaSocEnvRatioMean2=1.0;

		/**
		 * Penalty initialization
		 */
		hasPenalty = 0;
		penaltyH = 0.4;
		penaltyS = 4.0;
		
		MASTER_STD=0.0;
	
		optimizer=0;
		cognitiveAbility=10;
		normalSearch = 0;
		normalSearchStd=0.01;
		useBurnInPeriod = 1;

		lowExtrThres = 0.50;
		highExtrThres = 0.75;
		iniExtrType = 1;															
		iniScaledExtrLevMean = 0.5;
		iniScaledExtrLevStd = 0.2;
		noHigh=300;
		noMed=300;

		instPerCapitaAllocation = 1;
		instFrequency = 12;

		// measurement variables

		totalExtractionSinceBeginning = 0.0;
		stepFirstHitMinResource = 9999999;
		numberTimesAtMinimumRes = 0;
		resourceExtinct = 0;
		hhMonthsLow = 0;
		hhMonthsMed = 0;
		hhMonthsHigh = 0;

	}


	///////////////////////////////////////////////////////////////////////////
	// buildModel
	// We build the "conceptual" parts of the model.
	// (vs the display parts, and the schedule)
	//
	// Create a 2D world, tell the Agents about it.
	// Create agents, put them in the world, one per cell,
	// and also add them to the agentList.

	public void buildModel () {
		if ( rDebug > 0 )
			System.out.printf( "==> buildModel...\n" );

		// CALL FIRST -- defined in super class -- it starts RNG, etc
		buildModelStart();

		///////////////////////////////////////////////////////////////////////
		// *** Change below here -- initialization	

		double alphaLeisure;
		double alphaConsumption;
		double alphaNormSust;
		double socEnvRatio;
		double alphaSum;                                                                     // Used as a temp variable to normalize alpha's
		

// Initialize alpha values for Cluster 1 
		//Initialize Means of distributions
		if (useStringForAlphaMeans==1) {
			String alphaMeansArray1[] = alphaMeans1.split("_");
			alphaConsumptionMean1 = Double.parseDouble(alphaMeansArray1[0].trim());
			alphaLeisureMean1 = Double.parseDouble(alphaMeansArray1[1].trim());
			alphaNormSustMean1 = Double.parseDouble(alphaMeansArray1[2].trim());
		}

		
		// Initialize standard deviations of distributions
		if (useIdenticalStds==1) {
			alphaConsumptionStd1=alphaLeisureStd1=alphaNormSustStd1=alphaSocEnvRatioStd1=MASTER_STD;
		}

		
// Initialize alpha values for Cluster 2 
		// Initialize means of distributions
		if (useStringForAlphaMeans==1) {
			String alphaMeansArray2[] = alphaMeans2.split("_");
			alphaConsumptionMean2 = Double.parseDouble(alphaMeansArray2[0].trim());
			alphaLeisureMean2 = Double.parseDouble(alphaMeansArray2[1].trim());
			alphaNormSustMean2 = Double.parseDouble(alphaMeansArray2[2].trim());
		}


		// Initialize standard deviations of distributions
		if (useIdenticalStds==1) {
			alphaConsumptionStd2=alphaLeisureStd2=alphaNormSustStd2=alphaSocEnvRatioStd2=MASTER_STD;
		}


		if (useIdenticalAserms==1) {
			alphaSocEnvRatioMean2=alphaSocEnvRatioMean1;
		}
		
		// To keep alphaLeisureMean and alphaNormSustMean between [0,1]
		alphaLeisureMean1 = Math.min(alphaLeisureMean1, 1-alphaConsumptionMean1);
		alphaLeisureMean2 = Math.min(alphaLeisureMean2, 1-alphaConsumptionMean2);
		alphaNormSustMean1 = Math.max(0, 1.0-alphaConsumptionMean1-alphaLeisureMean1);
		alphaNormSustMean2 = Math.max(0, 1.0-alphaConsumptionMean2-alphaLeisureMean2);
				

		world = new Object2DGrid(sizeX,sizeY);                                              // instantiate the agent world using the x,y dimensions specified above
	   
		                                                                                    // create distributions for randomly drawing agent characteristics 
		Normal hhDist;                                                                      // household distribution used for household size
//		Normal resVarDist;                                                                  // resource variation distribution
//		Normal iniScaledExtrLevDistNorm;                                                          // normal distribution of possible extraction level preferences
//		Uniform iniScaledExtrLevDistUni;                                                          // uniform distribution of possible extraction level preferences

		Normal conPrefDist1;                                                                 // distribution of consumption preference weight values BAND 1
		Normal leisPrefDist1;                                                                // distribution of leisure preference weight values
		Normal normSustPrefDist1;                                                            // distribution of norm+sustainability preference weight
		Normal socEnvRatioPrefDist1;                                                         // distribution of possible preference weights for social norm/institutional factor

		Normal conPrefDist2;                                                                 // distribution of consumption preference weight values
		Normal leisPrefDist2;                                                                // distribution of leisure preference weight values BAND 2
		Normal normSustPrefDist2;                                                            // distribution of norm+sustainability preference weight
		Normal socEnvRatioPrefDist2;                                                         // distribution of possible preference weights for social norm/institutional factor

		
		MersenneTwister generator1 = new MersenneTwister((int)seed);                        // create a random number generator that takes a seed specified by the user
        hhDist = new Normal(hhSizeMean, hhSizeStd, generator1);     						// the normal function uses the standard deviation and not the variance


		// BAND 1
		conPrefDist1 = new Normal(alphaConsumptionMean1, alphaConsumptionStd1, generator1);    // that all follow from the initial seed and therefore allow for run replication
        leisPrefDist1 = new Normal(alphaLeisureMean1, alphaLeisureStd1, generator1);
		normSustPrefDist1 = new Normal(alphaNormSustMean1, alphaNormSustStd1, generator1);
		socEnvRatioPrefDist1 = new Normal(alphaSocEnvRatioMean1, alphaSocEnvRatioStd1, generator1);

		// BAND 2
        conPrefDist2 = new Normal(alphaConsumptionMean2, alphaConsumptionStd2, generator1);    // that all follow from the initial seed and therefore allow for run replication
        leisPrefDist2 = new Normal(alphaLeisureMean2, alphaLeisureStd2, generator1);
		normSustPrefDist2 = new Normal(alphaNormSustMean2, alphaNormSustStd2, generator1);
		socEnvRatioPrefDist2 = new Normal(alphaSocEnvRatioMean2, alphaSocEnvRatioStd2, generator1);


		// object to control distribution of growth rate changes and initial quantity of the resource
		initialResQuantity = 0.5*4.1*(initialPerCapitaForestArea*sizeX*sizeY*hhSizeMean)*10000/600; 
		// = 0.5 proportion in branches * 4.1 kg / m2 * initialForestArea ha * 10^4 m2/ha /  600 kg/m3
        forestResource = new Resource( this, initialResQuantity, resGrowthMean, resGrowthStd, gatheringTimeFactor);
		Resource.setResMinimumValue( resMinimumValue );   // the "reserve" level of the resource

		// create institution with the time lapse (months) between institutional annoucements
		// and the type of calculation (either per capita [1] or per households [0])
        institution = new Institution( this, instFrequency, instPerCapitaAllocation );


		// Tell the Agents about the other major objects
		Agent.setModel( this );
		Agent.setResource( forestResource );
		Agent.setWorld ( world );
		Agent.setInstitution( institution );
		Agent.setHHAgentList( agentList );
        
        int hhSize;  

		// populate the entire grid, starting with the nonEnv agents.

		for (int j = 0; j < world.getSizeX(); j++){
			for (int i = 0; i < world.getSizeY(); i++){
				hhSize = hhDist.nextInt();
				while (hhSize <= 0){
					hhSize = hhDist.nextInt();
				}

				double initialScaledExtractionLevel;  // get this depending on the iniExtrType
				if (iniExtrType == 0 ){
					// uniform distribution draw
					// *** KLUGE ALERT -- note we use the mean and some constants to get a range-- argh!
					double minExtr = Math.max( 0.0, iniScaledExtrLevMean / 2.0 );
					double maxExtr = Math.min( 1.0, iniScaledExtrLevMean * 1.5 );
					initialScaledExtractionLevel = getUniformDoubleFromTo( minExtr, maxExtr );
				}
				else if (iniExtrType==1) {
					// normal distribution draw (make sure that the random draw is in [0,1])
					initialScaledExtractionLevel = getNormalDoubleProb ( iniScaledExtrLevMean, iniScaledExtrLevStd );
				}
				
				// ***>> NOTE:  Won't this next bit 'segregate' the levels spatially??? <<****
				else {   // Specify the exact numbers of each type
					if ( j*world.getSizeY()+i < noHigh ) 
						initialScaledExtractionLevel = getUniformDoubleFromTo( highExtrThres, 1.0 );
					else if ( j*world.getSizeY()+i < noHigh+noMed ) 
						initialScaledExtractionLevel = getUniformDoubleFromTo( lowExtrThres, highExtrThres );
					else 
						initialScaledExtractionLevel = getUniformDoubleFromTo( 0.0, lowExtrThres ); 
						
				}
				
				// Initialize our weights to an unacceptable value
				alphaLeisure = -1.0;
				alphaConsumption = -1.0;
				alphaNormSust = -1.0;
				socEnvRatio = -1.0;

	
				// TWO STRIPES (INCLUDING ONE GIANT "STRIPE")
				if (clusterType==1) {     
					if ( j<=clusterSize*sizeX ) {
						clusterID=1;						
						while (alphaLeisure < 0 ||  alphaLeisure > 1){
							alphaLeisure = leisPrefDist1.nextDouble();
						}
						while (alphaConsumption < 0 || alphaConsumption > 1 ){
							alphaConsumption = conPrefDist1.nextDouble();
						}
						while (alphaNormSust < 0 || alphaNormSust >1) {
							alphaNormSust = normSustPrefDist1.nextDouble();
						}
						while (socEnvRatio < 0 || socEnvRatio > 1){
							socEnvRatio = socEnvRatioPrefDist1.nextDouble();
						}				
					}

					else {
						clusterID=2;
						while (alphaLeisure < 0 ||  alphaLeisure > 1){
							alphaLeisure = leisPrefDist2.nextDouble();
						}
						while (alphaConsumption < 0 || alphaConsumption > 1 ){
							alphaConsumption = conPrefDist2.nextDouble();
						}
						while (alphaNormSust < 0 || alphaNormSust >1) {
							alphaNormSust = normSustPrefDist2.nextDouble();
						}
						while (socEnvRatio < 0 || socEnvRatio > 1){
							socEnvRatio = socEnvRatioPrefDist2.nextDouble();
						}									
					}
				}
				
				
//				SETS UP BOX-WITHIN-BOX CLUSTERS 				
				else if (clusterType==2){
//					 This inelegant calculation sets up the cluster: TRUE -> outer box; FALSE -> inner box
					
					if ( (j<(sizeX-clusterSize*sizeX)/2 || j>(sizeX+sizeX*clusterSize)/2) || (i<(sizeY-sizeY*clusterSize)/2 || i>(sizeY+sizeY*clusterSize)/2)) 	{
						clusterID=1;
						while (alphaLeisure < 0 ||  alphaLeisure > 1){
							alphaLeisure = leisPrefDist1.nextDouble();
						}
						while (alphaConsumption < 0 || alphaConsumption > 1 ){
							alphaConsumption = conPrefDist1.nextDouble();
						}
						while (alphaNormSust < 0 || alphaNormSust >1) {
							alphaNormSust = normSustPrefDist1.nextDouble();
						}
						while (socEnvRatio < 0 || socEnvRatio > 1){
							socEnvRatio = socEnvRatioPrefDist1.nextDouble();
						}
					}
					
					else {    // inner box
						clusterID=2;
						while (alphaLeisure < 0 ||  alphaLeisure > 1){
							alphaLeisure = leisPrefDist2.nextDouble();
						}
						while (alphaConsumption < 0 || alphaConsumption > 1 ){
							alphaConsumption = conPrefDist2.nextDouble();
						}
						while (alphaNormSust < 0 || alphaNormSust >1) {
							alphaNormSust = normSustPrefDist2.nextDouble();
						}
						while (socEnvRatio < 0 || socEnvRatio > 1){
							socEnvRatio = socEnvRatioPrefDist2.nextDouble();
						}									
					}
				}



				// Normalization -- making sure alpha's sum to 1
				alphaSum=alphaLeisure+alphaConsumption+alphaNormSust;
				alphaLeisure = alphaLeisure/alphaSum;
				alphaConsumption = alphaConsumption/alphaSum;
				alphaNormSust = alphaNormSust/alphaSum;
				
				// create a new agent
				Agent agent = new Agent (i, j, optimizer, cognitiveAbility, normalSearch, normalSearchStd, socEnvRatio, alphaNormSust, 
						alphaConsumption, alphaLeisure, searchType, hhSize, initialScaledExtractionLevel, clusterID, decreasingReturnsFactor, useBurnInPeriod );

				// put the agent into the agent world (grid), add the agent to the list of households, 
				// and add the popualtion of the household agent to the population count
				world.putObjectAt(i,j,agent);
				agentList.add(agent);
				population = population + hhSize;

				recordAgentExtractionLevel ( agent );	// incr number*Extractors and hhMonths* vars

			}
		}		

		hhMonthsLow = hhMonthsMed = hhMonthsHigh = 0;  // set these counters to zero

		// Now we are ready for agents to set up their neighborList (if searchType != random)
		if ( searchType != searchTypeRandom ) {
			if ( rDebug > 0 )  System.out.printf( "Setup non-random neighborList for each agent:\n" );
			for ( Agent a : agentList )
				a.setUpNeighborList ();
			if ( rDebug > 0 )
				checkNeighborLists();
		}
		
		// *** End Change -- initialization
		///////////////////////////////////////////////////////////////////////

		// some post-load finishing touches
		startReportFile();

		// you probably don't want to remove any of the following
		// calls to process parameter changes and write the
		// initial state to the report file.
		// NB -> you might remove/add more agentChange processing
        applyAnyStoredChanges();
        stepReport();
        getReportFile().flush();
        getPlaintextReportFile().flush();

		if ( rDebug > 0 )
			System.out.printf( "<== buildModel done.\n" );
	}

	///////////////////////////////////////////////////////////////////////////
	// checkNeighborLists
	// This does a quick check to see how clustered the neighbor are.
	// For each neighbor of agent A, count how many of those are neighbors of
	// each other.  (Divide by 2 since they are symmetric.)
	//
	public void checkNeighborLists () {
		double totalNborsNbors = 0;

		System.out.printf( "\n\n-> checkNeighborsLists:\n" );

		for ( Agent a : agentList ) {
			int agentsNborsNbors = 0;
			Vector<Agent> neighborList = a.getNeighborList();
			// note nbor doesn't have to run to last entry, since there are none after that to x-check
			for ( int nbor = 0; nbor < neighborList.size()-1; ++nbor ) {
				Agent nborAgent = neighborList.get( nbor );
				if ( rDebug > 1 ) {
					System.out.printf( "   Checking agent %d@%d,%d, neighbors: ", 
								   nborAgent.getID(), nborAgent.getX(), nborAgent.getY() );
					nborAgent.printNeighbors();
					System.out.printf( "\n" );
				}
				for ( int nbor1 = nbor + 1; nbor1 < neighborList.size(); ++nbor1 ) {
					if ( nborAgent.isNeighbor( neighborList.get( nbor1 ) ) )
						++agentsNborsNbors;
				}
			}
			if ( rDebug > 1 ) 
				System.out.printf( "   Agent %d agentsNborsNbors=%d.\n", a.getID(), agentsNborsNbors );
			totalNborsNbors += agentsNborsNbors;
		}

		// divide by 2 because we counted neighbors of neighbors both ways (A <-> B)
		double avgNborsNbors = totalNborsNbors / ( agentList.size() * 2 );

		System.out.printf( "<- checkNeighbrsLists done.  avgNborsNbors=%.2f.\n\n", avgNborsNbors );

	}


	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	// step
	// The top of the "conceptual" model's main dynamics
	// - zero out some counters
	// - send step() to forestResource 
	// - send step() to institution
	// - for each agent
	//      step, then do some counting
	// - for each agent
	//      postStep()
	// - update some totals
	// - call stepReport()

	public void step () {

		totalExtractionInPeriod = 0.0;

		forestResource.step();
		institution.step();

		uchicago.src.sim.util.SimUtilities.shuffle(agentList);
		numberHighExtractors = 0;
		numberMedExtractors = 0;
		numberLowExtractors = 0;
		
		for (int i = 0; i < agentList.size(); i ++) {
			Agent agent = (Agent)agentList.elementAt(i);
 			agent.step();
 			recordAgentExtractionLevel ( agent );	// incr number*Extractors and hhMonths* vars
 		}

		// do any poststep clean up for agents 
		for ( Agent a: agentList )
			a.postStep();
		
   		forestResource.postStep();

        // increase our tally of the total amount of extraction since the beginning of the simulation
		totalExtractionSinceBeginning = totalExtractionSinceBeginning + totalExtractionInPeriod;

		if ( getNumberTimesAtMinumumRes() == 1 )
			stepFirstHitMinResource = (int) getTickCount();

		stepReport();

	}

	/////////////////////////////////////////////////////////////////////////////////
	// recordAgentExtractionLevel
	// depending on the extractionlevel put it in one of the low,med,high categories
	// and increase the hhmonths variables
	public void recordAgentExtractionLevel ( Agent agent ) {
		if (agent.getScaledCurrentExtractionLevel() <= lowExtrThres){
			numberLowExtractors = numberLowExtractors + 1;     
			hhMonthsLow = hhMonthsLow + 1;
		}
		else if (agent.getScaledCurrentExtractionLevel() >= highExtrThres)
			{
 				numberHighExtractors = numberHighExtractors + 1;
				hhMonthsHigh = hhMonthsHigh + 1;
			}
		else {
			numberMedExtractors = numberMedExtractors + 1;
			hhMonthsMed = hhMonthsMed + 1;
		}
	}

	/////////////////////////////////////////////////////////////////////////////////
	// stepReport
	// each step write out: 
    //   timeStep  ...data...
	// currently, just writes the time step numbers!

	public void stepReport () {
		String s;
		if ( rDebug > 0 )
			System.out.printf( "==> Model stepReport %.0f:\n", getTickCount() );

		// set up a string with the values to write -- start with time step
	   	s = String.format( "%5.0f  ", getTickCount() );

		s += String.format( " %5.2f  %1.6f %3d  %3d  %3d   %3d", totalExtractionInPeriod,
							(getCurrentForestResourceQuantity()/initialResQuantity), getNumberTimesAtMinumumRes(),
							numberLowExtractors, numberMedExtractors, numberHighExtractors );

		s += String.format( " %6d %6d %6d", hhMonthsLow, hhMonthsMed, hhMonthsHigh );
		
		//UtilityMean	UtilityStdDev	ExtractionStdDev	DeviationMean	PenaltyMean
		s += String.format(" %f %f %f %f %f", getUtilityMean(), getUtilityStd(), getExtractionStd(), getDeviationMean(), getPenaltyMean()); 

		// write it to the xml and plain text report files
		//writeLineToReportFile ( "<stepreport>" + s + "</stepreport>" );
		writeLineToPlaintextReportFile( s );
		System.out.println(s);

		// flush the buffers so the data is not lost in a "crash"
		//getReportFile().flush();
		//getPlaintextReportFile().flush();

	}

	/////////////////////////////////////////////////////////////////////////////////
	// endReportFile
	// write a line
	//   # end time steps
	// and write a few ending numbers to the end of the report file
	// then call the super endReportFile to finish closing it, etc
	
	public void endReportFile () {
		String s;

		writeLineToPlaintextReportFile( "# end time steps" );

		// get the summary number of agents who spent most of their time
		// in each of the extraction categories. Store in the number*Extractors variables
		getNumberOfAgentsWhoSpentMajorityTimeInExtractionCategories();

		// calc the mean/sd of the amount extracted (abs) across agents.
		calcMeanStdDevOfAgentScaledTotalExtractionLevels();
		

		s = String.format( "# \n" );
		s += String.format( "# Summary measurements:\n" );

		s += String.format( "# stepFirstHitMinResource = %d\n", stepFirstHitMinResource );
		s += String.format( "# No. agents spent majority of steps low/med/high extractors:   %3d  %3d  %3d\n",
						numberLowExtractors, numberMedExtractors, numberHighExtractors	);
		s += String.format( "# Total abs overall (and per agent):       %9.2f   %7.2f \n",
							totalExtractionSinceBeginning,   
							totalExtractionSinceBeginning / agentList.size() ); 

		s += String.format( "# Total abs Extracted Per Agent (Mean/SD):  %9.2f   %7.2f\n",
						meanAgentAbsTotalExtractionAmount, sdAgentAbsTotalExtractionAmount );
		
		if (alphaSocEnvRatioStd1 > 0) {
			s += String.format("# For Cluster No 1 \n");
			double binCenter, binWidth;
			binWidth = alphaSocEnvRatioStd1/2.000;
			binCenter = alphaSocEnvRatioMean1 - 3.500*binWidth;
				while(binCenter <= alphaSocEnvRatioMean1 + 4.0*binWidth) {
					calcMeanStdDevOfAgentAbsTotalExtractionLevelsByASER(1 , binCenter, binWidth);
					s += String.format("# aSER = %1.2f ExtractionMean = %4.4f  ExtractionSD = %4.4f \n", 
							binCenter, meanAgentAbsTotalExtractionAmountByASERM, sdAgentAbsTotalExtractionAmountByASERM);
					binCenter += binWidth;
				}
		}
		if (alphaSocEnvRatioStd2 > 0) {
			s += String.format("# For Cluster No 2 \n");
			double binCenter, binWidth;
			binWidth = alphaSocEnvRatioStd2/2;
			binCenter = alphaSocEnvRatioMean2 - 3.5*binWidth;
				while(binCenter <= alphaSocEnvRatioMean2 + 4.0*binWidth) {
					calcMeanStdDevOfAgentAbsTotalExtractionLevelsByASER(2 , binCenter, binWidth);
					s += String.format("# aSER = %1.2f ExtractionMean = %4.4f  ExtractionSD = %4.4f \n", 
							binCenter, meanAgentAbsTotalExtractionAmountByASERM, sdAgentAbsTotalExtractionAmountByASERM);
					binCenter += binWidth;
				}
		}
		


		s += String.format( "# \n" );
		writeLineToPlaintextReportFile( s );

		super.endReportFile();

	}
	



	public void getNumberOfAgentsWhoSpentMajorityTimeInExtractionCategories() {
		int stepsLo, stepsMed, stepsHi;

		numberLowExtractors = numberMedExtractors = numberHighExtractors = 0;

		for ( Agent a: agentList ) {
			stepsLo = a.getTimeStepsAsLowExtractor();
			stepsMed = a.getTimeStepsAsMedExtractor();
			stepsHi = a.getTimeStepsAsHighExtractor();

			if ( stepsLo >= stepsMed ) {  	// must be low or high
				if ( stepsLo >= stepsHi )
					++numberLowExtractors;
				else
					++numberHighExtractors;
			} else  { 						// most be Med or High
				if ( stepsMed >= stepsHi )
					++numberMedExtractors;
				else
					++numberLowExtractors;
			}

		}

	}
	
	
	
	// calcMeanStdDevOfAgentAbsTotalExtractionLevels
	// loop over agents, asking each for its total absolute extraction amount.
	// Put into a common maths desc. stats object, then ask for mean/std
	// and store in Model iv's for printing in stepReport().
	//
	public void calcMeanStdDevOfAgentScaledTotalExtractionLevels() {
		DescriptiveStatistics avgDStats;
		avgDStats = DescriptiveStatistics.newInstance(); 

		for ( Agent a: agentList ) {
			avgDStats.addValue( a.getAbsTotalAmountExtracted() );
		}
		
		meanAgentAbsTotalExtractionAmount = avgDStats.getMean();
		sdAgentAbsTotalExtractionAmount = avgDStats.getStandardDeviation();
	}




	public void calcMeanStdDevOfAgentAbsTotalExtractionLevelsByASER(int clusterIndex, double binCenter, double binWidth) {
		DescriptiveStatistics extrByAserStats; 
		extrByAserStats = DescriptiveStatistics.newInstance();

		for ( Agent a: agentList ) {
			if (a.getClusterID()== clusterIndex) {
				if (a.getSocEnvRatio() > binCenter - binWidth/2 && a.getSocEnvRatio() < binCenter + binWidth/2) { 
					extrByAserStats.addValue( a.getAbsTotalAmountExtracted());
				}
			}

		}
		meanAgentAbsTotalExtractionAmountByASERM = extrByAserStats.getMean();
		sdAgentAbsTotalExtractionAmountByASERM = extrByAserStats.getStandardDeviation();
	}
		
		
		
		


	
	
	
	// writeHeaderCommentsToReportFile
	// customize to match what you are writing to the report files in stepReport.
	
	public void writeHeaderCommentsToReportFile () {
		writeLineToReportFile( "<comment>" );
		writeLineToReportFile( "         " );
		writeLineToReportFile( "</comment>" );

		
		writeLineToPlaintextReportFile( "# " );
		writeLineToPlaintextReportFile( "# " );
		writeLineToPlaintextReportFile( "# #Mn = number of times resource total <= minResVal." );
		writeLineToPlaintextReportFile( "# #hhMonths lo/med/hi -> cummulative #agents-steps in each" );
		writeLineToPlaintextReportFile( "# " );

		writeLineToPlaintextReportFile( "#             Resource       Number Agents        hhMonths	" );
		writeLineToPlaintextReportFile( "# time  Extr/T   Total #Mn  Low  Med  High    Low  Medium  High	UtilityMean	UtilityStdDev	ExtractionStdDev	DeviationMean	PenaltyMean" );
	}

	//////////////////////////////////////////////////////////////////////////////////
	// printProjectHelp
	// this could be filled in with some help to get from running with -help parameter
	
   	public void printProjectHelp() {
		// print project help

		System.out.printf( "\n%s -- Belief1 \n", getName() );

		System.out.printf( "\n" );

		System.out.printf( "Settable Parameterts:\n" );
		System.out.printf( "  SizeX, SizeY -- world size\n" );

		System.out.printf( "\n" );
		System.out.printf( "To compile:  bin/compile.sh\n" );
		System.out.printf( "or try this to get warnings:\n" );
		System.out.printf( "   bin/compile.sh -Xlint:unchecked\n" );
		System.out.printf( "Note there are a lot from ModelParameters.java that are ok.\n" );
		System.out.printf( "\n" );
		System.out.printf( "To run:\n" );
		System.out.printf( "   bin/guirun.sh\n" );
		System.out.printf( "   bin/batchrun.sh T=500 \n" );
		System.out.printf( "\n" );

		setParameterDefaultValues();

		printParametersMap();

		System.exit( 0 );

	}

	/////////////////////////////////////////////////////////////////////////////
	// processEndOfRun
	// called once, at end of run.
	// writes some final info, closes report files, etc.
	public void processEndOfRun ( ) {
		if ( rDebug > 0 )  
			System.out.printf("\n\n===== processEndOfRun =====\n\n" );
		applyAnyStoredChanges();
		endReportFile();

		this.fireStopSim();
	}

	//////////////////////////////////////////////////////////////////////////////
	public Schedule getSchedule () {	return schedule; }

	public String getName () { return "Model"; }

	// setters and getters
	// notes:
	// - we use the schedule != null to indicated model has been initialized
	// - some things can't be changed after model initialization
	//   (which things just depends on how the model is implemented)
	// - if we set something after model initialization,
	//   we need to write an change entry to the report file.
	// - some things need to send messages to update class variables.
	// 
	// NOTE: if you want changes a user makes to parameter like numBugs
	//       to be used after a restart (vs going back to defaults), 
	// you probably have to change setup() to not reinitialize IVs.

	public int getSizeX () { return sizeX; }
	public void setSizeX (int sizeX) { 
		this.sizeX = sizeX; 
		if (  schedule != null ) {
			System.err.printf("\nCan't change sizeX mid-run.\n");
			System.err.printf( "\nChange will not take effect until re-init.\n" );
		}
	}
	public int getSizeY () { return sizeY; }
	public void setSizeY (int sizeY) { 
		this.sizeY = sizeY;  
		if (  schedule != null ) {
			System.err.printf("\nCan't change sizeY mid-run.\n");
			System.err.printf( "\nChange will not take effect until re-init.\n" );
		}
	}

	public int getNumRandomOthers() { return numRandomOthers; }
	public void setNumRandomOthers( int val ) { numRandomOthers = val; }
	
	public int getSearchType(){ return searchType; }
	public void setSearchType(int val){  	searchType = val; }
  
	public double getProbPickRandomNeighbor () { return probPickRandomNeighbor; }
	public void setProbPickRandomNeighbor ( double val ) { probPickRandomNeighbor = val; }


  public int getOptimizer(){
	  	return optimizer;
	 }
	  
  public void setOptimizer(int val){
	  	optimizer = val;
   }
	  
  public int getCognitiveAbility(){
	  	return cognitiveAbility;
  }
		  
  public void setCognitiveAbility(int val){
	  	cognitiveAbility = val;
  }
		  
  public void setNormalSearch(int val){
	  	normalSearch = val;
 }
	  
public int getNormalSearch(){
	  	return normalSearch;
}
  
public double getNormalSearchStd(){
  	return normalSearchStd;
 }
  
public void setNormalSearchStd(double val){
	normalSearchStd = val;
}
  
  public Vector<Agent> getAgentList(){
  	return agentList;
  }
  
  public Resource getForestResource(){
  	return forestResource;
  }
  
  public int getPopulation(){
  	return population;
  }
  
  public double getResGrowthMean (){
  	return resGrowthMean;
  }
  
  public void setResGrowthMean(double val){
  	resGrowthMean = val;
  }
  
  public double getResGrowthStd(){
  	return resGrowthStd;
  }
  
  public void setResGrowthStd(double val){
  	resGrowthStd = val;
  }
  
  public double getHHSizeMean(){
  	return hhSizeMean;
  }
  
  public void setHHSizeMean(double val){
  	hhSizeMean = val;
  }
  
  public double getHHSizeStd(){
  	return hhSizeStd;
  }
  
  public void setHHSizeStd(double val){
  	hhSizeStd = val;
  }
  
  public double getLowExtrThres() {
  	return lowExtrThres;
  }
  public void setLowExtrThres(double val){
  	lowExtrThres = val;
  }
  
  public double getHighExtrThres() {
  	return highExtrThres;
  }
  public void setHighExtrThres(double val){
  	highExtrThres = val;
  }

  public double getIniScaledExtrLevMean() {
  	return iniScaledExtrLevMean;
  }
  public void setIniScaledExtrLevMean(double val){
  	iniScaledExtrLevMean = val;
  }
  
  public double getIniScaledExtrLevStd() {
  	return iniScaledExtrLevStd;
  }
  public void setIniScaledExtrLevStd(double val){
  	iniScaledExtrLevStd = val;
  }
  
  public int getIniExtrType() {
  	return iniExtrType;
  }
  public void setIniExtrType(int val){
  	iniExtrType = val;
  }
  
  public int getUseStringForAlphaMeans() {
	  	return useStringForAlphaMeans;
  }
  
  public void setUseStringForAlphaMeans(int val){
	  useStringForAlphaMeans = val;
  }

  public double getDecreasingReturnsFactor() {
	  	return decreasingReturnsFactor;
  }

  public void setDecreasingReturnsFactor(double val){
	  decreasingReturnsFactor = val;
  }

  public int getUseIdenticalStds() {
	  	return useIdenticalStds;
  }

public void setUseIdenticalStds(int val){
	useIdenticalStds = val;
}  
  
public int getUseIdenticalAserms() {
  	return useIdenticalAserms;
}

public void setUseIdenticalAserms(int val){
useIdenticalAserms = val;
}  

  public String getAlphaMeans1(){
	  	return alphaMeans1;
  }
  
  public void setAlphaMeans1(String val){
	 	alphaMeans1 = val;
  }
	  
  public String getAlphaMeans2(){
	  	return alphaMeans2;
  }

  public void setAlphaMeans2(String val){
	 	alphaMeans2 = val;
  }	  
  
  
  public double getAlphaConsumptionMean1(){
  	return alphaConsumptionMean1;
  }
  public void setAlphaConsumptionMean1(double val){
  	alphaConsumptionMean1 = val;
  }
  
  public double getAlphaLeisureMean1(){
  	return alphaLeisureMean1;
  }
  public void setAlphaLeisureMean1(double val){
  	alphaLeisureMean1 = val;
  }
  
  public double getAlphaNormSustMean1(){
  	return alphaNormSustMean1;
  }
  public void setAlphaNormSustMean1(double val){
  	alphaNormSustMean1 = val;
  }

  public double getAlphaSocEnvRatioMean1(){
	return alphaSocEnvRatioMean1;
  }
 public void setAlphaSocEnvRatioMean1(double val){
  	alphaSocEnvRatioMean1 = val;
  }
  
  public double getAlphaConsumptionStd1(){
 	return alphaConsumptionStd1;
  }
  public void setAlphaConsumptionStd1(double val){
  	alphaConsumptionStd1 = val;
  }
  
  public double getAlphaLeisureStd1(){
  	return alphaLeisureStd1;
  }
  public void setAlphaLeisureStd1(double val){
  	alphaLeisureStd1 = val;
  }
  
  public double getAlphaNormSustStd1(){
  	return alphaNormSustStd1;
  }
  public void setAlphaNormSustStd1(double val){
  	alphaNormSustStd1 = val;
  }
   
 public double getAlphaSocEnvRatioStd1(){
  	return alphaSocEnvRatioStd1;
  }

 public void setAlphaSocEnvRatioStd1(double val){
  	alphaSocEnvRatioStd1 = val;
  }


 public double getAlphaConsumptionMean2(){
	  	return alphaConsumptionMean2;
	  }
	  public void setAlphaConsumptionMean2(double val){
	  	alphaConsumptionMean2 = val;
	  }
	  
	  public double getAlphaLeisureMean2(){
	  	return alphaLeisureMean2;
	  }
	  public void setAlphaLeisureMean2(double val){
	  	alphaLeisureMean2 = val;
	  }
	  
	  public double getAlphaNormSustMean2(){
	  	return alphaNormSustMean2;
	  }
	  public void setAlphaNormSustMean2(double val){
	  	alphaNormSustMean2 = val;
	  }

	  public double getAlphaSocEnvRatioMean2(){
		return alphaSocEnvRatioMean2;
	  }
	 public void setAlphaSocEnvRatioMean2(double val){
	  	alphaSocEnvRatioMean2 = val;
	  }
	  
	  public double getAlphaConsumptionStd2(){
	 	return alphaConsumptionStd2;
	  }
	  public void setAlphaConsumptionStd2(double val){
	  	alphaConsumptionStd2 = val;
	  }
	  
	  public double getAlphaLeisureStd2(){
	  	return alphaLeisureStd2;
	  }
	  public void setAlphaLeisureStd2(double val){
	  	alphaLeisureStd2 = val;
	  }
	  
	  public double getAlphaNormSustStd2(){
	  	return alphaNormSustStd2;
	  }
	  public void setAlphaNormSustStd2(double val){
	  	alphaNormSustStd2 = val;
	  }
	   
	 public double getAlphaSocEnvRatioStd2(){
	  	return alphaSocEnvRatioStd2;
	  }

	 public void setAlphaSocEnvRatioStd2(double val){
	  	alphaSocEnvRatioStd2 = val;
	  } 
  
	 public double getMASTER_STD(){
	  	return MASTER_STD;
	  }

	 public void setMASTER_STD(double val){
	  	MASTER_STD = val;
	  } 

  public int getClusterType(){
  	return clusterType;
  }

  public void setClusterType(int val){
  	clusterType = val;
  }

  public double getClusterSize(){
	  	return clusterSize;
	  }

	  public void setClusterSize(double val){
	  	clusterSize = val;
	  }


  public double getInitialPerCapitaForestArea() {
  	return initialPerCapitaForestArea;
  }

  public void setInitialPerCapitaForestArea(double val){
  	initialPerCapitaForestArea = val;
  }


  public double getResMinimumValue(){
  	return resMinimumValue;
  }
  public void setResMinimumValue(double val){
  	resMinimumValue = val;
  }
  
  public double getCurrentForestResourceQuantity(){
	  return forestResource.getAbsoluteQuantity();
  }
  
  public double getGatheringTimeFactor(){
	  	return gatheringTimeFactor;
	  }
	  public void setGatheringTimeFactor(double val){
		  gatheringTimeFactor = val;
	  }

  public void setInstPerCapitaAllocation(int val){
      instPerCapitaAllocation = val;
  }
	  
  public int getInstPerCapitaAllocation(){
      return instPerCapitaAllocation;
  }
  
  public void setInstFrequency(int val){
	  instFrequency = val;
  }
		  
  public int getInstFrequency(){
	  return instFrequency;
  }	  

	
  public void setNoHigh(int val){
     noHigh = val;
  }
			  
  public int getNoHigh(){
	  return noHigh;
  }	    
	  
  public void setNoMed(int val){
	noMed = val;
  }
			  
  public int getNoMed(){
	return noMed;
  }	  	

  public void setNormType(int val){
	normType = val;
  }
				  
 public int getNormType(){
	return normType;
  }	  	
	
 
 public void setUseBurnInPeriod(int val){
	 useBurnInPeriod= val;
  }
					  
 public int getUseBurnInPeriod(){
	return useBurnInPeriod;
  }	  	
		
	
 public void incTotalExtractionInPeriod(double val){
	totalExtractionInPeriod = totalExtractionInPeriod + val;
 }

	// Called by resource when resource quantity goes down to minimum allowed.
	public void setResourceExtinct() {   
		resourceExtinct = 1;
		++numberTimesAtMinimumRes;
	}

	public int getNumberTimesAtMinumumRes() {
		return numberTimesAtMinimumRes;
	}
	
	/**
     * Calculate the average utility across households.
     *
     * @return
     */
    public double getUtilityMean() {
        double averageUtility = 0.0;

        for (Agent a : agentList) {
            averageUtility += a.getUtility(a.getScaledCurrentExtractionLevel());
        }

        return averageUtility / agentList.size();
    }

    /**
     * Calculate the standard deviation of utility across households.
     *
     * @return
     */
    public double getUtilityStd() {
        double averageUtility = getUtilityMean();
        double totalSS = 0.0;

        for (Agent a : agentList) {
            totalSS += Math.pow(averageUtility
                    - a.getUtility(a.getScaledCurrentExtractionLevel()), 2);
        }

        return totalSS / (agentList.size() - 1);
    }


    /**
     * Calculate the average amount of deviation from the institution's level.
     */
    public double getDeviationMean() {
        double difference = 0.0;

        for (Agent a : agentList) {
            difference += Math.abs(a.getScaledCurrentExtractionLevel()
                    - a.getScaledCurrentSustainabilityLevel());
        }

        return difference / agentList.size();
    }

    /**
     * Calculate the average extraction level across households.
     *
     * @return
     */
    public double getExtractionMean() {
        double averageExtraction = 0.0;

        for (Agent a : agentList) {
            averageExtraction += a.getScaledCurrentExtractionLevel();
        }

        return averageExtraction / agentList.size();
    }

    /**
     * Calculate the standard deviation of extraction level across households.
     *
     * @return
     */
    public double getExtractionStd() {
        double averageExtraction = getExtractionMean();

        double totalSS = 0.0;

        for (Agent a : agentList) {
            totalSS += Math.pow(averageExtraction
                    - a.getScaledCurrentExtractionLevel(), 2);

        }

        return totalSS / (agentList.size() - 1);
    }
    
    /**
     * Calculate the average penalty incurred by households.
     */
    public double getPenaltyMean() {
    	if (getHasPenalty() == 0) {
    		return 0.0;
    	}
    	
    	double averagePenalty = 0.0;
    	
    	for (Agent a: agentList) {
    		averagePenalty += getPenaltyCost(a.getScaledCurrentExtractionLevel());
    	}
    	
    	return averagePenalty /(float)agentList.size();
    }
}

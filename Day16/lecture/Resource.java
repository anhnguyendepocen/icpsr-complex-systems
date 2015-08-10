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
//import java.util.Vector;


import uchicago.src.sim.gui.Drawable;
import uchicago.src.sim.gui.SimGraphics;
//import uchicago.src.sim.util.Random;
import cern.jet.random.*;

public class Resource implements Drawable {

	public static double   resMinimumValue = 1.0;  // min value, to allow regrowth.

  	private Color agentShade;
  	private Model model;

  	private double absoluteQuantity;
  	private double actualGrowthQuantity;
  	private double headloadGatheringTime;
  	private double gatheringTimeFactor;
  	
	// private double initialForestArea;
//	private double initialForestDistance;
//	private double newForestDistance;
	private double initialAbsoluteResourceQuantity;
			
	private double resGrowthMean, resGrowthStd;

	/**
	 * Constructor: Initializes the resource object
	 * @param mod       Model
	 * @param quantity  initial amount of resource (absolute)
	 * @param mean      mean of normal growth rate
	 * @param std		std dev of normal growth rate
	 */
  	public Resource ( Model mod, double quantity, double mean, double std,  double gFactor){
    	
    	agentShade = new Color(1,0,0);
    	model = mod;
      	absoluteQuantity = Math.max( quantity, resMinimumValue );
      	actualGrowthQuantity = 0; 
		resGrowthMean = mean; 
		resGrowthStd = std;
	//	initialForestArea = area;
		initialAbsoluteResourceQuantity = absoluteQuantity;
		gatheringTimeFactor = gFactor;

  	}
  	
  	/**
  	 * What is written to the screen if the programmer asks that the agent be printed 
  	 * to the console. Here it is just the agents x,y coordinates.
  	 */
   	public String toString() {
		return "Resource";
  	}
	

	//////////////////////////////////////////////////////////////////
	// step
	// Called from Model.step for main resource dynamics
	// - make sure value is >= minimum
	// - grow at a rate in [0,1]
	//
  	public void step() {

		if ( absoluteQuantity < resMinimumValue ) {
			System.err.printf( "===> Resource-step() WARNING@T=%.0f: absoluteQ %.3f < resMinVal %.3f !!\n",
					model.getTickCount(), absoluteQuantity,  resMinimumValue );
			absoluteQuantity = resMinimumValue;
		}
		
		actualGrowthQuantity = -1.0;
  		while (actualGrowthQuantity<0) {
			actualGrowthQuantity =  Model.getNormalDoubleProb( resGrowthMean, resGrowthStd ) * absoluteQuantity;
		}
  		absoluteQuantity = absoluteQuantity + actualGrowthQuantity;	
//		newForestDistance = initialForestDistance + (1-absoluteQuantity/initialAbsoluteResourceQuantity) * Math.pow(initialForestArea, 0.5);
		headloadGatheringTime = 2 + Math.pow(initialAbsoluteResourceQuantity/absoluteQuantity, gatheringTimeFactor);
  	}
  	

	//////////////////////////////////////////////////////////////////
	// step
	// Called from Model.step at end, after all agents have acted.
	// Just checks to see if Resource hit the minimum this step.
	//
	public void postStep() {
		if ( absoluteQuantity <=  resMinimumValue ) {
			model.setResourceExtinct();
		}
	}

	///////////////////////////////////////////////////////////////////////////////
	// decrementResource
	// decrement resource by amount reguested, but don't go below min allowed.
	// tell model about and return amount actually decremented.
	//
  	public double decrementResource( double val ){
		double actualDecr;
		double origQ = absoluteQuantity;

  		if ( model.getRDebug() > 2 ) 
			System.out.printf( "-decrementResource by %.3f\n", absoluteQuantity );

		if ( val <= absoluteQuantity - resMinimumValue )
			actualDecr = val;
		else
			actualDecr = absoluteQuantity - resMinimumValue;

  		absoluteQuantity = absoluteQuantity - actualDecr;



		model.incTotalExtractionInPeriod( actualDecr );

		if ( absoluteQuantity < resMinimumValue ) {
			System.err.printf( "===> Resource-decrementResource() WARNING@T=%.0f: absoluteQ %f < resMinVal %.2f !!\n",
					model.getTickCount(), absoluteQuantity,  resMinimumValue );
			System.err.printf( "===> val=%f  actualDecr=%f  origQ=%f\n",
							   val, actualDecr, origQ );
	
		}

		return actualDecr;
   	}
  	
	///////////////////////////////////////////////////////////////////////////////
	// setters and getters


	public static void setResMinimumValue ( double d ) {
		resMinimumValue = d;
	}

  	public double getAbsoluteQuantity(){
  		return absoluteQuantity;
  	}
  	
  	public double getActualGrowthQuantity(){
  		return actualGrowthQuantity;
  	}

//  	public double getNewForestDistance(){
//  		return newForestDistance;
//  	}


  	public double getHeadloadGatheringTime(){
  		return headloadGatheringTime;
  	}  	
  	
	/* (non-Javadoc)
	 * @see uchicago.src.sim.gui.Drawable#getX()
	 */
	public int getX() {
		return 0;
	}

	/* (non-Javadoc)
	 * @see uchicago.src.sim.gui.Drawable#getY()
	 */
	public int getY() {
		return 0;
	}
	
  	public void draw(SimGraphics g) {
  		//int shade = (int)(belief*255);
  		//agentShade = new Color(shade,0,0);
  		//g.drawFastCircle ( agentShade);
  		agentShade = new Color(255,0,0);
  		g.drawFastCircle(agentShade);
  	}
}

package Belief2;

/**
 * Import statements: import libraries that aid in the development of this class
 * java. libraries are provided by java.
 * uchicago. libraries are part of the RePast simulation toolkit
 */
import java.awt.Color;
import java.util.Vector;

import cern.jet.random.Normal;
import cern.jet.random.engine.MersenneTwister;

import uchicago.src.sim.gui.Drawable;
import uchicago.src.sim.gui.SimGraphics;
import uchicago.src.sim.util.Random;

public class Institution implements Drawable{
		 
  	private Vector<Agent> popList;
  	private Color agentShade;
  	private Model model;
  	private Resource forestResource;
  	private double sustainableExtraction;
  	private int numberHouseholds;
  	private int instFrequency;
  	private int perCapitaAllocation;
  	private int clock;
  	private double growth;
  	
  
	/**
	 * Constructor: Initializes the agent location and space it will be located in
	 * @param mod model controling the step method of this institution
	 * @param instFreq frequency of the institution reporting value to the model or agentList *********************????
	 * @param perCap the type of sustainability calculation the institution is making  ************************????
	 */
  	public Institution (Model mod, int instFreq, int perCap){
    	popList = mod.getAgentList();
    	numberHouseholds = popList.size();
    	agentShade = new Color(1,0,0);
    	model = mod;
    	forestResource = model.getForestResource();
    	perCapitaAllocation = perCap;
    	instFrequency = instFreq;
    	clock = 0;
    	growth = 0;
  	}
  	
  	/**
  	 * What is written to the screen if the programmer asks that the agent be printed 
  	 * to the console. Here it is just the agents x,y coordinates.
  	 */
   	public String toString() {
		return "Institution";
  	}
	
	/////////////////////////////////////////////////////////////////////////
	/**
	 * step -- called each Model step(), implements institution dynamics
	 * 
	 * Every instFrequency steps recalculate a sustainability recommendation,
	 * and tell all the agents about it.
	 *    perCapitaAllocation == 1 -> every individual get an equal share
	 *    otherwise                -> every household gets equal share
	 * clock is how long since last recalculation
	 * Notes:
	 * - growth is the cummulative since last calculation, so when we 
	 *   do the calculation, divide by instFrequency (== to clock at that point).
	 */
  	public void step() {
  		clock = clock + 1;
  		growth = growth + forestResource.getActualGrowthQuantity();
  		
  		if ( clock >= instFrequency ){
			sustainableExtraction = Math.max( 0, growth/clock );
			for ( Agent a : popList ) {
  	  			if ( perCapitaAllocation==1 ) {
					int population = model.getPopulation();
					a.setAbsSustainabilityLevel( sustainableExtraction * a.getHHSize()/population );   // note that a.setAbsSustainabilityLevel() also internally calculates and sets a.scaledSustainabilityLevel
				}
  	  			else 
					a.setAbsSustainabilityLevel( sustainableExtraction / numberHouseholds );	// note that a.setAbsSustainabilityLevel() also internally calculates and sets a.scaledSustainabilityLevel
  	  		}  
  	  		growth = 0;
  	  		clock = 0;
  		}
  	}
  	
  	public double getSustainableExtraction(){
  		return sustainableExtraction;
  	}
  	
  	
  	public void draw(SimGraphics g) {
  		//int shade = (int)(belief*255);
  		//agentShade = new Color(shade,0,0);
  		//g.drawFastCircle ( agentShade);
  		agentShade = new Color(255,0,0);
  		g.drawFastCircle(agentShade);
  	}

	/* (non-Javadoc)
	 * @see uchicago.src.sim.gui.Drawable#getX()
	*/
  	public int getX() {
		// 
		return 0;
	}

	/* (non-Javadoc)
	 * @see uchicago.src.sim.gui.Drawable#getY()
	 */
	public int getY() {
		// 
		return 0;
	}
	
}

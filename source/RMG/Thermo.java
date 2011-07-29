////////////////////////////////////////////////////////////////////////////////
//
//	RMG - Reaction Mechanism Generator
//
//	Copyright (c) 2002-2011 Prof. William H. Green (whgreen@mit.edu) and the
//	RMG Team (rmg_dev@mit.edu)
//
//	Permission is hereby granted, free of charge, to any person obtaining a
//	copy of this software and associated documentation files (the "Software"),
//	to deal in the Software without restriction, including without limitation
//	the rights to use, copy, modify, merge, publish, distribute, sublicense,
//	and/or sell copies of the Software, and to permit persons to whom the
//	Software is furnished to do so, subject to the following conditions:
//
//	The above copyright notice and this permission notice shall be included in
//	all copies or substantial portions of the Software.
//
//	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//	DEALINGS IN THE SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////


import java.util.*;
import java.io.*;

import jing.chem.*;
import jing.chemParser.*;
import jing.param.*;
    import jing.chemUtil.*;
//import bondGroups.*;
import jing.rxn.*;
import jing.rxnSys.*;
import jing.mathTool.*;


public class Thermo {

  //## configuration RMG::RMG
public static void main(String[] args) {
//  initializeSystemProperties();

	RMG.globalInitializeSystemProperties();
        Logger.initialize();
        ReactionModelGenerator rmg = new ReactionModelGenerator();
	LinkedHashSet speciesSet = new LinkedHashSet();
    String thermo_output = "";
    Temperature systemTemp = new Temperature();
    LinkedHashMap speciesFromInputFile = new LinkedHashMap();

 try {
          FileReader in = new FileReader("thermo_input.txt");
          BufferedReader data = new BufferedReader(in);
          
          // Read the first line of thermo_input.txt
          String line = ChemParser.readMeaningfulLine(data, true);
          StringTokenizer st = new StringTokenizer(line);
          // The first line should start with "Solvation", otherwise do nothing and display a message to the user
          if (st.nextToken().startsWith("Solvation")) {
        	  line = st.nextToken().toLowerCase();
        	  // The options for the "Solvation" field are "on" or "off" (as of 18May2009), otherwise do nothing and display a message to the user
        	  // Note: I use "Species.useInChI" because the "Species.useSolvation" updates were not yet committed.
        	  if (line.equals("on")) {
                      rmg.setUseSolvation(true);
                      Species.useSolvation = true;
                      rmg.readAndMakePAL();
                    String solventname = st.nextToken().toLowerCase();
                    rmg.readAndMakeSL(solventname);

        	  } else if (line.equals("off")) {
        		  Species.useSolvation = false;
        		  thermo_output += "Gas-phase Thermochemistry.\n\n";
        	  } else {
        		  System.out.println("Error in reading thermo_input.txt file:\nThe field 'Solvation' has the options 'on' or 'off'.");
        		  return;
        	  }
              // Read in the temperature of the system
              line = ChemParser.readMeaningfulLine(data, true);
              st = new StringTokenizer(line);
              if (!st.nextToken().startsWith("Temperature")) {
                  System.out.println("Error in reading thermo_input.txt file:\n The field 'Temperature' should follow the 'Solvation' field.");
              }
              double tempValue = Double.parseDouble(st.nextToken());
              String tempUnits = st.nextToken();
              systemTemp = new Temperature(tempValue,tempUnits);
              thermo_output += "System Temperature = " + systemTemp.getK() + "K" + "\n";

            
            line = ChemParser.readMeaningfulLine(data, true);
            if (line.toLowerCase().startsWith("primarythermolibrary")) {
            	rmg.readAndMakePTL(data);
            }
            else {
            	System.err.println("ThermoDataEstimator: Could not locate the PrimaryThermoLibrary field." +
            			"Line read was: " + line);
            	System.exit(0);
            }

        	  // Read in the ChemGraphs and compute their thermo, while there are ChemGraphs to read in
        	  line = ChemParser.readMeaningfulLine(data, true);
        	  while (line != null) {
        		  String speciesName = line;
        		  Graph g = ChemParser.readChemGraph(data);
        		  ChemGraph cg = null;
        		  try {
        			  cg = ChemGraph.make(g);
        		  } catch (ForbiddenStructureException e) {
        			  System.out.println("Error in reading graph: Graph contains a forbidden structure.\n" + g.toString());
        			  System.exit(0);
        		  }

                          ReactionModelGenerator.addChemGraphToListIfNotPresent_ElseTerminate(speciesFromInputFile,cg,speciesName);

        		  Species species = Species.make(speciesName,cg);
        		  speciesSet.add(species);
        		  line = ChemParser.readMeaningfulLine(data, true);
        	  }
          } else
        	  System.out.println("Error in reading thermo_input.txt file:\nThe first line must read 'Solvation: on/off'.");

          in.close();

          thermo_output += "Output: Name" + "\t" + "deltaHsolv298 [=] kcal/mol"+ "\t" + "deltaSsolv298 [=] kcal/mol/K" + "\n" + "\n";
          
          Iterator iter = speciesSet.iterator();       
          while (iter.hasNext()){
        	  Species spe = (Species)iter.next();
 
            double deltaH0 = spe.getChemGraph().getdeltaHsolv298();
            double deltaS0 = spe.getChemGraph().getdeltaSsolv298();
            thermo_output += spe.getName() + "\t" + deltaH0 + "\t" + deltaS0 + "\n";

          }
          
          try {
        	  File thermoOutput = new File("thermo_output.txt");
        	  FileWriter fw = new FileWriter(thermoOutput);
        	  fw.write(thermo_output);
        	  fw.close();
          } catch (IOException e) {
        	  System.out.println("Error in writing thermo_output.txt file.");
        	  System.exit(0);
          }
 }
 catch (FileNotFoundException e) {
   System.err.println("File was not found!\n");
 }
 catch(IOException e){
   System.err.println("Something wrong with ChemParser.readChemGraph");
 }

System.out.println("Done!\n");

};

}

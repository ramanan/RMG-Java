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
import jing.rxnSys.*;

/**
 * Author: AJ
 * Function: Takes a list of NASA polynomials as input and returns user specified thermo data.
 *
 */
public class NASAThermoEstimator {

    public static void main(String[] args) {

        RMG.globalInitializeSystemProperties();
        try {
            File file = new File(args[0]);
            BufferedReader reader = new BufferedReader(new FileReader(file));
            readNASAFromFile(reader);
        } catch (InvalidChemGraphException e) {
            e.printStackTrace();
        } catch (ForbiddenStructureException e) {
            e.printStackTrace();
        } catch (IOException e) {
            System.out.println(e.toString());
        }
        Logger.info("Done!\n");
    }

    private static void readNASAFromFile(BufferedReader reader) throws IOException, ForbiddenStructureException {

        String line = ChemParser.readMeaningfulLine(reader, true);
        String ThermoBurcat = "";
        ThermoBurcat += "Values of 0.0 imply temperatures which were outside the recommended temperature range in the PrIMe data" + "\n" + "\n";
        ThermoBurcat += "Name" + "\t" + "G300" + "\t" + "G400" + "\t" + "G500" + "\t" + "G600" + "\t" + "G700" + "\t" + "G800" + "\t" + "G900" + "\t" + "G1000" + "\t" + "G1100" + "\t" + "G1200" + "\n" + "\n";
        String ThermoRMG = "";
        ThermoRMG += "Name" + "\t" + "G300" + "\t" + "G400" + "\t" + "G500" + "\t" + "G600" + "\t" + "G700" + "\t" + "G800" + "\t" + "G900" + "\t" + "G1000" + "\t" + "G1100" + "\t" + "G1200" + "\n" + "\n";
        String JTherGasInChIList = "";

        while (line != null) {
             if (line.toLowerCase().startsWith("inchi")){
//            if (line.startsWith("s000")) {
                String AdjList = "";
                String name = line;
//                String inchi = ChemParser.readMeaningfulLine(reader, true);
                String inchi = line;
                AdjList = Species.inchi2AdjList(inchi);
                Graph g = ChemParser.readAdjList(AdjList);
                ChemGraph cg = ChemGraph.make(g);
                Species spe = Species.make(inchi, cg);

                // Check if the current species has more than 4 carbon atoms. If yes the write it to an InChI list which will then be used as input to JTherGas
                if (spe.getChemGraph().getCarbonNumber() >= 4) {
                    JTherGasInChIList += inchi + "\t" + name + "\n";
                    //}
                    ThermoBurcat += inchi;
                    ThermoRMG += inchi;

                    String firstline = ChemParser.readMeaningfulLine(reader, true);
                    String secondline = ChemParser.readMeaningfulLine(reader, true);
                    String thirdline = ChemParser.readMeaningfulLine(reader, true);
                    String fourthline = ChemParser.readMeaningfulLine(reader, true);

                    String dataString = firstline + "\n" + secondline + "\n" + thirdline + "\n" + fourthline;
                    NASAThermoData nasaThermoData = null;
                    nasaThermoData = new NASAThermoData(dataString);

                    double lowT = Double.parseDouble(dataString.substring(48, 54));
                    Temperature Tlow = new Temperature(lowT, "K");
                    double highT = Double.parseDouble(dataString.substring(57, 64));
                    Temperature Thigh = new Temperature(highT, "K");

                    NASAThermoData RMGnasaThermoData = spe.getNasaThermoData();

                    double Trange[] = {300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200};
                    for (double temp : Trange) {
                        Temperature T = new Temperature(temp, "K");
                        double Gprime;
                        double Grmg;

                        // Check if the current temperature is within the valid range for the current species
                        if (T.getK() > Tlow.getK() && T.getK() < Thigh.getK()) {
                            Gprime = nasaThermoData.calculateFreeEnergy(T);
                        } else {
                            //Assign an arbitrary value of 0.0 to G values that lie outside the valid T range of the given species
                            System.out.println("Outside valid temperature range for species: " + name);
                            Gprime = 0.0;
                        }
                        Grmg = RMGnasaThermoData.calculateFreeEnergy(T);
                        ThermoBurcat += "\t" + Gprime;
                        ThermoRMG += "\t" + Grmg;
                    }
                    ThermoBurcat += "\n";
                    ThermoRMG += "\n";
                }
            }
            line = ChemParser.readMeaningfulLine(reader, true);
        }
        try {
            File primedat = new File("ThermoPrime.txt");
            FileWriter fw = new FileWriter(primedat);
            fw.write(ThermoBurcat);
            fw.close();
            System.out.println("Results written to ThermoPrime.txt");
        } catch (IOException e) {
            System.out.println("Could not write ThermoPrime.txt");
            System.exit(0);
        }
        try {
            File rmgdat = new File("ThermoRMG.txt");
            FileWriter fw = new FileWriter(rmgdat);
            fw.write(ThermoRMG);
            fw.close();
            System.out.println("Results written to ThermoRMG.txt");
        } catch (IOException e) {
            System.out.println("Could not write ThermoRMG.txt");
            System.exit(0);
        }
        try {
            File inchidat = new File("JTherGasInChI.txt");
            FileWriter fw = new FileWriter(inchidat);
            fw.write(JTherGasInChIList);
            fw.close();
            System.out.println("Results written to JTherGasInChI.txt");
        } catch (IOException e) {
            System.out.println("Could not write JTherGasInChI.txt");
            System.exit(0);
        }
    }
}

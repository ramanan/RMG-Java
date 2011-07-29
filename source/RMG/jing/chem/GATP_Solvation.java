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

package jing.chem;

import java.util.*;
import jing.param.*;
import jing.chemUtil.*;
import jing.param.*;
import jing.chemUtil.*;
import jing.mathTool.*;
import jing.param.Temperature;
import jing.rxnSys.*;

public class GATP_Solvation implements GeneralSolvationGAPP {

    private static GATP_Solvation INSTANCE = new GATP_Solvation();

    // Constructors
    private  GATP_Solvation() {
    }

    //Reading in solvent information
    ReactionModelGenerator rmg = new ReactionModelGenerator();
    SolventData solvent = ReactionModelGenerator.getSolvent();

    //Assigning values of commonly used constants
    double R=8.314;                                        // Gas constant units J/mol K
    double T=298;                                          // Standard state temperature
	
	public ThermoData generateSolvThermoData(ChemGraph p_chemGraph) {
        
            double deltaH0 = calculatedeltaHsolv(p_chemGraph);
            double deltaS0 = calculatedeltaSsolv(p_chemGraph);

            // Generation of Gas Phase data to add to the solution phase quantities
            ThermoData solvationCorrection = new ThermoData(deltaH0, deltaS0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,"Solvation correction");
            return solvationCorrection;
    }
	
	
	protected static GATP_Solvation getINSTANCE() {
        return INSTANCE;
    }

        public double calculatedeltaHsolv(ChemGraph p_chemGraph){

            AbramData descriptors = new AbramData();
            descriptors = p_chemGraph.getAbramData();

            double S = descriptors.S;
            double B = descriptors.B;
            double E = descriptors.E;
            double L = descriptors.L;
            double A = descriptors.A;

            double c_h = solvent.c_h;
            double e_h = solvent.e_h;
            double s_h = solvent.s_h;
            double a_h = solvent.a_h;
            double b_h = solvent.b_h;
            double l_h = solvent.l_h;

            double deltaH0 = c_h + s_h*S + b_h*B + e_h*E + l_h*L + a_h*A;    // Implementation of Mintz model for calculation of solution phase enthalpy (kJ/mol)
            deltaH0 = deltaH0/4.18;                                            // Conversion from kJ/mol to kcal/mol
            return deltaH0;
        }

        public double calculatedeltaGsolv(ChemGraph p_chemGraph){

            AbramData descriptors = new AbramData();
            descriptors = p_chemGraph.getAbramData();

            double S = descriptors.S;
            double B = descriptors.B;
            double E = descriptors.E;
            double L = descriptors.L;
            double A = descriptors.A;

            double c_g = solvent.c_g;
            double e_g = solvent.e_g;
            double s_g = solvent.s_g;
            double a_g = solvent.a_g;
            double b_g = solvent.b_g;
            double l_g = solvent.l_g;

            double logK = c_g + s_g*S + b_g*B + e_g*E + l_g*L + a_g*A;    // Implementation of Abraham Model for calculation of partition coefficient
            double deltaG0 = -8.314*298*2.303*logK;                       // J/mol
            deltaG0 = deltaG0/4180;                                       // conversion from kJ/mol to kcal/mol
            return deltaG0;
        }

        public double calculatedeltaSsolv(ChemGraph p_chemGraph){

            double deltaH0 = calculatedeltaHsolv(p_chemGraph);
            double deltaG0 = calculatedeltaGsolv(p_chemGraph);
            double deltaS0 = (deltaH0 - deltaG0)/298;                         // This quantity is in terms of kcal/mol/K
            deltaS0 = deltaS0*1000;                                             // conversion from kcal/mol/K to cal/mol/K which are the default RMG units
            return deltaS0;
        }
       
}

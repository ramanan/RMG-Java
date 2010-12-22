/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jing.chem;
import java.io.*;
import java.util.*;

import jing.chemUtil.*;
import jing.chemParser.*;

/**
 *
 * @author amrit
 */
public class SolventData {


//Declaration of solvent attributes
public String name;
public String comments;

public double c_g;
public double s_g;
public double e_g;
public double l_g;
public double a_g;
public double b_g;


public double c_h;
public double s_h;
public double e_h;
public double l_h;
public double a_h;
public double b_h;


public double viscosity;

//Constructors
	public  SolventData() {
        super();
		}   

    public  SolventData(double p_cg, double p_eg, double p_sg, double p_ag, double p_bg, double p_lg, double p_ch, double p_eh, double p_sh, double p_ah, double p_bh, double p_lh) {

        c_g = p_cg;
        s_g = p_sg;
        e_g = p_eg;
        l_g = p_lg;
        a_g = p_ag;
        b_g = p_bg;


        c_h = p_ch;
        s_h = p_sh;
        e_h = p_eh;
        l_h = p_lh;
        a_h = p_ah;
        b_h = p_bh;
       
		}

    public  SolventData(String p_name, SolventData p_solvent, String p_comments) {

        c_g = p_solvent.c_g;
        s_g = p_solvent.s_g;
        e_g = p_solvent.e_g;
        l_g = p_solvent.l_g;
        a_g = p_solvent.a_g;
        b_g = p_solvent.b_g;

        c_h = p_solvent.c_h;
        s_h = p_solvent.s_h;
        e_h = p_solvent.e_h;
        l_h = p_solvent.l_h;
        a_h = p_solvent.a_h;
        b_h = p_solvent.b_h;

        comments = p_comments;
        name = p_name;


		}

    public  SolventData(String p_name) {

        name = p_name;

    }




protected double getc_g() {
       return c_g;
    }

protected double gets_g() {
       return s_g;
    }

protected double gete_g() {
       return e_g;
    }

protected double getl_g() {
       return l_g;
    }

protected double geta_g() {
       return a_g;
    }

protected double getb_g() {
       return b_g;
    }

protected double getc_h() {
       return c_h;
    }

protected double gets_h() {
       return s_h;
    }

protected double geta_h() {
       return a_h;
    }

protected double getb_h() {
       return b_h;
    }

protected double getl_h() {
       return l_h;
    }

protected String getName() {
       return name;
    }

}

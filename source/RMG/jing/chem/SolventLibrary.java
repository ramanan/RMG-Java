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
public class SolventLibrary {

protected static HashMap library;

	public SolventLibrary(String name, String location) {
		library = new HashMap();
		readSolventLibrary(name, location);
	}


//Reading in the Solvent library
	public void readSolventLibrary(String name, String path) {

		String solventLibraryDirectory = System.getProperty("jing.chem.SolventLibrary.pathName");
		String libraryFile = solventLibraryDirectory + "/" + path + "/Library.txt";
		try {
			read(libraryFile, name);
		}
		catch (IOException e) {
			System.err.println("RMG cannot read Solvent Library: " + name + "\n" + e.getMessage());
		}
	}

	public void read(String p_library, String p_name) throws IOException, FileNotFoundException {
		String source = "Solvent Library: " + p_name;
		System.out.println("Reading " + source);
		library = readLibrary(p_library, source);
	}

	public HashMap readLibrary(String p_transportFileName, String source) throws IOException {
		try{
			FileReader in = new FileReader(p_transportFileName);
			BufferedReader data = new BufferedReader(in);
			HashMap tempLibrary = new HashMap();
			String line = ChemParser.readMeaningfulLine(data, true);
			while (line != null){
				StringTokenizer token = new StringTokenizer(line);
				String name = token.nextToken();
				String abram = token.nextToken();
				try{
					SolventData solventData = new SolventData(
							Double.parseDouble(abram),
							Double.parseDouble(token.nextToken()),
							Double.parseDouble(token.nextToken()),
							Double.parseDouble(token.nextToken()),
							Double.parseDouble(token.nextToken()),
                            Double.parseDouble(token.nextToken()),
                            Double.parseDouble(token.nextToken()),
                            Double.parseDouble(token.nextToken()),
                            Double.parseDouble(token.nextToken()),
                            Double.parseDouble(token.nextToken()),
                            Double.parseDouble(token.nextToken()),
                            Double.parseDouble(token.nextToken())
                            );
					String comments = "";
					while (token.hasMoreTokens()) {
						comments = comments + " " + token.nextToken();
					}
					SolventData newSolventData = new SolventData(name, solventData, comments);
                    tempLibrary.put(name, newSolventData);
					library.put(name, newSolventData);

				}
				catch (NumberFormatException e) {
				}
				line = ChemParser.readMeaningfulLine(data, true);
			}
			in.close();
		    return library;
		}
		catch (Exception e){
		   throw new IOException("Can't read solvent descriptors in solvent library!");
		}
	}

	public static SolventData getSolventData(String solv_name){
		if (library == null) return null;

		SolventData solv = (SolventData)library.get(solv_name);
		if (solv != null){
			return solv;
		}
		Iterator iter = library.keySet().iterator();
		while (iter.hasNext()){
			String name = (String)iter.next();
			if (name.equals(solv_name)){
				solv = (SolventData)library.get(name);
				return solv;
			}
		}
        if (solv==null){
        System.out.println("Solvent does not exist in RMG solvent library");
        System.out.println("Stopping.");
		throw new Error();
        }
        return null;
	}

}

//!********************************************************************************
//!
//!    RMG: Reaction Mechanism Generator
//!
//!    Copyright: Jing Song, MIT, 2002, all rights reserved
//!
//!    Author's Contact: jingsong@mit.edu
//!
//!    Restrictions:
//!    (1) RMG is only for non-commercial distribution; commercial usage
//!        must require other written permission.
//!    (2) Redistributions of RMG must retain the above copyright
//!        notice, this list of conditions and the following disclaimer.
//!    (3) The end-user documentation included with the redistribution,
//!        if any, must include the following acknowledgment:
//!        "This product includes software RMG developed by Jing Song, MIT."
//!        Alternately, this acknowledgment may appear in the software itself,
//!        if and wherever such third-party acknowledgments normally appear.
//!
//!    RMG IS PROVIDED "AS IS" AND ANY EXPRESSED OR IMPLIED
//!    WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//!    OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//!    DISCLAIMED.  IN NO EVENT SHALL JING SONG BE LIABLE FOR
//!    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
//!    OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
//!    OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//!    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//!    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
//!    THE USE OF RMG, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//!
//!******************************************************************************



package jing.chem;


import java.io.*;
import java.util.*;
import jing.chemUtil.*;
import jing.chemParser.*;
import jing.chemUtil.Arc;
import jing.chemUtil.Node;
import jing.chemUtil.Graph;
import jing.param.Temperature;

//## package jing::chem

//----------------------------------------------------------------------------
// jing\chem\ChemGraph.java
//----------------------------------------------------------------------------

//## class ChemGraph
public class ChemGraph implements Matchable {

    protected int MAX_OXYGEN_NUM = 6;		//## attribute MAX_OXYGEN_NUM
	protected int MAX_CARBON_NUM = 8;       //SS

	/**
    Maximal radical number allowed in a ChemGraph.
    */
    protected static int MAX_RADICAL_NUM = 3;		//## attribute MAX_RADICAL_NUM

    /**
    Chemical Formula of a ChemGraph.
    */
    protected String chemicalFormula = null;		//## attribute chemicalFormula

    /**
    The overall forbidden structure.  When any new ChemGraph instance is generated, RMG check if it has any of the forbidden structure.  If it has, it wont be generated.
    */
    protected static HashSet forbiddenStructure = new HashSet();		//## attribute forbiddenStructure

    protected int internalRotor = 0;		//## attribute internalRotor

    /**
    A collection of all the possible symmetry Axis in a ChemGraph.
    For example: the C=C=O skeleton in (CH3)2C=C=O
    */
    protected HashSet symmetryAxis = null;		//## attribute symmetryAxis

    /**
    Symmetry number of a ChemGraph.  Used in calculating entropy.
    */
    protected int symmetryNumber = -1;		//## attribute symmetryNumber

    /**
    It is a unique string representation for this chem structure.  The method to generating it has not implemented yet.
    */
    protected String uniqueString;		//## attribute uniqueString

    protected Graph graph;
    protected Species species;
    protected ThermoData thermoData;
    protected GeneralGAPP thermoGAPP;

    // Constructors

    //## operation ChemGraph()
    private  ChemGraph() {
        //#[ operation ChemGraph()
        //#]
    }
    //## operation ChemGraph(Graph)
    private  ChemGraph(Graph p_graph) throws ForbiddenStructureException {
        //#[ operation ChemGraph(Graph)
        graph = p_graph;

        if (isForbiddenStructure(p_graph) || getRadicalNumber() > MAX_RADICAL_NUM || getOxygenNumber() > MAX_OXYGEN_NUM) {
        	graph = null;
        	throw new ForbiddenStructureException(p_graph.toString());
        }


        //#]
    }

    /**
    Requires:
    Effects: saturate all the node atom's undefined valence by adding Hydrogen.
    Modifies: this.graph.nodeList
    */
    //## operation addMissingHydrogen()
    public void addMissingHydrogen() {
        //#[ operation addMissingHydrogen()
        Atom H = Atom.make(ChemElement.make("H"), FreeElectron.make("0"));
        Bond S = Bond.make("S");
        HashMap addedH = new HashMap();

        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	Atom atom = (Atom)node.getElement();
        	int val = (int)atom.getValency();

        	double bondOrder = 0;
        	Iterator neighbor_iter = node.getNeighbor();
        	while (neighbor_iter.hasNext()) {
        		Arc arc = (Arc)neighbor_iter.next();
        		Bond bond = (Bond)arc.getElement();
        		bondOrder += bond.getOrder();
        	}
//        	if (bondOrder > val) throw new InvalidConnectivityException();
//        	else if (bondOrder < val) {
//        		addedH.put(node, new Integer(val-bondOrder));
//        	}
        	if (bondOrder < val) {
        		addedH.put(node, new Integer(val-(int)(bondOrder+1.0e-8)));
        	}
        }
        Graph g = getGraph();
        iter = addedH.keySet().iterator();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	int Hnum = ((Integer)addedH.get(node)).intValue();
        	for (int i=0;i<Hnum; i++) {
        		Node newNode = g.addNode(H);
        		g.addArcBetween(node, S, newNode);
        	}
        	node.updateFgElement();
        }

        return;
        //#]
    }

    public static ChemGraph saturate(ChemGraph p_chemGraph) {
       //#[ operation saturate(ChemGraph)
       int max_radNum_molecule = ChemGraph.getMAX_RADICAL_NUM();
       int max_radNum_atom = Math.min(8,max_radNum_molecule);

       ChemGraph result = null;
       try {
               result = ChemGraph.copy(p_chemGraph);
       }
       catch (Exception e) {
               System.out.println(e.getMessage());
               System.exit(0);
       }

       FreeElectron satuated = FreeElectron.make("0");
       Atom H = Atom.make(ChemElement.make("H"),satuated);
       Bond S = Bond.make("S");
       Graph g = result.getGraph();
       int nn = g.getHighestNodeID();
       for (int i = 0 ; i < nn; i++) {
               Node node = g.getNodeAt(i);
               if (node != null) {
                       Atom atom = (Atom)node.getElement();
                       int HNum = atom.getRadicalNumber();
                       if (atom.isRadical()) {
                               Atom newAtom = new Atom(atom.getChemElement(),satuated);
                               node.setElement(newAtom);
                               node.updateFeElement();
                               for (int j = 0; j < HNum; j++) {
                                       Node n = g.addNode(H);
                                       g.addArcBetween(node,S,n);
                               }
                               node.updateFgElement();
                       }
               }
       }

       return result;
       /*
       int max_radNum_molecule = ChemGraph.getMAX_RADICAL_NUM();
       int max_radNum_atom = Math.min(8,max_radNum_molecule);
       int [] idArray = new int[max_radNum_molecule];
       Atom []  atomArray = new Atom[max_radNum_molecule];
       Node [][] newnode = new Node[max_radNum_molecule][max_radNum_atom];

       int radicalSite = 0;
       Iterator iter = p_chemGraph.getNodeList();
       FreeElectron satuated = FreeElectron.make("0");
       while (iter.hasNext()) {
               Node node = (Node)iter.next();
               Atom atom = (Atom)node.getElement();
               if (atom.isRadical()) {
                       radicalSite ++;
                       // save the old radical atom
                       idArray[radicalSite-1] = node.getID().intValue();
                       atomArray[radicalSite-1] = atom;
                       // new a satuated atom and replace the old one
                       Atom newAtom = new Atom(atom.getChemElement(),satuated);
                       node.setElement(newAtom);
                       node.updateFeElement();
               }
       }

       // add H to satuate chem graph
       Atom H = Atom.make(ChemElement.make("H"),satuated);
       Bond S = Bond.make("S");
       for (int i=0;i<radicalSite;i++) {
               Node node = p_chemGraph.getNodeAt(idArray[i]);
               Atom atom = atomArray[i];
               int HNum = atom.getRadicalNumber();
               for (int j=0;j<HNum;j++) {
                       newnode[i][j] = g.addNode(H);
                       g.addArcBetween(node,S,newnode[i][j]);
               }
               node.updateFgElement();
       }
       */
       //#]
   }

    /**
    Requires: acyclic ChemGraph
    Effects: calculate and return the symmetry number centered at p_node atom
    Modifies:
    */
    //## operation calculateAtomSymmetryNumber(Node)
    public int calculateAtomSymmetryNumber(Node p_node) {
        //#[ operation calculateAtomSymmetryNumber(Node)
        // note: acyclic structure!!!!!!!!!!!!!
        int sn = 1;

        // if no neighbor or only one neighbor, sigma = 1, return 1;
        int neighborNumber = p_node.getNeighborNumber();
        if (neighborNumber < 2) return sn;

        Atom atom = (Atom)p_node.getElement();
        Iterator neighbor_iter = p_node.getNeighbor();
        FGElement fge = (FGElement)p_node.getFgElement();

        if (!atom.isRadical()) {
        // satuated atom symmetric number calculation
        	if (fge.equals(FGElement.make("Cs"))) {
        	// Cs:
               	Arc a1 = (Arc)neighbor_iter.next();
               	Arc a2 = (Arc)neighbor_iter.next();
               	Arc a3 = (Arc)neighbor_iter.next();
               	Arc a4 = (Arc)neighbor_iter.next();
               	if (p_node.isSymmetric(a1,a2)) {
               		if (p_node.isSymmetric(a1,a3)) {
          				if (p_node.isSymmetric(a1,a4)) {
        				// AAAA
        					sn *= 12;
               			}
               			else {
               			// AAAB
         					sn *= 3;
               			}
               		}
               		else {
          				if (p_node.isSymmetric(a1,a4)) {
        				// AAAB
               				sn *= 3;
               			}
               			else if (p_node.isSymmetric(a3,a4)) {
               			// AABB
               				sn *= 2;
               			}
               			else {
        				// AABC
               			}
               		}
               	}
               	else {
               		if (p_node.isSymmetric(a1,a3)) {
               			if (p_node.isSymmetric(a1,a4)) {
               			// AAAB
               				sn *= 3;
               			}
               			else if (p_node.isSymmetric(a2,a4)) {
               			// AABB
               				sn *= 2;
               			}
               			else {
               			// AABC
               			}
               		}
               		else if (p_node.isSymmetric(a2,a3)) {
               			if (p_node.isSymmetric(a1,a4)) {
               			// AABB
               				sn *= 2;
               			}
               			else if (p_node.isSymmetric(a2,a4)) {
               			// AAAB
               				sn *= 3;
               			}
               			else {
               			// AABC
               			}
               		}
               		else {
               			// AABC or ABCD
               		}
               	}
        	}
        	else if (fge.equals(FGElement.make("Os"))) {
        	// Os:
               	Arc a1 = (Arc)neighbor_iter.next();
               	Arc a2 = (Arc)neighbor_iter.next();
               	if (p_node.isSymmetric(a1,a2)) {
               		sn *= 2;
               	}
        	}
        	else if (fge.equals(FGElement.make("Cdd"))) {
        	// Cdd:
               	Arc a1 = (Arc)neighbor_iter.next();
               	Arc a2 = (Arc)neighbor_iter.next();
               	if (p_node.isSymmetric(a1,a2)) {
               		sn *= 2;
                }
          	}
        }
        else {
        // radical symmetric number calculation
        	if (fge.equals(FGElement.make("Cs"))) {
        	// only consider Cs. and Cs..
        		FreeElectron fe = atom.getFreeElectron();
            	if (fe.getOrder() == 1) {
        		// mono-radical Cs.
                	Arc a1 = (Arc)neighbor_iter.next();
                	Arc a2 = (Arc)neighbor_iter.next();
                	Arc a3 = (Arc)neighbor_iter.next();
                	if (p_node.isSymmetric(a1,a2)) {
                		if (p_node.isSymmetric(a1,a3))
                			sn *= 6;
                		else
                			sn *= 2;
                	}
                	else {
                		if (p_node.isSymmetric(a1,a3) || p_node.isSymmetric(a2,a3))
                			sn *= 2;
                	}
        		}
        		else if (fe.getOrder() == 2) {
        		// bi-radical Cs..
                	Arc a1 = (Arc)neighbor_iter.next();
                	Arc a2 = (Arc)neighbor_iter.next();
                	if (p_node.isSymmetric(a1,a2))
                		sn *= 2;
        		}
        	}
        }

        return sn;
        //#]
    }

    /**
    Requires: acyclic ChemGraph
    Effects: calculate and return the symmetry number by all the possible symmetry axis in this ChemGraph.
    Modifies:
    */
    //## operation calculateAxisSymmetryNumber()
    public int calculateAxisSymmetryNumber() {
        //#[ operation calculateAxisSymmetryNumber()
        int sn = 1;
        // note: acyclic structure!!!!!!!!!!!!!
        symmetryAxis = new HashSet();

        Iterator iter = getArcList();
        while (iter.hasNext()) {
        	Arc arc = (Arc)iter.next();
        	Bond bond = (Bond)arc.getElement();
        	if (bond.isDouble()) {
        		Iterator neighbor_iter = arc.getNeighbor();
        		Node n1 = (Node)neighbor_iter.next();
        		Node n2 = (Node)neighbor_iter.next();

        		FGElement Cd = FGElement.make("Cd");
        		FGElement Cdd = FGElement.make("Cdd");

        		FGElement fge1 = (FGElement)n1.getFgElement();
        		FGElement fge2 = (FGElement)n2.getFgElement();

                HashSet axis = new HashSet();
                axis.add(arc);
                if (fge1.equals(Cdd)) n1 = getToEndOfAxis(arc,n1,axis);
                if (fge2.equals(Cdd)) n2 = getToEndOfAxis(arc,n2,axis);

               	Atom atom1 = (Atom)n1.getElement();
               	Atom atom2 = (Atom)n2.getElement();

               	if (atom1.isRadical() || atom2.isRadical()) return sn;

               	Bond D = Bond.make("D");

               	if (!symmetryAxis.contains(axis)) {
               		symmetryAxis.add(axis);
               		boolean l1 = n1.isLeaf();
               		boolean l2 = n2.isLeaf();
               		if (!l1 && !l2) {
        				Iterator i1 = n1.getNeighbor();
        				Iterator i2 = n2.getNeighbor();
        				Arc a1 = (Arc)i1.next();
        				if (((Bond)a1.getElement()).equals(D)) a1 = (Arc)i1.next();
        				Arc a2 = (Arc)i1.next();
        				if (((Bond)a2.getElement()).equals(D)) a2 = (Arc)i1.next();

        				Arc a3 = (Arc)i2.next();
        				if (((Bond)a3.getElement()).equals(D)) a3 = (Arc)i2.next();
        				Arc a4 = (Arc)i2.next();
        				if (((Bond)a4.getElement()).equals(D)) a4 = (Arc)i2.next();

        				if (n1.isSymmetric(a1,a2) && n2.isSymmetric(a3,a4)) {
        					sn *= 2;
        				}
        			}
        			else if (!l1 && l2) {
        				Iterator i = n1.getNeighbor();
        				Arc a1 = (Arc)i.next();
        				if (((Bond)a1.getElement()).equals(D)) a1 = (Arc)i.next();
        				Arc a2 = (Arc)i.next();
        				if (((Bond)a2.getElement()).equals(D)) a2 = (Arc)i.next();

        				if (n1.isSymmetric(a1,a2)) {
        					sn *= 2;
        				}
        			}
        			else if (l1 && !l2) {
        				Iterator i = n2.getNeighbor();
        				Arc a1 = (Arc)i.next();
        				if (((Bond)a1.getElement()).equals(D)) a1 = (Arc)i.next();
        				Arc a2 = (Arc)i.next();
        				if (((Bond)a2.getElement()).equals(D)) a2 = (Arc)i.next();

        				if (n2.isSymmetric(a1,a2)) {
        					sn *= 2;
        				}
        			}
        		}
        	}
        }

        return sn;
        //#]
    }

    /**
    Requires: acyclic ChemGraph
    Effects: calculate and return the symmetry number centered at p_arc bond
    Modifies:
    */
    //## operation calculateBondSymmetryNumber(Arc)
    public int calculateBondSymmetryNumber(Arc p_arc) {
        //#[ operation calculateBondSymmetryNumber(Arc)
        // note: acyclic structure!!!!!!!!!!!!!

        int sn = 1;

        // if no neighbor or only one neighbor, sigma = 1, return 1;
        int neighborNumber = p_arc.getNeighborNumber();
        if (neighborNumber != 2) throw new InvalidNeighborException("bond has " + neighborNumber + " neighbor!");

        Bond bond = (Bond)p_arc.getElement();
        Iterator neighbor_iter = p_arc.getNeighbor();
        Node n1 = (Node)neighbor_iter.next();
        Node n2 = (Node)neighbor_iter.next();

        if (bond.isSingle() || bond.isDouble() || bond.isTriple()) {
        	if (p_arc.isSymmetric(n1,n2)) {
        		boolean opt = checkOpticalIsomer(p_arc);
        		if (!opt) sn *= 2;
        	}
        }

        return sn;
        //#]
    }

    /**
    Requires:
    Effects: return Cp(T)
    Modifies:
    */
    //## operation calculateCp(Temperature)
    public double calculateCp(Temperature p_temperature) {
        //#[ operation calculateCp(Temperature)
        return getThermoData().calculateCp(p_temperature);
        //#]
    }

    /**
Requires:
Effects: calculate and return the symmetry number of the cyclic portion of this ChemGraph
Modifies:
svp
*/
//## operation calculateCyclicSymmetryNumber()
public int calculateCyclicSymmetryNumber(){
  //#[ operation calculateCyclicSymmetryNumber()
 int sn = 1;
 LinkedList ring_structures = new LinkedList();//list of ring structures
 LinkedList cycle_list = new LinkedList();//list of cycles
 Iterator cycle_iter = getGraph().getCycle();
 while (cycle_iter.hasNext()){
   LinkedList current_cycle = (LinkedList)cycle_iter.next();
   cycle_list.add(current_cycle);
 }
 //if 2 graph components share at least one cycle, they belong to the same ring structure
 for (int i = 0; i <= cycle_list.size()-1; i++){
   LinkedList current_ring = (LinkedList)cycle_list.get(i);
   if (ring_structures.isEmpty()){
     ring_structures.add(current_ring);
   }
   else{
     int same_gc = 0;
     Iterator ring_structure_iter = ring_structures.iterator();
     Iterator current_ring_iter = current_ring.iterator();
     while (ring_structure_iter.hasNext()){
       LinkedList ring_structure = (LinkedList) ring_structure_iter.next();
       while (current_ring_iter.hasNext()) {
         GraphComponent gc = (GraphComponent) current_ring_iter.next();
         if (ring_structure.contains(gc)) {
           Iterator current_iter = current_ring.iterator();
           while (current_iter.hasNext()){
             GraphComponent current_gc = (GraphComponent)current_iter.next();
             if (!ring_structure.contains(current_gc)){
               ring_structure.add(current_gc);
               ring_structures.set(ring_structures.indexOf(ring_structure), ring_structure);
             }
           }
           same_gc++;
         }
       }
     }
     if (same_gc == 0){
       ring_structures.add(current_ring);
     }
   }
   if (i != cycle_list.size()-1){
     for (int j = 1; j <= cycle_list.size()-1; j++){
       LinkedList next_cycle = (LinkedList)cycle_list.get(j);
       Iterator ring_structures_iter = ring_structures.iterator();
       Iterator next_cycle_iter = next_cycle.iterator();
       while (ring_structures_iter.hasNext()){
         LinkedList current_ring_structure = (LinkedList)ring_structures_iter.next();
         while (next_cycle_iter.hasNext()){
           GraphComponent gc = (GraphComponent)next_cycle_iter.next();
           if (current_ring_structure.contains(gc)){
             Iterator ring_iter = next_cycle.iterator();
             while(ring_iter.hasNext()){
               GraphComponent current_gc = (GraphComponent)ring_iter.next();
               if (!current_ring_structure.contains(current_gc)){
                 current_ring_structure.add(current_gc);
                 ring_structures.set(ring_structures.indexOf(current_ring_structure),current_ring_structure);
               }
             }
             break;
           }
         }
       }
     }
   }
 }
 Iterator iter = ring_structures.iterator();
 while (iter.hasNext()){
   LinkedList current_ring_structure = (LinkedList)iter.next();
   Iterator gc_iter = current_ring_structure.iterator();
   LinkedList node_list = new LinkedList(); //list of all cyclic nodes
   LinkedList arc_list = new LinkedList(); //list of all cyclic arcs
   while (gc_iter.hasNext()){
     GraphComponent gc = (GraphComponent)gc_iter.next();
     gc.setVisited(false);
     if (gc instanceof Node){
       node_list.add(gc);
     }
     else {
       arc_list.add(gc);
     }
   }
   //find all sets of equal nodes
   LinkedList equal_node_list = new LinkedList();
   for (int i = 0; i <= node_list.size() - 2; i++) {
     Node current = (Node) node_list.get(i);
     if (!current.isVisited()) {
       LinkedList equal_list = new LinkedList(); //list of equivalent nodes
       current.setVisited(true);
       equal_list.add(current);
       for (int j = i + 1; j <= node_list.size() - 1; j++) {
         Node next = (Node) node_list.get(j);
         Iterator list_iter = equal_list.iterator();
         while (list_iter.hasNext()) {
           current = (Node) list_iter.next();
           if (isSymmetric(current, next)) {
             equal_list.add(next);
             next.setVisited(true);
             break;
           }
         }
       }
       equal_node_list.add(equal_list); //add list of equivalent nodes to list of sets of equivalent nodes
     }
   }
   //find all sets of equal arcs
   LinkedList equal_arc_list = new LinkedList();
   for (int i = 0; i <= arc_list.size() - 2; i++) {
     Arc current = (Arc) arc_list.get(i);
     if (!current.isVisited()) {
       LinkedList equal_list = new LinkedList(); //list of equivalent arcs
       current.setVisited(true);
       equal_list.add(current);
       for (int j = i + 1; j <= arc_list.size() - 1; j++) {
         Arc next = (Arc) arc_list.get(j);
         Iterator list_iter = equal_list.iterator();
         while (list_iter.hasNext()) {
           current = (Arc) list_iter.next();
           if (isSymmetric(current, next)) {
             equal_list.add(next);
             next.setVisited(true);
             break;
           }
         }
       }
       equal_arc_list.add(equal_list); //add list of equivalent arcs to list of sets of equivalent arcs
     }
   }
   //find largest set of equal nodes
   int node_sn = 1;
   Iterator node_list_iter = equal_node_list.iterator();
   while (node_list_iter.hasNext()) {
     LinkedList current = (LinkedList) node_list_iter.next();
     if (current.size() > node_sn) {
       node_sn = current.size(); //node symmetry number = size of largest set of equivalent nodes
     }
   }
   //find largest set of equal arcs
   int arc_sn = 1;
   Iterator arc_list_iter = equal_arc_list.iterator();
   while (arc_list_iter.hasNext()) {
     LinkedList current = (LinkedList) arc_list_iter.next();
     if (current.size() > arc_sn) {
       arc_sn = current.size(); //arc symmetry number = size of largest set of equivalent arcs
     }
   }
   if (node_sn == node_list.size() && arc_sn == arc_list.size()) { //all nodes equal and all arcs equal
     sn *= node_sn;
     sn *= 2;
     Node first_node = (Node)node_list.getFirst();
     FGElement fge = (FGElement)first_node.getFgElement();
     if (fge.equals(FGElement.make("Cs"))){
       LinkedList acyclic_neighbor = new LinkedList();
       Iterator neighbor_iter = first_node.getNeighbor();
       while (neighbor_iter.hasNext()){
         Arc arc = (Arc)neighbor_iter.next();
         if (!arc.getInCycle()){
           acyclic_neighbor.add(arc);
         }
       }
       if (acyclic_neighbor.size() == 2){
         Arc a1 = (Arc) acyclic_neighbor.getFirst();
         Arc a2 = (Arc) acyclic_neighbor.getLast();
         if (!first_node.isSymmetric(a1, a2)) {
           sn /= 2;
         }
       }
     }
   }
   else {
     if (node_sn >= arc_sn) {
       sn *= node_sn;
     }
     else {
       sn *= arc_sn;
     }
   }
   //if (sn >= 2 && sn%2 == 0){//added by Sally for non-planar PAH's
     //sn = correctSymmetryNumber(sn);
   //}
 }
return sn;
//#]
}


    /**
    Requires:
    Effects: return G(T)
    Modifies:
    */
    //## operation calculateG(Temperature)
    public double calculateG(Temperature p_temperature) {
        //#[ operation calculateG(Temperature)
        return getThermoData().calculateG(p_temperature);
        //#]
    }

    /**
    Requires:
    Effects:return H(T)
    Modifies:
    */
    //## operation calculateH(Temperature)
    public double calculateH(Temperature p_temperature) {
        //#[ operation calculateH(Temperature)
        return getThermoData().calculateH(p_temperature);
        //#]
    }

    //## operation calculateInternalRotor()
    public void calculateInternalRotor() {
        //#[ operation calculateInternalRotor()
        // add more check for resonance axis!!
        int rotor = 0;
        Graph g = getGraph();
        for (Iterator iter = g.getArcList(); iter.hasNext();) {
        	Arc a = (Arc)iter.next();
        	Bond bond = (Bond)a.getElement();
        	if (bond.isSingle() && !a.getInCycle()) {
        		Iterator atomIter = a.getNeighbor();
        		Node n1 = (Node)atomIter.next();
        		Node n2 = (Node)atomIter.next();
        		if (!n1.isLeaf() && !n2.isLeaf()) {
        			rotor++;
        		}
        	}
        }
        internalRotor = rotor;

        //#]
    }

//	## operation isLinear() 
    public boolean isLinear() {
        //#[ operation isLinear() 
        // only check for linearity in molecules with at least two atoms
        if (getAtomNumber() == 1) return false;
        
        // cyclic molecules are not linear
        if (!isAcyclic()) return false;
        
        // biatomic molecules are always linear
        if (getAtomNumber() == 2) return true;
        
        // molecules with only double bonds are linear (e.g. CO2)
        boolean allDouble = true;
        Iterator iter = getArcList();
        while (iter.hasNext()) {
        	Arc arc = (Arc)iter.next();
        	Bond bond = (Bond)arc.getElement();
        	if (!bond.isDouble()) allDouble = false;
        }
        if (allDouble) return true;
        
        // molecule with alternating single and triple bonds are linear (e.g. acetylene)
        boolean alternatingSingleTriple = true;
        Iterator node_iter = getNodeList();
        while (node_iter.hasNext()) {
        	Node node = (Node)node_iter.next();
        	int neighborNumber = node.getNeighborNumber();
        	if (neighborNumber == 2) {
        		Iterator neighbor_iter = node.getNeighbor();
        		Arc a1 = (Arc)neighbor_iter.next();
        		Bond b1 = (Bond)a1.getElement();
        		Arc a2 = (Arc)neighbor_iter.next();
        		Bond b2 = (Bond)a2.getElement();
        		if (! ((b1.isTriple() && b2.isSingle()) || (b1.isSingle() && b2.isTriple())))
        			alternatingSingleTriple = false;
        	}
        	else if (neighborNumber > 2)
        		alternatingSingleTriple = false;
        }
            
            if (alternatingSingleTriple) return true;
        
        // if none of the above are true, it's nonlinear
        return false;
        
        //#]
    }
	
    /**
    Requires:
    Effects: return S(T)
    Modifies:
    */
    //## operation calculateS(Temperature)
    public double calculateS(Temperature p_temperature) {
        //#[ operation calculateS(Temperature)
        return getThermoData().calculateS(p_temperature);




        //#]
    }

    //## operation calculateSymmetryNumber()
   public int calculateSymmetryNumber() {
       //#[ operation calculateSymmetryNumber()
       try {
         getGraph().identifyCycle();//svp
               int sn = 1;
               Iterator iter = getNodeList();
               while (iter.hasNext()) {
                       Node node = (Node)iter.next();
                       if (!node.getInCycle()){//svp
                        sn *= calculateAtomSymmetryNumber(node);
                        }

               }
               iter = getArcList();
               while (iter.hasNext()) {
                       Arc arc = (Arc)iter.next();
                       if (!arc.getInCycle()){//svp
                       sn *= calculateBondSymmetryNumber(arc);
                     }

               }
               sn *= calculateAxisSymmetryNumber();
               if (!isAcyclic()) {//svp
                sn *= calculateCyclicSymmetryNumber();
              }

               symmetryNumber = sn;
               return sn;
       }
       catch (ClassCastException e) {
               throw new InvalidChemGraphException();
       }

       //#]
   }

    /**
    Requies:
    Effects: check if p_arc bond is the single bond between two Oxygens, which is considered as optical isomers
    Modifies:
    */
    //## operation checkOpticalIsomer(Arc)
    public boolean checkOpticalIsomer(Arc p_arc) {
        //#[ operation checkOpticalIsomer(Arc)
        // check if the p_arc is -O-O-
        Bond b = (Bond)p_arc.getElement();
        if (b.isSingle()) {
        	Iterator neighbor_iter = p_arc.getNeighbor();
        	Node n1 = (Node)neighbor_iter.next();
        	Node n2 = (Node)neighbor_iter.next();
        	Atom a1 = (Atom)n1.getElement();
        	Atom a2 = (Atom)n2.getElement();
        	FGElement fge1 = (FGElement)n1.getFgElement();
        	FGElement fge2 = (FGElement)n2.getFgElement();
        	if (!a1.isRadical() && fge1.equals(FGElement.make("Os")) && !a2.isRadical() && fge2.equals(FGElement.make("Os"))) {
        		return true;
        	}
        }

        return false;


        //#]
    }

    /**
    Requires:
    Effects: clear the central node list
    Modifies:
    */
    //## operation clearCentralNode()
    public void clearCentralNode() {
        //#[ operation clearCentralNode()
        getGraph().clearCentralNode();
        //#]
    }

    /**
    Requires:
    Effects: return a new instance identical to this ChemGraph
    Modifies:
    */
    //## operation copy(ChemGraph)
    public static ChemGraph copy(ChemGraph p_chemGraph) throws ForbiddenStructureException {
        //#[ operation copy(ChemGraph)
        try {
        Graph g = Graph.copy(p_chemGraph.getGraph());

        ChemGraph cg = new ChemGraph(g);
        cg.uniqueString = p_chemGraph.getUniqueString();
        cg.chemicalFormula = p_chemGraph.getChemicalFormula();
        cg.species = p_chemGraph.getSpecies();
        cg.symmetryNumber = p_chemGraph.symmetryNumber;
        cg.thermoData = p_chemGraph.thermoData;
        cg.thermoGAPP = p_chemGraph.thermoGAPP;

        HashSet oldSymmetryAxis = p_chemGraph.getSymmetryAxis();
        if (oldSymmetryAxis != null) {
        	cg.symmetryAxis = new HashSet();
        	for (Iterator iAxis = oldSymmetryAxis.iterator(); iAxis.hasNext(); ) {
        		HashSet newAxis = new HashSet();
        		HashSet oldAxis = (HashSet)iAxis.next();
        		for (Iterator iArc = oldAxis.iterator(); iArc.hasNext(); ) {
        			Arc arc = (Arc)iArc.next();
        			Iterator iNode = arc.getNeighbor();
         			int n1 = ((Node)iNode.next()).getID().intValue();
           			int n2 = ((Node)iNode.next()).getID().intValue();
           			Arc newArc = cg.getArcBetween(n1,n2);
           			newAxis.add(newArc);
           		}
        		cg.symmetryAxis.add(newAxis);
           	}
        }

        return cg;
        }
        catch (ForbiddenStructureException e) {
        	throw new ForbiddenStructureException(e.getMessage());
        }
        //#]
    }

    /**
    Requires:
    Effects: return true iff two chemgraph have equivalent graph structures
    Modifies:
    */
    //## operation equals(Object)
    public boolean equals(Object p_chemGraph) {
        //#[ operation equals(Object)
        if (this == p_chemGraph) return true;

        return isEquivalent((ChemGraph)p_chemGraph);
        //#]
    }

    /**
    Requires:
    Effects: generate the chemical formula of this chem graph and return it.  if the graph is not initialized, return null.
    Modifies:
    */
    //## operation generateChemicalFormula()
    public String generateChemicalFormula() {
        //#[ operation generateChemicalFormula()
        if (getGraph() == null) return null;

        /*int cap = ChemElementDictionary.size();
        Vector type = new Vector(cap);
        int[] number = new int[cap];
        */
        int C_number = 0;
        int H_number = 0;
        int O_number = 0;
        int radical = 0;

        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	Atom atom = (Atom)node.getElement();
        	radical += atom.getRadicalNumber();

        	if (atom.isCarbon()) {
        		C_number++;
        	}
        	else if (atom.isHydrogen()) {
        		H_number++;
        	}
        	else if (atom.isOxygen()) {
        		O_number++;
        	}
        	else {
        		throw new InvalidChemNodeElementException();
        	}
        }

        String s = "";
        if (C_number>0) {
        	s = s + "C";
        	if (C_number >1) {
        		s = s + String.valueOf(C_number);
        	}
        }
        if (H_number>0) {
        	s = s + "H";
        	if (H_number >1) {
        		s = s + String.valueOf(H_number);
        	}
        }
        if (O_number>0) {
        	s = s + "O";
        	if (O_number >1) {
        		s = s + String.valueOf(O_number);
        	}
        }

        chemicalFormula = s;
        if (radical == 1) {
        	chemicalFormula = chemicalFormula + ".";
        }
        else if (radical == 2) {
        	chemicalFormula = chemicalFormula + "..";
        }
        else if (radical == 3) {
        	chemicalFormula = chemicalFormula + "...";
        }

        return chemicalFormula;


        //#]
    }

    /**
    Requires:
    Effects: if the thermoGAPP is not set, set default GAPP.  Use it to calculate the thermoData of this chem graph.  if there is any exception during this process, throw FailGenerateThermoDataException.
    Modifies: this.thermoData
    */
    //## operation generateThermoData()
    public ThermoData generateThermoData() throws FailGenerateThermoDataException {
        //#[ operation generateThermoData()
        // use GAPP to generate Thermo data
        try {
        	if (thermoGAPP == null) setDefaultThermoGAPP();
        	thermoData = thermoGAPP.generateThermoData(this);
        	return thermoData;
        }
        catch (Exception e) {
        	throw new FailGenerateThermoDataException();
        }
        //#]
    }

    /**
    Requires:
    Effects: return the Arc between two positions in this ChemGraph
    Modifies:
    */
    //## operation getArcBetween(int,int)
    public Arc getArcBetween(int p_position1, int p_position2) {
        //#[ operation getArcBetween(int,int)
        return getGraph().getArcBetween(p_position1,p_position2);
        //#]
    }

    /**
    Requires:
    Effects: return an arc iterator of this ChemGraph
    Modifies:
    */
    //## operation getArcList()
    public Iterator getArcList() {
        //#[ operation getArcList()
        return getGraph().getArcList();
        //#]
    }

    /**
    Requires:
    Effects: return the atom at the p_position in this chem graph; if p_position is empty, return null;
    Modifies:
    */
    //## operation getAtomAt(int)
    public Atom getAtomAt(int p_position) throws EmptyAtomException {
        //#[ operation getAtomAt(int)
        try {
        	return (Atom)(getNodeAt(p_position).getElement());
        }
        catch (NotInGraphException e) {
        	return null;
        }



        //#]
    }

    /**
    Requires:
    Effects: return the total atom number in this ChemGraph.
    Modifies:
    */
    //## operation getAtomNumber()
    public int getAtomNumber() {
        //#[ operation getAtomNumber()
        return getGraph().getNodeNumber();
        //#]
    }

    /**
    Requires:
    Effects: if there is a bond connecting p_position1 and p_position2, return that bond; otherwise, return null.
    Modifies:
    */
    //## operation getBondBetween(int,int)
    public Bond getBondBetween(int p_position1, int p_position2) throws EmptyAtomException {
        //#[ operation getBondBetween(int,int)
        try {
        	return (Bond)(getArcBetween(p_position1,p_position2).getElement());
        }
        catch (ClassCastException e) {
        	return null;
        }
        //#]
    }

    //## operation getCarbonNumber()
    public int getCarbonNumber() {
        //#[ operation getCarbonNumber()
        int cNum = 0;
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	Atom atom = (Atom)node.getElement();

        	if (atom.isCarbon()) {
        		cNum++;
        	}
        }
        return cNum;
        //#]
    }

    /**
    Requires:
    Effects: return the hashMap of centralNode in this ChemGraph
    Modifies:
    */
    //## operation getCentralNode()
    public HashMap getCentralNode() {
        //#[ operation getCentralNode()
        return getGraph().getCentralNode();
        //#]
    }

    /**
    Requires:
    Effects: return the node whose centralID equals p_position
    Modifies:
    */
    //## operation getCentralNodeAt(int)
    public Node getCentralNodeAt(int p_position) {
        //#[ operation getCentralNodeAt(int)
        return getGraph().getCentralNodeAt(p_position);


        //#]
    }

    /**
    Requires:
    Effects: return the number of the central nodes in this chem graph
    Modifies:
    */
    //## operation getCentralNodeNumber()
    public int getCentralNodeNumber() {
        //#[ operation getCentralNodeNumber()
        return getGraph().getCentralNodeNumber();
        //#]
    }

    //## operation getChemicalFormula()
    public String getChemicalFormula() {
        //#[ operation getChemicalFormula()
        if (chemicalFormula == null || chemicalFormula.length() == 0) generateChemicalFormula();
        return chemicalFormula;
        //#]
    }

    /**
    Requires:
    Effects: return the iterator loop over the graph's cycle list.
    Modifies:
    */
    //## operation getCycle()
    public Iterator getCycle() {
        //#[ operation getCycle()
        return getGraph().getCycle();
        //#]
    }

    //## operation getHydrogenNumber()
    public int getHydrogenNumber() {
        //#[ operation getHydrogenNumber()
        int hNum = 0;
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	Atom atom = (Atom)node.getElement();

        	if (atom.isHydrogen()) {
        		hNum++;
        	}
        }
        return hNum;
        //#]
    }

    /**
    Requires:
    Effects: add the weight of all the atoms in this graph to calculate the molecular weight of this chem graph
    Modifies:
    */
    //## operation getMolecularWeight()
    public double getMolecularWeight() {
        //#[ operation getMolecularWeight()
        double MW = 0;
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	MW += ((Atom)(node.getElement())).getWeight();
        }
        return MW;
        //#]
    }

    /**
    Requires:
    Effects: return the "legal" name of this ChemGraph.  The order of priority
    (1) uniqueString
    (2) species.name
    (3) chemicalFormula
    Modifies:
    */
    //## operation getName()
    public String getName() {
        //#[ operation getName()
        if (uniqueString != null && uniqueString.length() > 0) return uniqueString;

        String name = getSpecies().getName();
        if (name != null && name.length() > 0) return name;

        return getChemicalFormula();
        //#]
    }

    /**
    Requires:
    Effects: return the node whose ID equals p_position.
    Modifies:
    */
    //## operation getNodeAt(int)
    public Node getNodeAt(int p_position) {
        //#[ operation getNodeAt(int)
        return getGraph().getNodeAt(p_position);
        //#]
    }

    /**
    Requires:
    Effects: return the node whose ID equals p_ID.
    Modifies:
    */
    //## operation getNodeAt(Integer)
    public Node getNodeAt(Integer p_ID) {
        //#[ operation getNodeAt(Integer)
        return getGraph().getNodeAt(p_ID);
        //#]
    }

    /**
    Requires:
    Effects: return an iterator over the node collection of this graph
    Modifies:
    */
    //## operation getNodeList()
    public Iterator getNodeList() {
        //#[ operation getNodeList()
        return getGraph().getNodeList();
        //#]
    }

    //## operation getOxygenNumber()
    public int getOxygenNumber() {
        //#[ operation getOxygenNumber()
        int oNum = 0;
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node node = (Node)iter.next();
        	Atom atom = (Atom)node.getElement();

        	if (atom.isOxygen()) {
        		oNum++;
        	}
        }
        return oNum;
        //#]
    }


    /**
    Requires:
    Effects: return a collection of all the radical sites
    Modifies:
    */
    //## operation getRadicalNode()
    public HashSet getRadicalNode() {
        //#[ operation getRadicalNode()
        HashSet radicalNode = new HashSet();
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Node n = (Node)iter.next();
        	Atom a = (Atom)n.getElement();
        	if (a.isRadical()) {
        		radicalNode.add(n);
        	}
        }
        return radicalNode;
        //#]
    }

    /**
    Requires:
    Effects: calculate the total radical number in this chem graph.
    Modifies:
    */
    //## operation getRadicalNumber()
    public int getRadicalNumber() {
        //#[ operation getRadicalNumber()
        int radicalNumber = 0;
        Iterator iter = getNodeList();
        while (iter.hasNext()) {
        	Object element = ((Node)(iter.next())).getElement();
        	radicalNumber += ((Atom)element).getRadicalNumber();
        }
        return radicalNumber;
        //#]
    }

    //## operation getSymmetryNumber()
    public int getSymmetryNumber() {
        //#[ operation getSymmetryNumber()
        if (symmetryNumber < 0) calculateSymmetryNumber();
        return symmetryNumber;
        //#]
    }

    //## operation getThermoData()
    public ThermoData getThermoData() {
        //#[ operation getThermoData()
        if (thermoData == null) generateThermoData();
        return thermoData;
        //#]
    }

    /**
    Requires:
    Effects: find out the end of C=C=C... pattern
    Modifies:
    */
    //## operation getToEndOfAxis(Arc,Node,HashSet)
    private static final Node getToEndOfAxis(Arc p_beginArc, Node p_beginNode, HashSet p_axis) {
        //#[ operation getToEndOfAxis(Arc,Node,HashSet)
        Arc nextArc = null;
        Iterator iter = p_beginNode.getNeighbor();
        while (iter.hasNext()) {
        	nextArc = (Arc)iter.next();
        	if (nextArc != p_beginArc) break;
        }

        p_axis.add(nextArc);

        Node nextNode = nextArc.getOtherNode(p_beginNode);
        FGElement fge = (FGElement)nextNode.getFgElement();
        FGElement Cdd = FGElement.make("Cdd");

        if (!fge.equals(Cdd)) return nextNode;
        else {
        	return getToEndOfAxis(nextArc,nextNode,p_axis);
        }
        //#]
    }

    /**
    Requires:
    Effects: check if the element of every node is an atom, and if the element of every node is a bond.
    Modifies:
    */
    //## operation graphContentsOk(Graph)
    public static boolean graphContentsOk(Graph p_graph) {
        //#[ operation graphContentsOk(Graph)
        Iterator iter = p_graph.getNodeList();
        while (iter.hasNext()) {
        	Object atom = ((Node)iter.next()).getElement();
        	if (!(atom instanceof Atom)) return false;
        }
        iter = p_graph.getArcList();
        while (iter.hasNext()) {
        	Object bond = ((Arc)iter.next()).getElement();
        	if (!(bond instanceof Bond)) return false;
        }
        return true;
        //#]
    }

    /**
    Requires:
    Effects: return chemicalFormula's hashcode.  i.e., all the isomers have the same hashcode
    Modifies:
    */
    //## operation hashCode()
    public int hashCode() {
        //#[ operation hashCode()
        if (chemicalFormula == null) generateChemicalFormula();
        return chemicalFormula.hashCode();
        //#]
    }

    /**
    Requires:
    Effects: check all the possible reacted sites in this chemgraph according to the pass-in functional group or functional group collection.  If there are any matches, return all the matches in a linked list; otheriwse, return an empty list.
    Modifies:
    */
    //## operation identifyReactionMatchedSite(Matchable)
    public HashSet identifyReactionMatchedSite(Matchable p_functionalGroup) {
        //#[ operation identifyReactionMatchedSite(Matchable)
        if (p_functionalGroup instanceof FunctionalGroup) {
        	FunctionalGroup fg = (FunctionalGroup)p_functionalGroup;
        	boolean thisIsRadical = this.isRadical();
        	boolean fgIsRadical = fg.isRadical();
        	if (thisIsRadical == fgIsRadical) {
        		return getGraph().identifyAllOrderedMatchedSites(fg.getGraph());
        	}
        	else {
        		return new HashSet();
        	}
        }
        else if (p_functionalGroup instanceof FunctionalGroupCollection) {
        	HashSet result = new HashSet();
        	FunctionalGroupCollection fgc = (FunctionalGroupCollection)p_functionalGroup;
        	Iterator iter = fgc.getFunctionalGroups();
        	while (iter.hasNext()) {
        		FunctionalGroup fg = (FunctionalGroup)iter.next();
        		boolean thisIsRadical = this.isRadical();
        		boolean fgIsRadical = fg.isRadical();
        		if (thisIsRadical == fgIsRadical) {
           			HashSet site = getGraph().identifyAllOrderedMatchedSites(fg.getGraph());
        			result.addAll(site);
        		}
        	}
        	return result;
        }
        else {
        	throw new InvalidFunctionalGroupException();
        }

        //#]
    }

    /**
    Requires:
    Effects: get a collection of all the thermo matched site, this for corrections in GAPP.
    Modifies:
    */
    //## operation identifyThermoMatchedSite(FunctionalGroup)
    public HashSet identifyThermoMatchedSite(FunctionalGroup p_functionalGroup) {
        //#[ operation identifyThermoMatchedSite(FunctionalGroup)
        return getGraph().identifyAllUnorderedMatchedSite(p_functionalGroup.getGraph());
        //#]
    }

    /**
    Requires:
    Effects: return true if there is no cycle in the chem graph
    Modifies:
    */
    //## operation isAcyclic()
    public boolean isAcyclic() {
        //#[ operation isAcyclic()
        return getGraph().isAcyclic();
        //#]
    }

    /**
    Requires:
    Effects: return true iff this and p_chemGraph are equivalent chemgraphs.
    Modifies:
    */
    //## operation isEquivalent(ChemGraph)
    public boolean isEquivalent(ChemGraph p_chemGraph) {
        //#[ operation isEquivalent(ChemGraph)
        if (chemicalFormula == null) generateChemicalFormula();
        if (p_chemGraph.chemicalFormula == null) p_chemGraph.generateChemicalFormula();

        if (!getChemicalFormula().equals(p_chemGraph.getChemicalFormula())) return false;

        if (!getGraph().isEquivalent(p_chemGraph.getGraph())) return false;

        return true;
        //#]
    }

    /**
    Requires:
    Effects: return true iff this chemGraph contains forbidden structure.
    Modifies:
    */
    //## operation isForbiddenStructure(Graph)
    public static boolean isForbiddenStructure(Graph p_graph) {
        //#[ operation isForbiddenStructure(Graph)
        for (Iterator iter = forbiddenStructure.iterator(); iter.hasNext(); ) {
        	FunctionalGroup fg = (FunctionalGroup)iter.next();
        	Graph g = fg.getGraph();
        	if (p_graph.isSub(g)) {
        		return true;
        	}
        }
        return false;
        //#]
    }

    /**
    Requires:
    Effects: return true iff this chemgraph contains radical site
    Modifies:
    */
    //## operation isRadical()
    public boolean isRadical() {
        //#[ operation isRadical()
        return (getRadicalNumber() > 0);
        //#]
    }

    /**
    Requires:
    Effects: if p_functionalGroup is a FunctionalGroup, return if this chemgraph is matched with it at the central nodes; if p_functionalGroup is a FunctionalGroupCollection, return if this chemgraph is matched with any of the functionalgroup in the collection at the central node.  for all other case, return false.
    Modifies:
    */
    //## operation isSubAtCentralNodes(Matchable)
    public boolean isSubAtCentralNodes(Matchable p_functional) {
        //#[ operation isSubAtCentralNodes(Matchable)
        if (this == p_functional) return false;
        if (p_functional instanceof FunctionalGroup) {
        	return isSubAtCentralNodes((FunctionalGroup)p_functional);
        }
        else if (p_functional instanceof FunctionalGroupCollection) {
        	Iterator iter = ((FunctionalGroupCollection)p_functional).getFunctionalGroups();
        	while (iter.hasNext()) {
        		FunctionalGroup fg = (FunctionalGroup)iter.next();
        		if (isSubAtCentralNodes(fg)) return true;
          	}
          	return false;
        }
        else {
        	return false;
        }



        //#]
    }

    /**
   Requires: both graph components belong to the same graph
   Effects: return true if two graph components are equal
   Modifies:
   svp
   */
   //## operation isSymmetric(GraphComponent, GraphComponent)
   public boolean isSymmetric(GraphComponent p_gc1, GraphComponent p_gc2) {
   //#[ operation isSymmetric(GraphComponent, GraphComponent)
   Stack s1 = new Stack();
   Stack s2 = new Stack();
   if (p_gc1.isEquivalent(p_gc2, s1, s2)) {
     resetStack(s1);
       resetStack(s2);
       getGraph().resetMatchedGC();
       return true;
     }
     else {
       resetStack(s1);
       resetStack(s2);
       getGraph().resetMatchedGC();
       return false;
     }
   //#]
 }


    /**
    Requires:
    Effects: return if this chem graph is matched with p_functionalGroup at the central nodes.
    Modifies:
    */
    //## operation isSubAtCentralNodes(FunctionalGroup)
    public boolean isSubAtCentralNodes(FunctionalGroup p_functionalGroup) {
        //#[ operation isSubAtCentralNodes(FunctionalGroup)
        return getGraph().isSubAtCentralNodes(p_functionalGroup.getGraph());



        //#]
    }

    /**
    Requires:
    Effects: factory method for ChemGraph. make a new instance of ChemGraph.  Do such things for it:
    (1) add missing H
    (2) generate chemical formula
    (3) calculate symmetry number
    (4) calculate thermal properties
    Modifies:
    */
    //## operation make(Graph)
    public static ChemGraph make(Graph p_graph) throws InvalidChemGraphException, ForbiddenStructureException {
        //#[ operation make(Graph)
        try {
        	ChemGraph cg = new ChemGraph(p_graph);
        	cg.addMissingHydrogen();

        	if (cg.repOk()){
        		cg.generateChemicalFormula();
        		cg.calculateSymmetryNumber();
        		cg.generateThermoData();
        		cg.calculateInternalRotor();
        		return cg;
        	}
        	else {
        		throw new InvalidChemGraphException();
        	}
        }
        catch (ForbiddenStructureException e) {
        	throw new ForbiddenStructureException(e.getMessage());
        }



        //#]
    }

    /**
    Requires:
    Effects: read in forbidden structure for ChemGraph
    Modifies: this.forbiddenStructure
    */
    //## operation readForbiddenStructure()
    public static void readForbiddenStructure() throws IOException {
        //#[ operation readForbiddenStructure()
        try {
        	String forbiddenStructureFile = System.getProperty("jing.chem.ChemGraph.forbiddenStructureFile");
        	if (forbiddenStructureFile == null) {
        		System.out.println("undefined system property: jing.chem.ChemGraph.forbiddenStructureFile!");
        		System.out.println("No forbidden structure defined!");
        		return;
        	}

        	FileReader in = new FileReader(forbiddenStructureFile);
        	BufferedReader data = new BufferedReader(in);

        	// step 1: read in structure
        	String line = ChemParser.readMeaningfulLine(data);
        	read: while (line != null) {
        		StringTokenizer token = new StringTokenizer(line);
        		String fgname = token.nextToken();
        		Graph fgGraph = null;
        		try {
        			fgGraph = ChemParser.readFGGraph(data);
        		}
        		catch (InvalidGraphFormatException e) {
        			throw new InvalidFunctionalGroupException(fgname + ": " + e.getMessage());
        		}
        		if (fgGraph == null) throw new InvalidFunctionalGroupException(fgname);
        		FunctionalGroup fg = FunctionalGroup.make(fgname, fgGraph);

        		forbiddenStructure.add(fg);

        		line = ChemParser.readMeaningfulLine(data);
        	}

            in.close();
        	return;
        }
        catch (Exception e) {
        	throw new IOException(e.getMessage());
        }




        //#]
    }

    /**
    Requires:
    Effects: reset centralNode list in this ChemGraph according to node's CentralID information
    Modifies: this.graph.centralNode
    */
    //## operation refreshCentralNode()
    public void refreshCentralNode() {
        //#[ operation refreshCentralNode()
        getGraph().refreshCentralNode();
        //#]
    }

    /**
    Requires:
    Effects: check four aspects:
    (1) if graph.repOk() defined in Graph class
    (2) if graph is connected
    (3) if the graph contents are atoms and bonds
    (4) if the valency are okay for all atom
    (5) if the radical number is in the limit
    */
    //## operation repOk()
    public boolean repOk() {
        //#[ operation repOk()
        // check if the graph is connected
        if (!getGraph().repOk()) return false;

        // a chemical species should be a connected graph
        if (!getGraph().isConnected()) return false;

        // check if the elements stored in graph are atomd/bonds
        if (!graphContentsOk(getGraph())) return false;

        // check if the valency of every atom is satuated
        //if (!valencyOk()) return false;

        // check if the radical number greater than MAX_RADICAL_NUM
        if (getRadicalNumber() > MAX_RADICAL_NUM) return false;

        // check if the oxygen atom number is too large
        if (getOxygenNumber() > MAX_OXYGEN_NUM) return false;
		if (getCarbonNumber() > MAX_CARBON_NUM) return false;

        return true;



        //#]
    }

    /**
    Requires:
    Effects: reset reacting site as the pass-in p_site.
    Modifies: this.graph.centralNode
    */
    //## operation resetReactedSite(HashMap)
    public void resetReactedSite(HashMap p_site) throws SiteNotInSpeciesException {
        //#[ operation resetReactedSite(HashMap)
        setCentralNode(p_site);
        //#]
    }

    //## operation resetStack(Stack)
//svp
  public void resetStack(Stack p_stack) {
         //#[ operation resetStack(Stack)
         while (!p_stack.empty()) {
                 GraphComponent gc = (GraphComponent)p_stack.pop();
                 gc.setMatchedGC(null);
         }
         return;
         //#]
     }


    /**
    Requires:
    Effects: reset the only center to the p_node in this chem graph for thermo calculation
    Modifies: this.graph.centralNode, and  centralIDs in its associated nodes.
    */
    //## operation resetThermoSite(Node)
    public void resetThermoSite(Node p_node) {
        //#[ operation resetThermoSite(Node)
        getGraph().clearCentralNode();
        getGraph().setCentralNode(1,p_node);
        return;
        //#]
    }

    /**
    Requires:
    Effects: reset centreNode list as the pass-in p_site.
    Modifies: this.graph.centralNode
    */
    //## operation setCentralNode(HashMap)
    protected void setCentralNode(HashMap p_site) {
        //#[ operation setCentralNode(HashMap)
        try {
        	Graph g = getGraph();
        	g.clearCentralNode();
        	g.setCentralNodes(p_site);
        }
        catch (NotInGraphException e) {
        	throw new SiteNotInSpeciesException();
        }
        //#]
    }

    /**
    Requires:
    Effects: set GTPP as the thermoGAPP of this chem graph
    Modifies: thermoGAPP
    */
    //## operation setDefaultThermoGAPP()
    public void setDefaultThermoGAPP() {
        //#[ operation setDefaultThermoGAPP()
        thermoGAPP = GATP.getINSTANCE();
        return;
        //#]
    }

    /**
    Requires:
    Effects: return a string of this chemgraph.  the string includes two parts:
    (1) chemical formula
    (2) the string for the graph
    Modifies:
    */
    //## operation toString()
    public String toString() {
        //#[ operation toString()
        String s = "ChemFormula: " + getChemicalFormula() + '\n';
        s = s + getGraph().toStringWithoutCentralID();
        return s;
        //#]
    }

	 public String toString(int i) {
	        //#[ operation toString()
	        String s ="";// "ChemFormula: " + getChemicalFormula() + '\n';
	        s = s + getGraph().toStringWithoutCentralID();
	        return s;
	        //#]
	    }

    /**
    Requires:
    Effects: return a short string for this ChemGraph.  it includes two parts:
    (1) chemical formula
    (2) the graph string without H
    Modifies:
    */
    //## operation toStringWithoutH()
    public String toStringWithoutH() {
        //#[ operation toStringWithoutH()
        String s = "ChemFormula: " + getChemicalFormula() + '\n';
        s = s + getGraph().toStringWithoutCentralIDAndH();
        return s;
        //#]
    }

	/**
    Requires:
    Effects: return a short string for this ChemGraph.  it includes two parts:
    (1) the graph string without H
    Modifies:
    */
    //## operation toStringWithoutH()
    public String toStringWithoutH(int i) {
        //#[ operation toStringWithoutH()
        String s = "";//= "ChemFormula: " + getChemicalFormula() + '\n';
        s = s + getGraph().toStringWithoutCentralIDAndH();
        return s;
        //#]
    }
    /**
    Requires:
    Effects: check all the atom to see if it satisfies:
    Val of atom = radical number + ion number + sum of bond orders.
    if it is satisfied, return true, otherwise return false.
    Modifies:
    */
    //## operation valencyOk()
    public boolean valencyOk() throws InvalidNodeElementException {
        //#[ operation valencyOk()
        Iterator node_iter = graph.getNodeList();
        while (node_iter.hasNext()) {
        	Node node = (Node)node_iter.next();
        	Atom atom = (Atom)node.getElement();
        	double val = atom.getValency();
        	Iterator arc_iter = node.getNeighbor();
        	while (arc_iter.hasNext()) {
        		Bond bond = (Bond)((Arc)(arc_iter.next())).getElement();
        		val -= bond.getOrder();
        	}
        	if (Math.abs(val) >= 1e-10) return false;
        }
        return true;


        //#]
    }

    public int getMAX_OXYGEN_NUM() {
        return MAX_OXYGEN_NUM;
    }

    public static int getMAX_RADICAL_NUM() {
        return MAX_RADICAL_NUM;
    }

    public static HashSet getForbiddenStructure() {
        return forbiddenStructure;
    }

    public int getInternalRotor() {
        return internalRotor;
    }

    public HashSet getSymmetryAxis() {
        return symmetryAxis;
    }

    protected String getUniqueString() {
        return uniqueString;
    }

    public Graph getGraph() {
        return graph;
    }

    public Species getSpecies() {
        return species;
    }

    public void setSpecies(Species p_Species) {
        species = p_Species;
    }

    public GeneralGAPP getThermoGAPP() {
        return thermoGAPP;
    }

    public void setThermoGAPP(GeneralGAPP p_GeneralGAPP) {
        thermoGAPP = p_GeneralGAPP;
    }

}
/*********************************************************************
	File Path	: RMG\RMG\jing\chem\ChemGraph.java
*********************************************************************/


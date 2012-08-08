/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.liv.proteoidviewer.util;

/**
 *
 * @author fghali
 */
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/*This code is an adaption to the IonPredictorNew Class, the only difference being:
 * in the report file, instead of creating 2 columns that contain pep_ionType and m/z,
 * four columns are reported to put pep,ion (b,y..),charge(1+,2+),m/z in separated column.
 *
 *
 *
 * 12-March-2010: code has been modified to only report 1+ and 2+ ions.
 *
 * 22-March-2010: "Carbamidomethyl" (Iodoacetamide derivative) for Cysteine added as permanent mod** delta mono mass: 57.021464
 * 	 "methionine oxidation (MOX)" added as command line argument for variable mod.
 * Ideas: add an array/ harsh map to link variable modification with changes to mass (e.g. MOX = -17.035; NON = 0) then add according value
 * to the Ion methods (i.e. aIons). Three columns needed to match between 1. amino acid; 2. type of modification; 3. delta mass.
 *
 * 24-March-2010: temporarily disabled other ions and only reports b and y ions
 *
 * 10-May-2010: quick modification for fixed PTM (letters can be used to change hashmap are bjoxz): pY = X, pS = J, pT = O.
 *
 * 03-Jun-2010: code modified to produce parent ion mass for MGFgen.java
 *
 * */
public class TheoreticalFragmentation {

    public String resultNM;
    //static String resultNM = ""; 03-06-2010
    public static HashMap<Character, Double> hm;
    private static double Hmass = 1.008;
    private static double carbonmonoxideMass = 28.010;
    private static double ammoniaMass = 17.031;
    private static double waterMass = 18.015;
    private static double HydroniumMass = 19.023;   //hydronium ion = H3O+

     public  TheoreticalFragmentation(String inputPep ) {
        try {
            Map(inputPep);
        } catch (IOException ex) {
            Logger.getLogger(TheoreticalFragmentation.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public void predictIons(String inputPep, String fName) throws IOException {

        //Map(inputPep, fName);
        Map(inputPep);
        aIons(inputPep, fName);
        bIons(inputPep, fName);
        cIons(inputPep, fName);
        xiIons(inputPep, fName);
        yIons(inputPep, fName);
        zIons(inputPep, fName);
    }

    public void MapFill() {  //creates instance of new HashMap hm which catalogues amino acid one letter code as key and mass as value
        hm = new HashMap<Character, Double>();    //values are mass of peptide residues NOT the neutral molecule
        hm.put('A', 71.037114);
        hm.put('R', 156.101111);
        hm.put('N', 114.042927);
        hm.put('D', 115.026943);
        //hm.put ('C', 103.009);
        hm.put('C', 160.030649); // value for fixed modification of Cysteine carbamidomethyl.
        hm.put('E', 129.042593);
        hm.put('Q', 128.058578);
        hm.put('G', 57.021464);
        hm.put('H', 137.058912);
        hm.put('I', 113.084064);
        hm.put('L', 113.084064);
        hm.put('K', 128.094963);
        hm.put('M', 131.040485);
        //hm.put ('MOX', 147.035);
        hm.put('F', 147.068414);
        hm.put('P', 97.052764);
        hm.put('S', 87.032028);
        hm.put('T', 101.047679);
        hm.put('U', 150.95363);
        hm.put('W', 186.079313);
        hm.put('Y', 163.063329);
        hm.put('V', 99.068414);
        //the next three hashes are artificially created to represent different phospho sites (delta mass: 79.966331).
        hm.put('X', 243.029);	//pY
        hm.put('J', 166.998);	//pS
        hm.put('O', 181.014);	//pT
    }
    //public String Map(String inputPep, String fName)throws IOException{

    public String Map(String inputPep) throws IOException {
        MapFill();  //calls MapFill to creat hm instance to work on

        double total = 0;
        char[] pepArray = inputPep.toCharArray(); //converts input string (peptide in amino acid one letter code) to a character array

        for (int i = 0; i < pepArray.length; i++) {   //finds all chars in array and sums mass values
            char ch = pepArray[i];
            if (hm.containsKey(ch)) {
                total += hm.get(ch);
            }
        }
        if (total > 0) {   //returns the neutral mass and m/z values for each peptide
            //resultNM = (inputPep + "\t" + (total + waterMass) + "(neutral)" + "\t" + (total + Hmass) + "(m/z value)" + "\t" + "\r\n"); 03-06-2010
            resultNM = Double.toString(((total + waterMass) + 2 * Hmass) / 2);
        }

        return resultNM;
    }

    public void aIons(String input, String fName) {  //this is to work out masses for a Ions
        try {
            double a = 0;
            double aMass;
            double a2Mass;    //short for two plus ions
            double a3Mass;    //short for three plus ions
            int ionNo = 1;

            char[] pepArray = input.toCharArray();    //converts input string to char array

            String localAion = "";
            String resultIonsA = "";
            String resultIonsA2 = "";
            String resultIonsA3 = "";
            if (hm.containsKey((pepArray[ionNo]))) {
                for (int x = 0; x <= ((pepArray.length) - 2); x++) {
                    a += (hm.get((pepArray[x]))); //sums the masses of peptide residues
                    localAion += pepArray[x];     //finds sequence of a ion in amino acid one letter code
                    aMass = a + Hmass - carbonmonoxideMass;   //calculates protonated (+1) mass of peptide a ions
                    a2Mass = (aMass + Hmass) / 2;
                    a3Mass = (aMass + 2 * Hmass) / 3;

                    resultIonsA += (input + "," + "a" + ionNo + "," + "(1+)" + "," + (aMass) + "\r\n");
                    resultIonsA2 += (input + "," + "a" + ionNo + "," + "(2+)" + "," + (a2Mass) + "\r\n");
                    resultIonsA3 += (input + "," + "a" + ionNo + "," + "(3+)" + "," + (a3Mass) + "\r\n");
                    ionNo++;
                }
                // writeOut(resultIonsA, fName);
                // writeOut(resultIonsA2, fName);
                //writeOut(resultIonsA3, fName);
            }
        } catch (NullPointerException a) {
        }

    }

    public ArrayList<Double> getBIons(String input, String charge) {
         ArrayList<Double> values = new ArrayList();
         double b = 0;
         double bMass;
         double b2Mass;
         char[] pepArray = input.toCharArray();
         if (hm.containsKey((pepArray[1]))) {
                for (int x = 0; x <= ((pepArray.length) - 2); x++) {
                    b += (hm.get((pepArray[x])));

                    bMass = b + Hmass;
                    b2Mass = (bMass + Hmass) / 2;
                    if(charge.equals("1"))
                        values.add(Double.valueOf(bMass));
                    if(charge.equals("2"))
                        values.add(Double.valueOf(b2Mass));


                }
         }

        return values;
    }


    public ArrayList<Double> getYIons(String input, String charge) {
       
         ArrayList<Double> values = new ArrayList();
            double y = 0;
            double yMass;
            double y2Mass;

         char[] pepArray = input.toCharArray();
         if (hm.containsKey((pepArray[1]))) {
                for (int x = 0; x <= ((pepArray.length) - 2); x++) {
                    y += (hm.get((pepArray[x])));  //sums the masses of peptide residues
                    yMass = y + HydroniumMass;  //calculates protonated (+1) mass of peptide y ions (H3O+ as hydrolysis fragmentation reaction)
                    y2Mass = (yMass + Hmass) / 2;
                    if(charge.equals("1"))
                        values.add(Double.valueOf(yMass));
                    if(charge.equals("2"))
                        values.add(Double.valueOf(y2Mass));


                }
         }

        return values;
    }

    public void bIons(String input, String fName) {  //this is to work out masses for b Ions
        try {
            double b = 0;
            double bMass;
            double b2Mass;
            double b3Mass;
            int ionNo = 1;

            char[] pepArray = input.toCharArray();    //converts input string to char array

            String localBion = "";
            String resultIonsB = "";
            String resultIonsB2 = "";
            String resultIonsB3 = "";
            if (hm.containsKey((pepArray[ionNo]))) {
                for (int x = 0; x <= ((pepArray.length) - 2); x++) {
                    b += (hm.get((pepArray[x]))); //sums the masses of peptide residues
                    localBion += pepArray[x];     //finds sequence of b ion in amino acid one letter code
                    bMass = b + Hmass;      //calculates protonated (+1) mass of peptide b ions
                    b2Mass = (bMass + Hmass) / 2;
                    b3Mass = (bMass + 2 * Hmass) / 3;

                    resultIonsB += (input + "," + "b" + ionNo + "," + "(1+)" + "," + (bMass) + "\r\n");
                    resultIonsB2 += (input + "," + "b" + ionNo + "," + "(2+)" + "," + (b2Mass) + "\r\n");
                    resultIonsB3 += (input + "," + "b" + ionNo + "," + "(3+)" + "," + (b3Mass) + "\r\n");
                    ionNo++;
                }
                writeOut(resultIonsB, fName);
                writeOut(resultIonsB2, fName);
                //writeOut(resultIonsB3, fName);
            }
        } catch (NullPointerException a) {
        }

    }

    public void cIons(String input, String fName) {  //this is to work out the masses for c Ions
        try {
            double c = 0;
            double cMass;
            double c2Mass;
            double c3Mass;
            int ionNo = 1;

            char[] pepArray = input.toCharArray();    //converts input string to char array

            String localCion = "";
            String resultIonsC = "";
            String resultIonsC2 = "";
            String resultIonsC3 = "";
            if (hm.containsKey((pepArray[ionNo]))) {
                for (int x = 0; x <= ((pepArray.length) - 2); x++) {
                    c += (hm.get((pepArray[x]))); //sums the masses of peptide residues
                    localCion += pepArray[x];     //finds sequence of c ion in amino acid one letter code
                    cMass = c + Hmass + ammoniaMass;    //calculates protonated (+1) mass of peptide c ions
                    c2Mass = (cMass + Hmass) / 2;
                    c3Mass = (cMass + 2 * Hmass) / 3;

                    resultIonsC += (input + "," + "c" + ionNo + "," + "(1+)" + "," + (cMass) + "\r\n");
                    resultIonsC2 += (input + "," + "c" + ionNo + "," + "(2+)" + "," + (c2Mass) + "\r\n");
                    resultIonsC3 += (input + "," + "c" + ionNo + "," + "(3+)" + "," + (c3Mass) + "\r\n");
                    ionNo++;
                }
                //  writeOut(resultIonsC, fName);
                // writeOut(resultIonsC2, fName);
                //writeOut(resultIonsC3, fName);
            }
        } catch (NullPointerException a) {
        }

    }

    public void xiIons(String input, String fName) { //this is to work out the masses for x Ions

        try {
            double xi = 0;
            double xiMass;    //for x ions
            double x2Mass;
            double x3Mass;
            int ionNo = 1;

            StringBuffer reverse = new StringBuffer();
            for (int x = (input.length() - 1); x >= 0; x--) {
                reverse.append(input.charAt(x));
            }
            String xiString = reverse.toString();

            char[] pepArray = xiString.toCharArray();  //converts input string (peptide in amino acid one letter code) to char array

            String localXiion = "";
            String resultIonsXi = "";
            String resultIonsX2 = "";
            String resultIonsX3 = "";
            if (hm.containsKey((pepArray[ionNo]))) {
                for (int x = 0; x <= ((pepArray.length) - 2); x++) {
                    xi += (hm.get((pepArray[x])));  //sums the masses of peptide residues
                    localXiion += pepArray[x];  //finds sequence of y ion in amino acid one letter code
                    xiMass = xi + waterMass + carbonmonoxideMass - Hmass;  //calculates protonated (+1) mass of peptide y ions (H3O+ as hydrolysis fragmentation reaction)
                    x2Mass = (xiMass + Hmass) / 2;
                    x3Mass = (xiMass + 2 * Hmass) / 3;

                    resultIonsXi += (input + "," + "x" + ionNo + "," + "(1+)" + "," + (xiMass) + "\r\n");
                    resultIonsX2 += (input + "," + "x" + ionNo + "," + "(2+)" + "," + (x2Mass) + "\r\n");
                    resultIonsX3 += (input + "," + "x" + ionNo + "," + "(3+)" + "," + (x3Mass) + "\r\n");
                    ionNo++;
                }
                //  writeOut(resultIonsXi, fName);
                //  writeOut(resultIonsX2, fName);
                //writeOut(resultIonsX3, fName);
            }
        } catch (NullPointerException a) {
        }

    }

    public void yIons(String input, String fName) {  //this is to work out the masses for y Ions

        try {
            double y = 0;
            double yMass;
            double y2Mass;
            double y3Mass;
            int ionNo = 1;

            StringBuffer reverse = new StringBuffer();
            for (int x = (input.length() - 1); x >= 0; x--) {
                reverse.append(input.charAt(x));
            }
            String yString = reverse.toString();

            char[] pepArray = yString.toCharArray();  //converts input string (peptide in amino acid one letter code) to char array

            String localYion = "";
            String resultIonsY = "";
            String resultIonsY2 = "";
            String resultIonsY3 = "";
            if (hm.containsKey((pepArray[ionNo]))) {
                for (int x = 0; x <= ((pepArray.length) - 2); x++) {
                    y += (hm.get((pepArray[x])));  //sums the masses of peptide residues
                    localYion += pepArray[x];  //finds sequence of y ion in amino acid one letter code
                    yMass = y + HydroniumMass;  //calculates protonated (+1) mass of peptide y ions (H3O+ as hydrolysis fragmentation reaction)
                    y2Mass = (yMass + Hmass) / 2;
                    y3Mass = (yMass + 2 * Hmass) / 3;

                    resultIonsY += (input + "," + "y" + ionNo + "," + "(1+)" + "," + (yMass) + "\r\n");
                    resultIonsY2 += (input + "," + "y" + ionNo + "," + "(2+)" + "," + (y2Mass) + "\r\n");
                    resultIonsY3 += (input + "," + "y" + ionNo + "," + "(3+)" + "," + (y3Mass) + "\r\n");
                    ionNo++;
                }
                writeOut(resultIonsY, fName);
                writeOut(resultIonsY2, fName);
                //writeOut(resultIonsY3, fName);
            }
        } catch (NullPointerException a) {
        }

    }

    public void zIons(String input, String fName) {  //this is to work out the masses for z Ions

        try {
            double z = 0;
            double zMass;
            double z2Mass;
            double z3Mass;
            int ionNo = 1;

            StringBuffer reverse = new StringBuffer();
            for (int x = (input.length() - 1); x >= 0; x--) {
                reverse.append(input.charAt(x));
            }
            String zString = reverse.toString();

            char[] pepArray = zString.toCharArray();  //converts input string (peptide in amino acid one letter code) to char array

            String localZion = "";
            String resultIonsZ = "";
            String resultIonsZ2 = "";
            String resultIonsZ3 = "";
            if (hm.containsKey((pepArray[ionNo]))) {
                for (int x = 0; x <= ((pepArray.length) - 2); x++) {
                    z += (hm.get((pepArray[x])));  //sums the masses of peptide residues
                    localZion += pepArray[x];
                    zMass = z + HydroniumMass - ammoniaMass;
                    z2Mass = (zMass + Hmass) / 2;
                    z3Mass = (zMass + 2 * Hmass) / 3;

                    resultIonsZ += (input + "," + "z" + ionNo + "," + "(1+)" + "," + (zMass) + "\r\n");
                    resultIonsZ2 += (input + "," + "z" + ionNo + "," + "(2+)" + "," + (z2Mass) + "\r\n");
                    resultIonsZ3 += (input + "," + "z" + ionNo + "," + "(3+)" + "," + (z3Mass) + "\r\n");
                    ionNo++;
                }
                // writeOut(resultIonsZ, fName);
                // writeOut(resultIonsZ2, fName);
                //writeOut(resultIonsZ3, fName);
            }
        } catch (NullPointerException a) {
        }


    }

    public void writeOut(String res, String fName) { //writes out the result string 'res' to the file 'name' in the specified folder
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(fName, true));
            out.append(res);
            out.close();
        } catch (IOException e) {
            System.out.println("Exception" + e);
        }
    }

    public static void main(String[] args) throws Exception {

        int fLength = args.length;

        if (fLength >= 1) {

            String pepIn = args[0];
            String fName = args[1];
            //String variMod = args[2];
            String[] pepSeq = pepIn.split(","); //separate input peptide sequences

            for (int i = 0; i < pepSeq.length; i++) {

                String inputPep = pepSeq[i];

           //     TheoreticalFragmentation peptide1 = new TheoreticalFragmentation();

          //      peptide1.predictIons(inputPep, fName); //call predictIons method

            }

        } else {
            System.out.println("Please enter" + "/r/n" + "1.	The amino acid sequences of your peptide using upper case one-letter code "
                    + "(eg ARND, separate with comma)" + "/r/n" + "2.  names and locations for output file (eg c:/Result/Sample1.csv)"
                    + "/r/n" + "3. variable modifications (eg MOX for methionine oxidation; NON for no vari mod needed)"
                    + "/r/n" + "as command-line arguments");
        }
    }
}

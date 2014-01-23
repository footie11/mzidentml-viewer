/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package uk.ac.liv.proteoidviewer.util;

import java.util.Comparator;

/**
 *
 * @author fghali
 */
public class DoubleComparator implements Comparator {

    @Override
    public int compare(Object o1, Object o2) {
        double P = Double.parseDouble(o1.toString());
        double Q = Double.parseDouble(o2.toString());

        if (P > Q) {
            return 1;
        } else if (P < Q) {
            return -1;
        } else {
            return 0;
        }

    }

    
    

    
}

/*****************************************************************************
*
* Authors: SciCom research group-E.S.I. Universidad de Castilla-La Mancha
*          Paseo de la Universidad 4, 13004-Ciudad Real. SPAIN
*
* Release date: July 7, 2016
*
* Purpose: To implement class SparseArray
*
*****************************************************************************/

package fuzzy_communities;

import java.io.*;
import java.util.*;

/**
 *
 * @author joseangel
 */
public class Fuzzy_comunities {

    public static void main(String[] args) {
        String file = "Graph1.dat";
        int N = 7;
        int c = 2;
        
        SparseArray A = new SparseArray(N);
        double [][] U = new double [c][N];

        System.out.println("Creating graph " + file + "...");
        try {
            File fd = new File("src/" + file);
            Scanner read = new Scanner(fd);
            
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    A.put(i, j, read.nextInt());
                }
            }
            read.close();
        } catch (FileNotFoundException e) {
            System.err.println("File " + file + " could not be opened.");
        }

        FuzzyClustering fc = new FuzzyClustering(A, c, N);
        U = fc.findCommunities_GC();
        
        // Print U
        int i, j;
        for(i = 0; i < c; i++){
            System.out.print("Class " + i + ": ");
            for(j = 0; j < N; j++){
                System.out.print(U[i][j] + "  ");
            }
            System.out.print("\n");
        }
        
        
    }

}

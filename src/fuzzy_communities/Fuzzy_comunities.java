/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
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
        int[][] A = new int[N][N];
        double [][] U = new double [c][N];

        System.out.println("Creating graph " + file + "...");
        try {
            File fd = new File("src/" + file);
            Scanner read = new Scanner(fd);
            
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    A[i][j] = read.nextInt();
                }
            }
            read.close();
        } catch (FileNotFoundException e) {
            System.err.println("File " + file + " could not be opened.");
        }

        FuzzyClustering fc = new FuzzyClustering(A, c);
        U = fc.getCommunities();
        
        //Print U
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

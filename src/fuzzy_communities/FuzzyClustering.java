/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fuzzy_communities;

import java.util.Random;

/**
 * 
 * @author joseangel
 */
public class FuzzyClustering {
    private int [][] A; // Adjacency matrix
    private int N; // Number of vertices
    private int c; // Number of communities
    
    /**
     * Constructor method of the class. Obtain the main data we are going to use
     * to make the rest of computations.
     * 
     * @param A The adjacency matrix we are working with.
     * @param c The number of communities we are interesting to calculate.
     */
    public FuzzyClustering(int[][] A, int c) {
        this.A = A;
	this.c = c;
        N = A.length;
    }
    
    /**
     * Get communities.
     * Applies the algorithm to calculate the fuzzy communities of each vertex
     * of the graph.
     * @return  An array U[c][N] with the pertenence of each vertex to each 
     *          class.
     */
    public double [][] getCommunities(){
        double [][] U, dD, s;
        boolean continueIteration = true;
        int t, n_max_iterations, i, j, k, l;
        double max, inverse_of_c, epsilon, sumatory, aux;
        s = new double [N][N]; 
        dD = new double [c][N]; // Partial derivative of D
        for(i = 0; i < c; i++){
            for(j = 0; j < N; j++){
                dD[i][j] = 0;
            }
        }
        epsilon = 0.001;
        n_max_iterations = 1000;
        
        // Step 1
        U = getU();
        t = 0;
        
        while(continueIteration){
            // Calculate s
            for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                    s[i][j] = 0;
                    for (k = 0; k < c; k++) {
                        s[i][j] += U[k][i] * U[k][j];
                    }
                }
            }
            System.out.println(" The value of D is " + D(s));
            
            
            // Step 2
            max = Double.NEGATIVE_INFINITY;
            inverse_of_c = 1.0/c;
            for(k = 0; k < c; k++){
                for(l = 0; l < N; l++){
                    sumatory = 0.0;
                    for(i = 0; i < N; i++){
                        sumatory += ((A[i][l] - s[i][l]) + (A[l][i] - s[l][i])) 
                                * (inverse_of_c - U[k][i]);
                    }
                    dD[k][l] = 2 * sumatory; 
                    if(Math.abs(dD[k][l]) > max){
                        max = Math.abs(dD[k][l]);
                    }
                }
            }
            System.out.println("Iteration " + t + ". Max abs val of dD " + max
                    + "\n");

            // Step 3
            if(max < epsilon || t > n_max_iterations){
                continueIteration = false;
            }else{
                // Step 4
                for(i = 0; i < c; i++){
                    for(j = 0; j < N; j++){
                        U[i][j] =  U[i][j] - alpha() * dD[i][j];
                    }
                }
                
                // Step 5
                t++;
                
            }
        }
        
        return U;
    }
    
    /**
     * Calculates the initial partition array U used in the algorithm.
     * It uses a Gamma distribution to generate random numbers that are going to
     * be used to fill the array. Then, the array is normalized so that
     * the sum of a node membership to all classes is 1.
     * @return U
     */
    private double [][] getU(){
        double [][] U = new double [c][N];
        int i, j;
        double sum_colum, d_constant, c_constant, alpha, Z_ran_variable, 
                U_ran_variable, V, aux;
        boolean continueIteration;
        Random ran = new Random();
        alpha = 1.0;
        
        // Gamma(alpha,1) Generator for alpha >= 1
        // Step 1 - Constant values
        d_constant = alpha - 1.0/3.0;
        c_constant = 1.0/Math.sqrt(9 * d_constant);
        
        for(j = 0; j < N; j++){
            sum_colum = 0;
            for (i = 0; i < c; i++) {
                continueIteration = true;
                while (continueIteration) {
                    // Step 2 - Generate random variables Z and U
                    Z_ran_variable = ran.nextGaussian(); // Random value from a 
                                                         // Normal distribution
                    U_ran_variable = ran.nextDouble(); // Random value from a 
                                                       //Uniform distribution

                    // Step 3 - Finish condition
                    aux = (1 + c_constant * Z_ran_variable);
                    V = aux * aux * aux;
                    if (Z_ran_variable > -1.0 / c_constant
                            && Math.log(U_ran_variable) < (1.0 / 2 * 
                            Z_ran_variable * Z_ran_variable + d_constant - 
                            d_constant * V + d_constant * Math.log(V))) {
                        U[i][j] = d_constant * V;
                        continueIteration = false;
                    }
                } // End while
                sum_colum += U[i][j];
            }
            
            for(i = 0; i < c; i++){
                U[i][j] = U[i][j] / sum_colum;
            }
        }
        return U;
    }

    /**
     * Calculates D.
     * @param s     The value of s in a given moment on the graph. 
     * @return D
     */
    private double D(double[][] s) {
        double aux, D = 0.0;
        int i, j;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                aux = A[i][j] - s[i][j];
                D += aux * aux;
            }
        }
        return D;
    }
    
    /**
     * Calculates alpha constant.
     * Alpha is a small step size constant used to calculate the next partition
     * in the iteration.
     * Uses line minimization bracketing the minimum by the golden rule and 
     * minimum location using inverse parabolic interpolation or, if it does 
     * not work, golden ratio search.
     * @return alpha
     */
    private double alpha(){
        double a, b, c, x, fa, fb, fc, fx, aux, tol, golden =  1.6180339887;
        tol = 3e-8; //Square root of machine double precision ??????????????????????????????????
        
        /* Bracketing the minimum using golden ratio */
        // Initializing
	a = 0.0;
	b = 1.0;
	fa = function(a);
	fb = function(b);
	
	if (fb > fa) {  // We always want fb < fa (we always 
		aux = a;      // search in the decreasing direction)
		a = b;
		b = aux;
		aux = fa;
		fa = fb;
		fb = aux;
	}
        
        c = b + golden * (b-a); // Increasing by the golden ratio of the
	                        // a-b segment
	fc = function(c) ;
	
	while (fc < fb) { // Searching until the function increases
		a = b;
		fa = fb;
		b = c;
		fb = fc;
		c = b + golden * (b-a);
		fc = function(c);
	}
        
        /* Minimum bracketed */
        
        // Computing minimum applying Brent's approach
        x = brent(a, b, c, tol);
        
        return x;
    }
    
    /**
     * This method performs a 1-dimensional minimization.
     * It implements Brent's method which combines a golden-section search and 
     * parabolic interpolation.
     * @param a     Left endpoint of initial interval
     * @param b     Point between a and c and f(b) is less than f(a) and f(c)
     * @param c     Right endpoint of initial interval
     * @param tol   Desired length of the interval in which the minimum will be 
     *              determined to lie
     * @return x    The abscissa of the minimum value of the function in the
     *              interval [a,c]
     */
    private double brent(double a, double b, double c, double tol){
        double CGOLD, d, e, eps, xm, p, q, r, tol1, t2, u, v, w, fu, fv , fw, 
                fx, x, tol3;
        int t, tmax = 100;
        CGOLD = .5 * (3.0 - Math.sqrt(5.0));
        d = 0.0;

        eps = 1.2e-16; //Machine precision
        tol1 = eps + 1.0;
        eps = Math.sqrt(eps);

        x = w = v = b;
        e = 0.0;
        fw = fv = fx = function(x);
        tol3 = tol / 3.0;

        xm = .5 * (a + b);
        tol1 = eps * Math.abs(x) + tol3;
        t2 = 2.0 * tol1;
       
        //main loop
        for (t = 0; Math.abs(x - xm) > (t2 - .5 * (b - a)) && t < tmax; t++) {

            // REST OF CODE
            
            
            xm = .5 * (a + b);
            tol1 = eps * Math.abs(x) + tol3;
            t2 = 2.0 * tol1;
        } // End main loop
        
        if(t == tmax){
            throw new RuntimeException("ERROR. Max number of iterations allowed"
                    + " in brent method.");
        }
        
        return x;
    }
    
    
}

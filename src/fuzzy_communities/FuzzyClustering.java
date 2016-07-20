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

import java.util.Random;

/**
 * 
 * @author joseangel
 */
public class FuzzyClustering {
    private SparseArray A; // Adjacency matrix
    private int N; // Number of vertices
    private int c; // Number of communities
    private int two_m; // Two times the number of edges in the network.
    
    /**
     * Constructor method of the class. Obtain the main data we are going to use
     * to make the rest of computations.
     * 
     * @param A The adjacency matrix we are working with.
     * @param c The number of communities we are interesting to calculate.
     * @param N The number of nodes of the Graph.
     * @param two_m Two times the number of edges in the network.
     */
    public FuzzyClustering(SparseArray A, int c, int N, int two_m) {
        this.A = A;
	this.c = c;
        this.N = N;
        this.two_m = two_m;
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
        double max, inverse_of_c, epsilon, sumatory, aux, alpha;
        s = new double [N][N]; 
        dD = new double [c][N]; // Partial derivative of D
        for(i = 0; i < c; i++){
            for(j = 0; j < N; j++){
                dD[i][j] = 0;
            }
        }
        epsilon = 0.001;
        n_max_iterations = 100;
        
        // Step 1
        U = getU();
        t = 0;
        
        while(continueIteration){
            // Calculate s
            for(i = 0; i < N; i++){
                for(j = 0; j < N; j++){
                    s[i][j] = 0;
                    for(k = 0; k < c; k++){
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
                        sumatory += ((A.get(i, l) - s[i][l]) + (A.get(l, i) - s[l][i])) 
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
                alpha = alpha(U, dD);
                for(i = 0; i < c; i++){
                    for(j = 0; j < N; j++){
                        U[i][j] =  U[i][j] - alpha * dD[i][j];
                    }
                }
                
                // Step 5
                t++;
                
            }
        }
        
        return U;
    }
    
    /**
     * Finding communities.
     * Applies the algorithm to calculate the fuzzy communities of each vertex
     * of the graph.
     * @return  An array U[c][N] with the pertenence of each vertex to each 
     *          class.
     */
    public double [][] findCommunities(){

      double [][] dD;
      double [][] U;
      double alpha, max, inverse_of_c, epsilon, aux, sumatory, zero, two;
      boolean continueIteration = true;
      int t, n_max_iterations, i, j, k, l, m;
      
      dD = new double [c][N]; // Partial derivative of D
      
      
      // General constants
      zero = 0.0;
      two  = 2.0;
      inverse_of_c = 1.0/c;
      epsilon = 0.001;         // Derivatives limit
      n_max_iterations = 100;  // Max number of iterations allowed
      
      
      // Step 1
      U = getU();
      t = 0;
      
      while(continueIteration){      
        
        // Step 2
        max = Double.NEGATIVE_INFINITY;
        for(l = 0; l < N; l++){
          for(k = 0; k < c; k++){
            sumatory = zero;
            for(i = 0; i < N; i++){
              aux = zero;
              for (m = 0; m < c; m++) {   // Computing s matrix elements 
                aux += U[m][i] * U[m][l]; // on the fly
              }
              sumatory += ((A.get(i, l) - aux) + (A.get(l, i) - aux)) 
                        * (inverse_of_c - U[k][i]);
            }
            dD[k][l] = two * sumatory; 
            if(Math.abs(dD[k][l]) > max) max = Math.abs(dD[k][l]);
          }
        }
        System.out.println("Iteration " + t + ". Max abs val of dD " + max
                           + "\n");

        // Step 3
        if(max < epsilon || t == n_max_iterations){ // End conditions
          continueIteration = false;
        } else {
          // Step 4
          alpha = alpha (U, dD);
          for(i = 0; i < c; i++){
            for(j = 0; j < N; j++){
              U[i][j] -=  alpha * dD[i][j];
            }
          }
          
          // Step 5
          t++;  
        }
      }
              
      if(t == n_max_iterations) throw new RuntimeException("ERROR. Max number "
                                 + "of iterations allowed in findCommunities "
                                 + "method exceeded");
      
      return U;
    }    
    
    /**
     * Finding communities using conjugate gradient.
     * Applies the algorithm to calculate the fuzzy communities of each vertex
     * of the graph.
     * @return  An array U[c][N] with the pertenence of each vertex to each 
     *          class.
     */
    public double [][] findCommunities_GC(){

      double [][] dD, g, h;
      double [][] U;
      double alpha, max, inverse_of_c, epsilon, aux, sumatory, zero, two, 
             gg, dgg, gam;
      boolean continueIteration = true;
      int t, n_max_iterations, i, j, k, l, m;
      
      dD = new double [c][N]; // Partial derivative of D
      g  = new double [c][N]; 
      h  = new double [c][N];     
      
      // General constants
      zero = 0.0;
      two  = 2.0;
      inverse_of_c = 1.0/c;
      epsilon = 0.001;         // Derivatives limit
      n_max_iterations = 100;  // Max number of iterations allowed
      
      
      // Initializing U
      U = getU();
      t = 0;

      // Computing gradients
      for(l = 0; l < N; l++){
        for(k = 0; k < c; k++){
          sumatory = zero;
          for(i = 0; i < N; i++){
            aux = zero;
            for (m = 0; m < c; m++) {   // Computing s matrix elements 
              aux += U[m][i] * U[m][l]; // on the fly
            }
            sumatory += ((A.get(i, l) - aux) + (A.get(l, i) - aux)) 
                      * (inverse_of_c - U[k][i]);
          }
          dD[k][l] = two * sumatory; 
        }
      }

      for (i = 0; i < c; i++){
        for (j =0; j < N; j++) {
          g[i][j] = -dD[i][j];
          dD[i][j] = h[i][j] = g[i][j];
        }
      }
     
      
      while(continueIteration){      
        
        // Computing minimum along the search direction (-dD)
        alpha = alpha (U, dD);
        for(i = 0; i < c; i++){
          for(j = 0; j < N; j++){
            U[i][j] -=  alpha * dD[i][j];
          }
        }
        
        // Computing gradients
        max = Double.NEGATIVE_INFINITY;
        for(l = 0; l < N; l++){
          for(k = 0; k < c; k++){
            sumatory = zero;
            for(i = 0; i < N; i++){
              aux = zero;
              for (m = 0; m < c; m++) {   // Computing s matrix elements 
                aux += U[m][i] * U[m][l]; // on the fly
              }
              sumatory += ((A.get(i, l) - aux) + (A.get(l, i) - aux)) 
                        * (inverse_of_c - U[k][i]);
            }
            dD[k][l] = two * sumatory; 
            if(Math.abs(dD[k][l]) > max) max = Math.abs(dD[k][l]);
          }
        }
        System.out.println("Iteration " + t + ". Max abs val of dD " + max
                           + "\n");

        // Step 3
        if(max < epsilon || t == n_max_iterations){ // End conditions
          continueIteration = false;
        } else {
          // Applying conjugate gradient
          
          gg = dgg = zero;
          for (i = 0; i < c; i++){
            for (j =0; j < N; j++) {
              gg += g[i][j] *g[i][j];
              dgg += (dD[i][j] + g[i][j]) * dD[i][j];
            }
          }
          gam = dgg / gg; 

          for (i = 0; i < c; i++){
            for (j =0; j < N; j++) {
              g[i][j] = -dD[i][j];
              dD[i][j] = h[i][j] = g[i][j] + gam * h[i][j];
            }
          }
          // Step 5
          t++;  
        }
      }
              
      if(t == n_max_iterations) throw new RuntimeException("ERROR. Max number "
                                 + "of iterations allowed in findCommunities_GC"
                                 + " method exceeded");
      
      return U;
    }    
    
    /**
     * Determines the modularity of the graph.
     * @param U     The current U matrix.
     * @return      Q
     */
    private double modularity(double[][] U) {
        int i, j, k, ki;
        double sij, Q = 0.0, one_over_2m = 1.0 / two_m;

        for(i = 0; i < N; i++){
            ki = A.degree(i);
            for(j = 0; j < N; j++){
                sij = 0.0;
                for(k = 0; k < c; k++){
                    sij += U[k][i] * U[k][j];
                }
                Q += (A.get(i, j) - ki * A.degree(j) * one_over_2m) * sij;
            }
        }
        Q = Q * one_over_2m;
        return Q;
    }
        
    /**
     * Calculates the initial partition array U used in the algorithm.
     * It uses a Gamma distribution to generate random numbers that are going to
     * be used to fill the array. Then, the array is normalized so that
     * the sum of a node membership to all classes is 1.
     * @return  U
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
     * @return      The value of D
     */
    private double D(double[][] s) {
        double aux, D = 0.0;
        int i, j;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                aux = A.get(i, j) - s[i][j];
                D += aux * aux;
            }
        }
        return D;
    }
    
    /**
     * Calculates alpha constant.
     * Alpha is a small step size constant used to calculate the next partition
     * in the iteration.
     * Uses a line minimization implementing Brent's method which brackes the 
     * minimum by the golden rule and minimum location using inverse parabolic 
     * interpolation or, if it does not work, golden ratio search.
     * @param U     The current U matrix.
     * @param dD    The current value of the derivate of D.
     * @return      alpha
     */
    private double alpha(double[][] U, double[][] dD){
        double ax, bx, cx, x, fa, fb, fc, fx, aux, tol, golden =  1.6180339887,
                CGOLD, d, e, eps, xm, p, q, r, tol1, t2, u, v, w, fu, fv , fw, 
                tol3, a, b;
        tol = 3e-8; //Square root of machine double precision
        int t, tmax = 100;
        CGOLD = .5 * (3.0 - Math.sqrt(5.0));
        d = 0.0;
        
        /* Bracketing the minimum using golden ratio */
        // Initializing
	ax = 0.0;
	bx = 1.0;
	fa = function(ax, U, dD);
	fb = function(bx, U, dD);
	
	if (fb > fa) {  // We always want fb < fa (we always 
		aux = ax;      // search in the decreasing direction)
		ax = bx;
		bx = aux;
		aux = fa;
		fa = fb;
		fb = aux;
	}
        
        cx = bx + golden * (bx-ax); // Increasing by the golden ratio of the
	                        // a-b segment
	fc = function(c, U, dD) ;
	
	while (fc < fb) { // Searching until the function increases
		ax = bx;
		fa = fb;
		bx = c;
		fb = fc;
		cx = bx + golden * (bx-ax);
		fc = function(c, U, dD);
	}
        
        /* Minimum bracketed */
        
        // Computing minimum applying Brent's approach
        
        a = (ax < cx ? ax : cx);
        b = (ax > cx ? ax : cx);

        eps = 1.2e-16; //Machine precision
        tol1 = eps + 1.0;
        eps = Math.sqrt(eps);

        x = w = v = bx;
        e = 0.0;
        fw = fv = fx = function(x, U, dD);
        tol3 = tol / 3.0;

        xm = .5 * (a + b);
        tol1 = eps * Math.abs(x) + tol3;
        t2 = 2.0 * tol1;
       
        //main loop
        for(t = 0; Math.abs(x - xm) > (t2 - .5 * (b - a)) && t < tmax; t++){
            p = q = r = 0.0;
            
            if(Math.abs(e) > tol1){ // fit the parabola
                r = (x - w) * (fx - fv);
                q = (x - v) * (fx - fw);
                p = (x - v) * q - (x - w) * r;
                q = 2.0 * (q - r);
                if(q > 0.0){
                    p = -p;
                }else{
                    q = -q;
                }
                r = e;
                e = d;
            }
            
            if((Math.abs(p) < Math.abs(.5 * q * r))
                    && (p > q * (a - x))
                    && (p < q * (b - x))){ // a parabolic interpolation step
                d = p / q;
                u = x + d;
                if(((u - a) < t2) || ((b - u) < t2)){ // f must not be evaluated
                                                      // too close to a or b
                    d = tol1;
                    if(x >= xm){
                        d = -d;
                    }
                }
            }else{ // a golden-section step        
                if(x < xm){
                    e = b - x;
                }else{
                    e = a - x;
                }
                d = CGOLD * e;
            }

            if(Math.abs(d) >= tol1){ // f must not be evaluated too close to x
                u = x + d;
            }else if(d > 0.0){ 
                u = x + tol1;
            }else{
                u = x - tol1;
            }
            fu = function(u, U, dD);
                          
            if(fu <= fx){
                if(u >= x){
                    a = x;
                }else{
                    b = x;
                }
                v = w;
                fv = fw;
                w = x;
                fw = fx;
                x = u;
                fx = fu;
            }else{
                if(u < x){
                    a = u;
                }else{
                    b = u;
                }
                if(fu <= fw || w == x){
                    v = w;
                    w = u;
                    fv = fw;
                    fw = fu;
                }else if(fu <= fv || v == x || v == w){
                    v = u;
                    fv = fu;
                }
            }
            
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
    
    /**
     * Function that calculates the error (D) used in the calculous of the 
     * alpha.
     * @param alpha     The value of alpha in a given moment.
     * @param U         The current U matrix.
     * @param dD        The current value of the derivate of D.
     * @return          The value of D for the given input.
     */
    private double function(double alpha, double[][] U, double[][] dD){
        double[][] Uaux = new double[c][N], s = new double[N][N];
        int i, j, k;
        double aux, D = 0.0;
        
        for(i = 0; i < c; i++){
            for(j = 0; j < N; j++){
                Uaux[i][j] = U[i][j] - alpha * dD[i][j];
            }
        }

        for(i = 0; i < N; i++){
            for(j = 0; j < N; j++){
                s[i][j] = 0.0;
                for (k = 0; k < c; k++) {
                    s[i][j] += Uaux[k][i] * Uaux[k][j];
                }
            }
        }
        
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                aux = A.get(i, j) - s[i][j];
                D += aux * aux;
            }
        }
        
        return D;    // Returns the current value of D.
    }
    
}


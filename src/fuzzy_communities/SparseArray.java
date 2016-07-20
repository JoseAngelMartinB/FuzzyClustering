/*****************************************************************************
*
* Authors: SciCom research group-E.S.I. Universidad de Castilla-La Mancha
*          Paseo de la Universidad 4, 13004-Ciudad Real. SPAIN
*
* Release date: July 6, 2016
*
* Purpose: To implement class SparseArray
*
*****************************************************************************/

/**
 * @file
 * Implementation of class SparseArray
*/

package fuzzy_communities;

import java.util.Hashtable;

/**
 * \brief Defines a sparse array
 *
 * This class represents a sparse array where only non-zero real values are 
 * stored. However, all the elements can be referenced and a 0.0 value
 * is returned for elements with zero value.
 * @author SciCom research group-E.S.I. Universidad de Castilla-La Mancha.
 *         Paseo de la Universidad 4, 13004-Ciudad Real. SPAIN
 * 
 * @param <E> The class to which the element stored in the node belongs
 *            to.
 * @date July 6, 2016
 */
public class SparseArray {
  private Hashtable <Integer, Double> [] row; // Array of hashtables
  private int N; // Number of rows in the array
  
  /**
   * Constructor
   * @param N The number of rows in the array
   */
  public SparseArray(int N) {
    this.N = N;
    row = new Hashtable[N];
    for (int i = 0; i < N; i++) {
      row [i] = new Hashtable<>(5);  // Initial size 5 elements
    }
  }
  
  /**
   * Method for adding the i, j element to the sparse array. Only non-zero
   * elements are added
   * @param i The row of the element
   * @param j The colummn of the element
   * @param value The value of the (i, j) element
   */
  public void put (int i, int j, double value) {
    if (value != 0.0) row[i].put(j, value); // In a row, elements are indexed by
                                            // colummn
  }
  
  /**
   * Method for retrieving the i, j element from the sparse array
   * @param i The row of the element
   * @param j The colummn of the element
   * @return The value of the (i, j) element
   */
  public double get (int i, int j) {
    Double value = row[i].get(j); // Retrieving value in row i with key j
    if (value == null) return 0.0;
    return value;
  }
  
  /**
   * Method for retrieving the degree of the i element from the sparse array
   * @param i The row of the element
   * @return The degree of the element.
   */
  public int degree(int i){
    int value = row[i].size();
    return value;
  }
}

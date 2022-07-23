/** ---------------------------------------------------------------------------
Functions for converting between different types of reporting triangles
---------------------------------------------------------------------------- */
 
/**
  * Takes a reporting triangle by date of reference, and reshapes it into a
  * reporting triangle by date of report, excluding left-truncated dates.
  *
  * @param triangle DxT matrix representing a reporting triangle, i.e. columns
  * represent the date of reference and rows represent the delay
  * 
  * @param D maximum delay
  * 
  * @return A Dx(T-D) matrix representing a reporting triangle by date
  * of report. Columns represent the date of report (since first day that is
  * not left-truncated) and rows represent the delay.
  */
  matrix reporting_triangle_by_report(matrix triangle, int D){
    int T_all = cols(triangle);
    matrix[D+1, T_all-D] triangle_rep = rep_matrix(0, D+1, T_all-D);
    for (t in 1:T_all) { // iterating over reference dates
      for (d in (max(1, D-t+2)):min(D+1, T_all-t+1)) { // iterating over delays
        triangle_rep[d, t+d-1-D] = triangle[d, t];
      }
    }
    return triangle_rep;
  }
  
  /** Takes a reporting triangle by date of reference, and reshapes it into a
  * padded reporting triangle by date of report, excluding left-truncated dates.
  * This is a version optimized for stan, i.e. matrices are only accessed in
  * column-major order and the loop over delays has been vectorized. Note that
  * the dimensions and orientation of the returned matrix are different from the
  * non-optimized version.
  *
  * @param triangle DxT matrix representing a reporting triangle, i.e. columns
  * represent the date of reference and rows represent the delay
  * 
  * @param D maximum delay
  * 
  * @return A (T-D)xT matrix representing a reporting triangle by date
  * of report. Rows represent the date of report (since first day that is
  * not left-truncated) and columns represent the delay, but ragged. That is,
  * for a given reporting date, the sequence of reports with delay 1:(D+1) has
  * some leading and trailing zeros in the matrix. To obtain the overall cases
  * for a given reference date, one can then pre-multiply this matrix with a
  * T-length row vector of ones.
  */
  matrix reporting_triangle_by_report_padded(matrix triangle, int D){
    int T_all = cols(triangle);
    matrix[T_all-D, T_all] triangle_rep = rep_matrix(0, T_all-D, T_all);
    for (t in 1:T_all) { // iterating over reference dates
      triangle_rep[max(t-D, 1):min(t, T_all-D), T_all-t+1]
        = triangle[(max(1, D-t+2)):min(D+1, T_all-t+1), t];
    }
    return triangle_rep;
  }
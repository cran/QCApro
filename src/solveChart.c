# include <R.h>
# include <Rinternals.h>
# include <R_ext/Rdynload.h>


SEXP solveChart(SEXP x, SEXP y) { 
    int *p_x, *p_y, *p_output, i, j, k, xrows, xcols, yrows, ycols, line;
    
    SEXP usage = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(usage, 0, x = coerceVector(x, INTSXP));
    SET_VECTOR_ELT(usage, 1, y = coerceVector(y, INTSXP));
    
    xrows = nrows(x);
    xcols = ncols(x);
    yrows = nrows(y);
    ycols = ncols(y);
    
    int colsums[ycols];
    
    p_x = INTEGER(x);
    p_y = INTEGER(y);
    
    // create the list to be returned to R
    SEXP root = PROTECT(allocVector(VECSXP, 1));
    
    SEXP output = SET_VECTOR_ELT(root, 0, allocVector(INTSXP, xrows));
    p_output = INTEGER(output);
    
    
    for (i = 0; i < xrows; i++) {
        for (k = 0; k < ycols; k++) {
            colsums[k] = 0;
        }
        
        for (j = 0; j < xcols; j++) {
            line = p_x[i + xrows * j];
            
            for (k = 0; k < ycols; k++) {
                colsums[k] += p_y[line + yrows * k];
            }
        }
        
        //Rprintf("i: %d\n", i);
        //for (k = 0; k < ycols; k++) {
        //    Rprintf(" %d", colsums[k]);
        //}
        //Rprintf("\n");
        
        p_output[i] = 1;
        for (k = 0; k < ycols; k++) {
            if (colsums[k] == 0) {
                p_output[i] = 0;
            }
        }
    }
    
    UNPROTECT(2);
    
    return(root);
}




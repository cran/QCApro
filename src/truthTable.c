# include <R.h>
# include <Rinternals.h>
# include <R_ext/Rdynload.h>


SEXP truthTable(SEXP x, SEXP y, SEXP fuz, SEXP vo) { 
    int i, j, k, index;
    double *p_x, *p_inclpri, *p_which, *p_vo, min, so, sumx, sumpmin, prisum, temp1, temp2;
    int xrows, xcols, yrows, ncut, *p_y, *p_fuz;
    
    SEXP usage = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(usage, 0, x = coerceVector(x, REALSXP));
    SET_VECTOR_ELT(usage, 1, y = coerceVector(y, INTSXP));
    SET_VECTOR_ELT(usage, 2, fuz = coerceVector(fuz, INTSXP));
    SET_VECTOR_ELT(usage, 3, vo = coerceVector(vo, REALSXP));
    
    xrows = nrows(x);
    yrows = nrows(y);
    xcols = ncols(x);
    
    double copyline[xcols];
    
    p_x = REAL(x);
    p_y = INTEGER(y);
    p_fuz = INTEGER(fuz);
    p_vo = REAL(vo);
    
    // create the list to be returned to R
    SEXP root = PROTECT(allocVector(VECSXP, 2));
    SEXP inclpri = SET_VECTOR_ELT(root, 0, allocMatrix(REALSXP, 3, yrows));
    SEXP which = SET_VECTOR_ELT(root, 1, allocVector(REALSXP, xrows));
    
    p_inclpri = REAL(inclpri);
    p_which = REAL(which);
   
    /* generate the expressions' line numbers
    int expressions[yrows];
    for (i = 0; i < yrows; i++) {
        expressions[i] = i + 2;
    }
    */
    
    so = 0;
    
    // sum of the outcome variable
    for (i = 0; i < length(vo); i++) {
        so += p_vo[i];
    }
    
    
    for (k = 0; k < yrows; k++) { // loop for every line of the truth table matrix
        
        sumx = 0;
        sumpmin = 0;
        prisum = 0;  
        ncut = 0;
        
        for (i = 0; i < xrows; i++) { // loop over every line of the data matrix
            
            min = 1000;
            for (j = 0; j < xcols; j++) { // loop over each column of the data matrix
                copyline[j] = p_x[i + xrows * j];
                
                index = k + yrows * j;
                
                if (p_fuz[j] == 1) { // for the fuzzy variables, invert those who have the 3k value equal to 1 ("onex3k" in R)
                    if (p_y[index] == 0) {
                        copyline[j] = 1 - copyline[j];
                    }
                }
                else {
                    if (p_y[index] != (copyline[j])) {
                        copyline[j] = 0;
                    }
                    else {
                        copyline[j] = 1;
                    }
                }
                
                if (copyline[j] < min) {
                    min = copyline[j];
                }
                
            } // end of j loop, over columns
            
            sumx += min;
            sumpmin += (min < p_vo[i])?min:p_vo[i];
            temp1 = (min < p_vo[i])?min:p_vo[i];
            temp2 = 1 - p_vo[i];
            prisum += (temp1 < temp2)?temp1:temp2;
            ncut += (min > 0.5)?1:0;
            
        } // end of i loop
        
        
        //p_inclpri[k*3] = (sumpmin == 0 && sumx == 0)?0:(sumpmin/sumx);
        p_inclpri[k*3] = sumpmin/sumx;
        
        //temp1 = sumpmin - prisum;
        //temp2 = sumx - prisum;
        //p_inclpri[k*3 + 1] = (temp1 == 0 && temp2 == 0)?0:(temp1/temp2);
        p_inclpri[k*3 + 1] = (sumpmin - prisum)/(sumx - prisum);
        p_inclpri[k*3 + 2] = ncut;
        
    } // end of k loop
    
    
    for (i = 0; i < xrows; i++) { // loop over every line of the data matrix
        p_which[i] = 0;
        for (k = 0; k < yrows; k++) { // loop for every line of the truth table matrix
            min = 1000;
            for (j = 0; j < xcols; j++) { // loop over each column of the data matrix
                copyline[j] = p_x[i + xrows * j];
                
                index = k + yrows * j;
                
                if (p_fuz[j] == 1) { // for the fuzzy variables, invert those who have the 3k value equal to 1 ("onex3k" in R)
                    if (p_y[index] == 0) {
                        copyline[j] = 1 - copyline[j];
                    }
                }
                else {
                    if (p_y[index] != (copyline[j])) {
                        copyline[j] = 0;
                    }
                    else {
                        copyline[j] = 1;
                    }
                }
                
                if (copyline[j] < min) {
                    min = copyline[j];
                }
            } // end of j loop, over columns
            
            if (min > 0.5) {
                p_which[i] = k + 1;
            }
        }
        // end of k loop
        
    }
    
    
    UNPROTECT(2);
    
    return(root);
}




SEXP truthTableMem(SEXP x, SEXP noflevels, SEXP mbase, SEXP fuz, SEXP vo) { 
    int i, j, k, index;
    double *p_x, *p_inclpri, *p_which, *p_vo, min, so, sumx, sumpmin, prisum, temp1, temp2;
    int xrows, xcols, yrows, ncut, *pnoflevels, *pmbase, *p_fuz;
    
    SEXP usage = PROTECT(allocVector(VECSXP, 5));
    SET_VECTOR_ELT(usage, 0, x = coerceVector(x, REALSXP));
    SET_VECTOR_ELT(usage, 1, noflevels = coerceVector(noflevels, INTSXP));
    SET_VECTOR_ELT(usage, 2, mbase = coerceVector(mbase, INTSXP));
    SET_VECTOR_ELT(usage, 3, fuz = coerceVector(fuz, INTSXP));
    SET_VECTOR_ELT(usage, 4, vo = coerceVector(vo, REALSXP));
    
    xrows = nrows(x);
    // yrows = nrows(y);
    xcols = ncols(x);
    
    so = 0;
    
    
    
    double copyline[xcols]; 
    
    p_x = REAL(x);
    pnoflevels = INTEGER(noflevels);
    pmbase = INTEGER(mbase);
    p_fuz = INTEGER(fuz);
    p_vo = REAL(vo);
    
    yrows = pnoflevels[0];
    for (i = 1; i < length(noflevels); i++) {
        yrows = yrows * (pnoflevels[i]);
    }
    
    // create the list to be returned to R
    
    SEXP root = PROTECT(allocVector(VECSXP, 2));
    SEXP inclpri = SET_VECTOR_ELT(root, 0, allocMatrix(REALSXP, 3, yrows));
    SEXP which = SET_VECTOR_ELT(root, 1, allocVector(REALSXP, xrows));
    
    p_inclpri = REAL(inclpri);
    p_which = REAL(which);
    
    // sum of the outcome variable
    for (i = 0; i < length(vo); i++) {
        so += p_vo[i];
    }
    
    
    for (k = 0; k < yrows; k++) { // loop for every line of the truth table matrix
        
        sumx = 0;
        sumpmin = 0;
        prisum = 0;  
        ncut = 0;
        index = 1000;
        
        for (i = 0; i < xrows; i++) { // loop over every line of the data matrix
            
            min = 1000;
            for (j = 0; j < xcols; j++) { // loop over each column of the data matrix
                copyline[j] = p_x[i + xrows * j];
                
                index = div(div(k, pmbase[j]).quot, pnoflevels[j]).rem;
                //Rprintf("k: %d  i: %d  j: %d  pmbase[j] %d  noflevels[j] %d      %d\n", k, i, j, pmbase[j], pnoflevels[j], index);
                
                if (p_fuz[j] == 1) { // for the fuzzy variables, invert those who have the 3k value equal to 1 ("onex3k" in R)
                    if (index == 0) {
                        copyline[j] = 1 - copyline[j];
                    }
                }
                else {
                    if (index != (copyline[j])) {
                        copyline[j] = 0;
                    }
                    else {
                        copyline[j] = 1;
                    }
                }
                
                if (copyline[j] < min) {
                    min = copyline[j];
                }
                
            } // end of j loop, over columns
            // Rprintf("\n");
            sumx += min;
            sumpmin += (min < p_vo[i])?min:p_vo[i];
            temp1 = (min < p_vo[i])?min:p_vo[i];
            temp2 = 1 - p_vo[i];
            prisum += (temp1 < temp2)?temp1:temp2;
            ncut += (min > 0.5)?1:0;
            
        } // end of i loop
        
        
        //p_inclpri[k*3] = (sumpmin == 0 && sumx == 0)?0:(sumpmin/sumx);
        p_inclpri[k*3] = sumpmin/sumx;
        
        //temp1 = sumpmin - prisum;
        //temp2 = sumx - prisum;
        //p_inclpri[k*3 + 1] = (temp1 == 0 && temp2 == 0)?0:(temp1/temp2);
        p_inclpri[k*3 + 1] = (sumpmin - prisum)/(sumx - prisum);
        p_inclpri[k*3 + 2] = ncut;
        
    } // end of k loop
    
    
    for (i = 0; i < xrows; i++) { // loop over every line of the data matrix
        p_which[i] = 0;
        for (k = 0; k < yrows; k++) { // loop for every line of the truth table matrix
            min = 1000;
            for (j = 0; j < xcols; j++) { // loop over each column of the data matrix
                copyline[j] = p_x[i + xrows * j];
                
                index = div(div(k, pmbase[j]).quot, pnoflevels[j]).rem;
                
                if (p_fuz[j] == 1) { // for the fuzzy variables, invert those who have the 3k value equal to 1 ("onex3k" in R)
                    if (index == 0) {
                        copyline[j] = 1 - copyline[j];
                    }
                }
                else {
                    if (index != (copyline[j])) {
                        copyline[j] = 0;
                    }
                    else {
                        copyline[j] = 1;
                    }
                }
                
                if (copyline[j] < min) {
                    min = copyline[j];
                }
            } // end of j loop, over columns
            
            if (min > 0.5) {
                p_which[i] = k + 1;
            }
        }
        // end of k loop
        
    }
    
    
    UNPROTECT(2);
    
    return(root);
}


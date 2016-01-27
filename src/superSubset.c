# include <R.h>
# include <Rinternals.h>
# include <R_ext/Rdynload.h>

SEXP superSubset(SEXP x, SEXP y, SEXP fuz, SEXP vo, SEXP nec) { 
    int i, j, k, index;
    double *p_x, *p_incovpri, *p_vo, min, max, so = 0.0, sumx_min, sumx_max, sumpmin_min, sumpmin_max, prisum_min, prisum_max, temp1, temp2;
    int xrows, xcols, yrows, *p_y, *p_fuz, *p_nec;
    
    SEXP usage = PROTECT(allocVector(VECSXP, 5));
    SET_VECTOR_ELT(usage, 0, x = coerceVector(x, REALSXP));
    SET_VECTOR_ELT(usage, 1, y = coerceVector(y, INTSXP));
    SET_VECTOR_ELT(usage, 2, fuz = coerceVector(fuz, INTSXP));
    SET_VECTOR_ELT(usage, 3, vo = coerceVector(vo, REALSXP));
    SET_VECTOR_ELT(usage, 4, nec = coerceVector(nec, INTSXP));
    
    xrows = nrows(x);
    yrows = nrows(y);
    xcols = ncols(x);
    
    double copyline[xcols];
    
    p_x = REAL(x);
    p_y = INTEGER(y);
    p_fuz = INTEGER(fuz);
    p_vo = REAL(vo);
    p_nec = INTEGER(nec);
    
    
    // create the list to be returned to R
    SEXP incovpri = PROTECT(allocMatrix(REALSXP, 6, yrows));
    p_incovpri = REAL(incovpri);
    
    
    // sum of the outcome variable
    for (i = 0; i < length(vo); i++) {
        so += p_vo[i];
    }
    
    
    min = 1000;
    max = 0;
    
    for (k = 0; k < yrows; k++) { // loop for every line of the truth table matrix
        
        sumx_min = 0;
        sumx_max = 0;
        sumpmin_min = 0;
        sumpmin_max = 0;
        prisum_min = 0;  
        prisum_max = 0;
        
        for (i = 0; i < xrows; i++) { // loop over every line of the data matrix
            
            for (j = 0; j < xcols; j++) { // loop over each column of the data matrix
                copyline[j] = p_x[i + xrows * j];
                
                index = k + yrows * j;
                
                if (p_fuz[j] == 1) { // for the fuzzy variables, invert those who have the 3k value equal to 1 ("onex3k" in R)
                    if (p_y[index] == 1) {
                        copyline[j] = 1 - copyline[j];
                    }
                }
                else {
                    if (p_y[index] != (copyline[j] + 1)) {
                        copyline[j] = 0;
                    }
                    else {
                        copyline[j] = 1;
                    }
                }
                
                if (p_y[index] != 0) {
                    
                    if (copyline[j] < min) {
                        min = copyline[j];
                    }
                    
                    if (copyline[j] > max) {
                        max = copyline[j];
                    }
                }
                
            } // end of j loop, over columns
            
            sumx_min += min;
            sumx_max += max;
            sumpmin_min += (min < p_vo[i])?min:p_vo[i];
            sumpmin_max += (max < p_vo[i])?max:p_vo[i];
            temp1 = (min < p_vo[i])?min:p_vo[i];
            temp2 = p_nec[0]?(1 - min):(1 - p_vo[i]);
            prisum_min += (temp1 < temp2)?temp1:temp2;
            temp1 = (max < p_vo[i])?max:p_vo[i];
            temp2 = 1 - max;
            prisum_max += (temp1 < temp2)?temp1:temp2;
            
            min = 1000; // re-initialize min and max values
            max = 0;
            
        } // end of i loop
        
        p_incovpri[k*6] = (sumpmin_min == 0 && sumx_min == 0)?0:(sumpmin_min/sumx_min);
        p_incovpri[k*6 + 1] = (sumpmin_min == 0 && so == 0)?0:(sumpmin_min/so);
        p_incovpri[k*6 + 2] = (sumpmin_max == 0 && sumx_max == 0)?0:(sumpmin_max/sumx_max);
        p_incovpri[k*6 + 3] = (sumpmin_max == 0 && so == 0)?0:(sumpmin_max/so);
        
        temp1 = sumpmin_min - prisum_min;
        temp2 = p_nec[0]?so:sumx_min - prisum_min;
        p_incovpri[k*6 + 4] = (temp1 == 0 && temp2 == 0)?0:(temp1/temp2);
        
        temp1 = sumpmin_max - prisum_max;
        temp2 = so - prisum_max;
        p_incovpri[k*6 + 5] = (temp1 == 0 && temp2 == 0)?0:(temp1/temp2);
        
    } // end of k loop
    
    
    UNPROTECT(2);
    
    return(incovpri);
}




SEXP superSubsetMem(SEXP x, SEXP noflevels, SEXP mbase, SEXP fuz, SEXP vo, SEXP nec) { 
    int i, j, k, index;
    double *px, *pincovpri, *pvo, min, max, so = 0.0, sumx_min, sumx_max, sumpmin_min, sumpmin_max, prisum_min, prisum_max, temp1, temp2;
    int xrows, xcols, yrows, *pnoflevels, *pmbase,  *pfuz, *pnec;
    
    SEXP usage = PROTECT(allocVector(VECSXP, 6));
    SET_VECTOR_ELT(usage, 0, x = coerceVector(x, REALSXP));
    SET_VECTOR_ELT(usage, 1, noflevels = coerceVector(noflevels, INTSXP));
    SET_VECTOR_ELT(usage, 2, mbase = coerceVector(mbase, INTSXP));
    SET_VECTOR_ELT(usage, 3, fuz = coerceVector(fuz, INTSXP));
    SET_VECTOR_ELT(usage, 4, vo = coerceVector(vo, REALSXP));
    SET_VECTOR_ELT(usage, 5, nec = coerceVector(nec, INTSXP));
    
    
    px = REAL(x);
    pnoflevels = INTEGER(noflevels);
    pmbase = INTEGER(mbase);
    pfuz = INTEGER(fuz);
    pvo = REAL(vo);
    pnec = INTEGER(nec);
    
    
    yrows = pnoflevels[0] + 1;
    for (i = 1; i < length(noflevels); i++) {
        yrows = yrows * (pnoflevels[i] + 1);
    }
    yrows = yrows - 1;
    
    
    xrows = nrows(x);
    //yrows = nrows(y);
    xcols = ncols(x);
    
    double copyline[xcols];
    
    
    
    // create the list to be returned to R
    SEXP incovpri = PROTECT(allocMatrix(REALSXP, 6, yrows));
    pincovpri = REAL(incovpri);
    
    
    // sum of the outcome variable
    for (i = 0; i < length(vo); i++) {
        so += pvo[i];
    }
    
    
    min = 1000;
    max = 0;
    
    for (k = 0; k < yrows; k++) { // loop for every line of the truth table matrix
        
        sumx_min = 0;
        sumx_max = 0;
        sumpmin_min = 0;
        sumpmin_max = 0;
        prisum_min = 0;  
        prisum_max = 0;
        
        for (i = 0; i < xrows; i++) { // loop over every line of the data matrix
            
            for (j = 0; j < xcols; j++) { // loop over each column of the data matrix
                copyline[j] = px[i + xrows * j];
                
                index = div(div(k + 1, pmbase[j]).quot, pnoflevels[j] + 1).rem;
                
                if (pfuz[j] == 1) { // for the fuzzy variables, invert those who have the 3k value equal to 1 ("onex3k" in R)
                    if (index == 1) {
                        copyline[j] = 1 - copyline[j];
                    }
                }
                else {
                    if (index != (copyline[j] + 1)) {
                        copyline[j] = 0;
                    }
                    else {
                        copyline[j] = 1;
                    }
                }
                
                if (index != 0) {
                    
                    if (copyline[j] < min) {
                        min = copyline[j];
                    }
                    
                    if (copyline[j] > max) {
                        max = copyline[j];
                    }
                }
                
            } // end of j loop, over columns
            
            sumx_min += min;
            sumx_max += max;
            sumpmin_min += (min < pvo[i])?min:pvo[i];
            sumpmin_max += (max < pvo[i])?max:pvo[i];
            temp1 = (min < pvo[i])?min:pvo[i];
            temp2 = pnec[0]?(1 - min):(1 - pvo[i]);
            prisum_min += (temp1 < temp2)?temp1:temp2;
            temp1 = (max < pvo[i])?max:pvo[i];
            temp2 = 1 - max;
            prisum_max += (temp1 < temp2)?temp1:temp2;
            
            min = 1000; // re-initialize min and max values
            max = 0;
            
        } // end of i loop
        
        pincovpri[k*6] = (sumpmin_min == 0 && sumx_min == 0)?0:(sumpmin_min/sumx_min);
        pincovpri[k*6 + 1] = (sumpmin_min == 0 && so == 0)?0:(sumpmin_min/so);
        pincovpri[k*6 + 2] = (sumpmin_max == 0 && sumx_max == 0)?0:(sumpmin_max/sumx_max);
        pincovpri[k*6 + 3] = (sumpmin_max == 0 && so == 0)?0:(sumpmin_max/so);
        
        temp1 = sumpmin_min - prisum_min;
        temp2 = pnec[0]?so:sumx_min - prisum_min;
        pincovpri[k*6 + 4] = (temp1 == 0 && temp2 == 0)?0:(temp1/temp2);
        
        temp1 = sumpmin_max - prisum_max;
        temp2 = so - prisum_max;
        pincovpri[k*6 + 5] = (temp1 == 0 && temp2 == 0)?0:(temp1/temp2);
        
    } // end of k loop
    
    
    UNPROTECT(2);
    
    return(incovpri);
}


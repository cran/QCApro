# include <R.h>
# include <Rinternals.h>
# include <R_ext/Rdynload.h>

SEXP removeRedundants(SEXP rowno, SEXP noflevels, SEXP mbase) {
    
    int *pointer_next, *pointer_final, *pointer_temp1, *pointer_temp2, *pointer_rowno, *pointer_noflevels, *pointer_mbase;
    int previous, lmbase, ltemp2, lrowno, lmbasei, i, j, k, rn, finalength, lungime, flag2, flag1, templung;
    SEXP next, final, temp1, temp2;
    
    SEXP usage = PROTECT(allocVector(VECSXP, 7));
    
    SET_VECTOR_ELT(usage, 0, rowno = coerceVector(rowno, INTSXP));
    SET_VECTOR_ELT(usage, 1, noflevels = coerceVector(noflevels, INTSXP));
    SET_VECTOR_ELT(usage, 2, mbase = coerceVector(mbase, INTSXP));
    
    pointer_rowno = INTEGER(rowno);
    pointer_noflevels = INTEGER(noflevels);
    pointer_mbase = INTEGER(mbase);
    lmbase = length(mbase);
    lrowno = length(rowno);
    
    
    // create the vector of next positions (similar to linked lists)
	SET_VECTOR_ELT(usage, 3, next = allocVector(INTSXP, lrowno));
	pointer_next = INTEGER(next);
	
	for (i = 0; i < lrowno; i++) {
	    pointer_next[i] = i + 1;
	}
	
    
	rn = 0;
	flag1 = 0;
	
	while (rn < lrowno) {
	    templung = 1;
	    previous = rn;
        
        SET_VECTOR_ELT(usage, 4, temp1 = allocVector(INTSXP, 1));
        pointer_temp1 = INTEGER(temp1);
        pointer_temp1[0] = pointer_rowno[rn];
        
        
        flag2 = 0;
            
        for (i = 0; i < lmbase; i++) {
            lmbasei = lmbase - i - 1;
            if (div(div(pointer_rowno[rn] - 1, pointer_mbase[lmbasei]).quot, pointer_noflevels[lmbasei] + 1).rem == 0) {
                flag2 = 1;
                lungime = templung * (pointer_noflevels[lmbasei] + 1);
                
                SET_VECTOR_ELT(usage, 5, temp2 = allocVector(INTSXP, lungime));
                pointer_temp2 = INTEGER(temp2);
                
                for (j = 0; j < length(temp1); j++) {
                    pointer_temp2[j] = pointer_temp1[j];
                    for (k = 0; k < pointer_noflevels[lmbasei]; k++) {
                        pointer_temp2[j + length(temp1)*(k + 1)] = pointer_temp1[j] + (k + 1)*pointer_mbase[lmbasei];
                    }
                }
                
                if (i < lmbase) {
                    SET_VECTOR_ELT(usage, 4, temp1 = allocVector(INTSXP, lungime));
                    pointer_temp1 = INTEGER(temp1);
                    
                    for (j = 0; j < lungime; j++) {
                        pointer_temp1[j] = pointer_temp2[j];
                    }
                    templung = lungime;
                }
            }
        }
        
        if (flag2 == 1) { // the current row number has subsets, stored in temp2
            ltemp2 = length(temp2);
            
            i = pointer_next[previous];
            j = 0;
            while (i < lrowno && j < ltemp2) {
                if (pointer_rowno[i] < pointer_temp2[j]) {
                    previous = i;
                    i = pointer_next[i];
                }
                else if (pointer_rowno[i] > pointer_temp2[j]) {
                    j++;
                }
                else { // if (pointer_rowno[i] == pointer_temp2[j])
                    flag1 = 1;
                    pointer_next[previous] = pointer_next[i];
                    i = pointer_next[i];
                    j++;
                }
            }
        }
        
		rn = pointer_next[rn];
    }
        
    
    if (flag1 == 0) { // no single row has been found redundant
        UNPROTECT(1);
        return(rowno);
    }
    else {
    
        finalength = 0;
        i = 0;
        while (i < lrowno) {
            i = pointer_next[i];
            finalength++;
        }
        
        SET_VECTOR_ELT(usage, 6, final = allocVector(INTSXP, finalength));
        pointer_final = INTEGER(final);
        
        i = 0;
        j = 0;
        while (i < lrowno) {
            pointer_final[j] = pointer_rowno[i];
            i = pointer_next[i];
            j += 1; 
        }
        
        UNPROTECT(1);
        return(final);
    }
}


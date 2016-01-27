# include <R.h>
# include <Rinternals.h>
# include <R_ext/Rdynload.h>
# include <stdlib.h>


SEXP findSubsets(SEXP rowno, SEXP noflevels, SEXP mbase, SEXP max) {
    int *prowno, *pnoflevels, *pmbase, *pmax, lmbase, lmbasei, i, j, k, lungime, flag, templung, *ptemp1, *ptemp2;
    SEXP temp1, temp2;
    
    
    
    SEXP usage = PROTECT(allocVector(VECSXP, 6));
    
    SET_VECTOR_ELT(usage, 0, rowno = coerceVector(rowno, INTSXP));
    SET_VECTOR_ELT(usage, 1, noflevels = coerceVector(noflevels, INTSXP));
    SET_VECTOR_ELT(usage, 2, mbase = coerceVector(mbase, INTSXP));
    
    prowno = INTEGER(rowno);
    pnoflevels = INTEGER(noflevels);
    pmbase = INTEGER(mbase);
    
    if (max == R_NilValue) {
        SET_VECTOR_ELT(usage, 3, max = allocVector(INTSXP, 1));
        pmax = INTEGER(max);
        pmax[0] = prowno[length(rowno) - 1];
    }
    else {
        SET_VECTOR_ELT(usage, 3, max = coerceVector(max, INTSXP));
        pmax = INTEGER(max);
    }
    
    SET_VECTOR_ELT(usage, 4, temp1 = allocVector(INTSXP, 1));
    ptemp1 = INTEGER(temp1);
    
    ptemp1[0] = prowno[0];
    flag = 0;
    lmbase = length(mbase);
    
    templung = 1;
    
    
    //Rprintf("length(mbase): %d\n", length(mbase));
    
    for (i = 0; i < lmbase; i++) {
        lmbasei = lmbase - i - 1;
        //Rprintf("rowno: %d; mbase.val: %d; noflevels.i: %d\n", prowno[0], pmbase[i], pnoflevels[i] + 1);
        //Rprintf("%d d/d %d = %d\n", prowno[0] - 1, pmbase[i], div(prowno[0] - 1, pmbase[i]).quot);
        //Rprintf("%d dd %d = %d\n", div(prowno[0] - 1, pmbase[i]).quot, pmbase[i], div(div(prowno[0] - 1, pmbase[i]).quot, pnoflevels[i] + 1).rem);
        if (div(div(prowno[0] - 1, pmbase[lmbasei]).quot, pnoflevels[lmbasei] + 1).rem == 0) {
            flag = 1;
            lungime = templung * (pnoflevels[lmbasei] + 1);
            //Rprintf("lungime: %d\n", lungime);
            SET_VECTOR_ELT(usage, 5, temp2 = allocVector(INTSXP, lungime));
            ptemp2 = INTEGER(temp2);
            
            for (j = 0; j < length(temp1); j++) {
                ptemp2[j] = ptemp1[j];
                for (k = 0; k < pnoflevels[lmbasei]; k++) {
                    ptemp2[j + length(temp1)*(k + 1)] = ptemp1[j] + (k + 1)*pmbase[lmbasei];
                }
            }
            
            if (i < length(mbase)) {
                SET_VECTOR_ELT(usage, 4, temp1 = allocVector(INTSXP, lungime));
                ptemp1 = INTEGER(temp1);
                
                for (j = 0; j < lungime; j++) {
                    ptemp1[j] = ptemp2[j];
                }
                templung = lungime;
            }
        }
    }
    
    //Rprintf("lungime: %d; templung: %d\n\n", lungime, templung);
    
    if (flag == 1) {
    
        templung = 0;
        
        for (i = 0; i < lungime; i++) {
            //Rprintf("%d ", ptemp2[i]);
            if (ptemp2[i] < (pmax[0] + 1)) {
                templung += 1;
            }
        }
        
        //Rprintf("lungime: %d; templung: %d; i: %d\n", lungime, templung, i);
        
        SET_VECTOR_ELT(usage, 4, temp1 = allocVector(INTSXP, templung - 1)); //"- 1" because the first element needs to be eliminated as well
        ptemp1 = INTEGER(temp1);
        
        j = 0;
        for (i = 1; i < lungime; i++) {
            if (ptemp2[i] < pmax[0] + 1) {
                ptemp1[j] = ptemp2[i];
                j += 1;
            }
        }
    }
    else {
        UNPROTECT(1);
        return(R_NilValue);
    }
    
    UNPROTECT(1);
    
    return(temp1);
}

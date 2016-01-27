# include <R.h>
# include <Rinternals.h>
# include <Rmath.h>
# include <R_ext/Rdynload.h>


SEXP allSol(SEXP k, SEXP mtrx) {
    
    SEXP root, temp1, temp2;
    int *p_k, *p_mtrx, *p_output, lmbasej, x, y, mtrxrows, mtrxcols, base2count, totalcount, flag;
    int *p_work, i, j, *p_temp1, *p_temp2, lungime, lengthtemp1, *p_rn; // long long int if possible
    
    // SEXP     Rf_allocVector(SEXPTYPE, R_xlen_t);
    
    PROTECT(root = Rf_allocVector(VECSXP, 7));
    
    SET_VECTOR_ELT(root, 0, k = coerceVector(k, INTSXP));
    SET_VECTOR_ELT(root, 1, mtrx = coerceVector(mtrx, INTSXP));
    p_k = INTEGER(k);
    p_mtrx = INTEGER(mtrx);
    
    mtrxrows = nrows(mtrx);
    mtrxcols = ncols(mtrx);
    
    int colsums[mtrxcols];
    
    // Rprintf("mtrxrows: %d\n", mtrxrows);
    
    // urmatorul cod imi calculeaza puterile lui 2, ex. pentru trei variable 4 2 1
    long long int power[mtrxrows];
    power[mtrxrows - 1] = 1;
    //Rprintf("power: %lld", power[mtrxrows - 1]);
    for (j = 1; j < mtrxrows; j++) {
        power[mtrxrows - j - 1] = R_pow_di(2, j);
        //Rprintf(", %lld", power[mtrxrows - j - 1]);
    }
    
    
    long long int totalrows = power[0]*2;
    //Rprintf("\ntotalrows: %lld\n", totalrows);
    
    SEXP work = SET_VECTOR_ELT(root, 2, allocVector(INTSXP, totalrows));
    p_work = INTEGER(work);
    
    // initialize the whole vector with the value of zero
    for (i = 0; i < totalrows; i++) {
        p_work[i] = 0;
    }
    
    totalcount = 0;
        
    int base2row[mtrxrows];
    
    // Rprintf("totalrows: %d\n", totalrows);
    
    // from i = 1, because the first element in any base 2 combination is always zero
    for (i = 1; i < totalrows; i++) {
        
        // some of the lines will already be redundant and flaged with -1
        if (p_work[i] >= 0) {
            
            base2count = 0;
            
            // reinitializez colsums (sumele de coloane din PI chart) cu valoarea 0
            for (j = 0; j < mtrxcols; j++) {
                colsums[j] = 0;
            }
            
            //calculate the base 2 row
            for (j = 0; j < mtrxrows; j++) {
                base2row[mtrxrows - j - 1] = div(div(i, power[mtrxrows - j - 1]).quot, 2).rem;
            }
            
            for (j = 0; j < mtrxrows; j++) {
                
                base2count += base2row[j];
                
                if (base2row[j] == 1) {
                    for (x = 0; x < mtrxcols; x++) {
                        colsums[x] += p_mtrx[j + mtrxrows * x];
                    }
                }
            }
            
            flag = 1;
            
            for (x = 0; x < mtrxcols; x++) {
                //Rprintf("%d ", colsums[x]);
                if (colsums[x] == 0) {
                    flag = 0;
                }
            }
            //Rprintf("\n");
            
            // the row number in base 2 must contain at least "k" positions equal to 1
            // "k" is received via R, it is the minimum  number of PIs necessary to cover all PI chart columns
            if (base2count < p_k[0]) {
                flag = 0;
            }
            
            if (flag == 0) {
                p_work[i] = 0;
            }
            else {
                p_work[i] = 1;
                
                //Rprintf("i: %d; work: %d\n", i, p_work[i]);
                //for (x = 0; x < mtrxrows; x++) {
                //    Rprintf("%d ", base2row[x]);
                //}
                
                totalcount += 1;
                
                // aici ar trebui sa fac cu -1 toate subseturile
                // incep cu base2row, unde are valoarea 0... stiu procedura
                
                lungime = 1;
                SET_VECTOR_ELT(root, 3, temp1 = allocVector(INTSXP, lungime));
                p_temp1 = INTEGER(temp1);
                p_temp1[0] = i;
                
                flag = 0; // no subsets
                
                for (j = 0; j < mtrxrows; j++) {
                    lmbasej = mtrxrows - j - 1;
                    if (base2row[lmbasej] == 0) {
                        flag = 1; // this line number HAS subsets
                        lungime = lungime * 2; // everything is base 2, therefore there are at MOST <one> other subset
                                               // the length is therefore multiplied by 2 (including the starting row number)
                        
                        /*
                        Rprintf("lmbasej: %d\ntemp1: ", lmbasej);
                        for (x = 0; x < length(temp1); x++) {
                            Rprintf("%d ", p_temp1[x]);
                        }
                        Rprintf("\n");
                        */
                        
                        SET_VECTOR_ELT(root, 4, temp2 = allocVector(INTSXP, lungime));
                        p_temp2 = INTEGER(temp2);
                        
                        lengthtemp1 = length(temp1);
                        
                        for (x = 0; x < lengthtemp1; x++) {
                            p_temp2[x] = p_temp1[x];
                            p_temp2[x + lengthtemp1] = p_temp1[x] + power[lmbasej];
                        }
                        
                        /*
                        Rprintf("temp2: ");
                        for (x = 0; x < length(temp2); x++) {
                            Rprintf("%d ", p_temp2[x]);
                        }
                        Rprintf("\n");
                        */
                        
                        if (j < mtrxrows) {
                            // copiaza toate valorile din temp2 in temp1
                            SET_VECTOR_ELT(root, 3, temp1 = allocVector(INTSXP, lungime));
                            p_temp1 = INTEGER(temp1);
                            for (x = 0; x < lungime; x++) {
                                p_temp1[x] = p_temp2[x];
                            }
                        }
                    }
                }
                
                
                if (flag == 1) {
                    /*
                    Rprintf("line: %d ( ", p_work[i]);
                            
                    for (x = 0; x < mtrxrows; x++) {
                        Rprintf("%d ", base2row[x]);
                    }
                    Rprintf(")\n");
                    
                    Rprintf("Linii: ");
                    for (x = 0; x < length(temp2); x++) {
                        Rprintf("%d ", p_temp2[x]);
                    }
                    Rprintf("\n");
                    */
                    
                    // x incepe de la 1 pentru ca prima linie este chiar cea de la care am plecat
                    for (x = 1; x < length(temp2); x++) {
                        p_work[p_temp2[x]] = -1;
                    }
                    
                    
                }
            }
            
            
            
        }
        else {
            p_work[i] = 0;
        }
        
        //Rprintf("i: %d; work: %d\n", i, p_work[i]);
    }
    
    
    // b2rs base2row sums
    int b2rs[totalcount];
    
    // rn row numbers
    SEXP rn = SET_VECTOR_ELT(root, 5, allocVector(INTSXP, totalcount));
    SEXP output = SET_VECTOR_ELT(root, 6, allocMatrix(INTSXP, mtrxrows, totalcount));
    
    p_rn = INTEGER(rn);
    p_output = INTEGER(output);
    
    if (totalcount > 1) { // all this (sorting, etc.) makes sense only if there are at least two solutions
        x = 0;
        for (i = 0; i < totalrows; i++) {
            if (p_work[i] == 1) {
                p_rn[x] = i;
                b2rs[x] = 0;
                for (j = 0; j < mtrxrows; j++) {
                    b2rs[x] += div(div(i, power[mtrxrows - j - 1]).quot, 2).rem; // this is either 0 or 1
                }
                //Rprintf("i: %d; rn: %d; b2rs: %d\n", i, p_rn[x], b2rs[x]);
                x += 1;
            }
        }
        
        // now sort the row numbers according to the base 2 counts
        // x is the actual position
        // y is the swap
        
        for (i = 0; i < totalcount - 1; i++) {
            //Rprintf("rn: %d\n", p_rn[i]);
            //Rprintf("b2rc: %d\n", b2rs[i]);
            
            x = i;
            for (j = i + 1; j < totalcount; j++) {
                if (b2rs[i] > b2rs[j]) {
                    x = j;
                }
            }
            
            if (x != i) {
                y = b2rs[i];
                b2rs[i] = b2rs[x];
                b2rs[x] = y;
                y = p_rn[i];
                p_rn[i] = p_rn[x];
                p_rn[x] = y;
            }
        }
        
        
        
        
        // sort the line numbers as well, according to the base 2 counts
        for (i = 0; i < totalcount - 1; i++) {
            x = i;
            for (j = i + 1; j < totalcount; j++) {
                if (b2rs[i] == b2rs[j] && p_rn[i] > p_rn[j]) {
                    x = j;
                }
            }
            
            if (x != i) {
                y = p_rn[i];
                p_rn[i] = p_rn[x];
                p_rn[x] = y;
            }
        }
        
    }
    else {
        for (i = 0; i < totalrows; i++) {
            if (p_work[i] == 1) {
                p_rn[0] = i; // there is only one solution, therefore [0]
            }
        }
    }
    
    
    x = 0;
    
    for (i = 0; i < totalcount; i++) {
        //Rprintf("rn: %d\n", p_rn[i]);
        for (j = 0; j < mtrxrows; j++) {
            p_output[x] = div(div(p_rn[i], power[j]).quot, 2).rem; // this is either 0 or 1
            if (p_output[x] == 1) {
                p_output[x] = j + 1;
            }
            x += 1;
        }
    }
    
    UNPROTECT(1);
    return(output);
    //return(R_NilValue);
}

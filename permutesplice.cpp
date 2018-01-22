#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <cmath>

int main(int argc, char *argv[]) {
  char *cfgline, *firstdatafile, *seconddatafile, *outfile, *infile, *left_abutment, *right_abutment, **names;
  char *varname, *value, *filler, *line, *inputformat, *currmarker, *arb, **chdata;
  int k, nperm, cllength, i, j, ia, ib, ic, id, ie, iF, ig, ih, ii, ij, ierror, getflag, plantsTotal, **nMissing;
  int nfac, m, n, *nf, ***coeffs, *savedk;
  double **recombfrac, minmaplength, currmaplength, *distances;
  FILE *CFG, *DAT, *OUT;
  cllength = 512;
  if (argc == 3) {
    k = (int)strtol(argv[2], (char**)NULL, 10);
    if (k > cllength) {
      cllength = k;
      printf("Maximum file name length cllength has been increased from 512 to %d\n", cllength);
    }
  }
  cfgline = (char *) calloc(cllength + 2, sizeof(char));
  line = (char *) calloc(cllength + 2, sizeof(char));
  arb = (char *) calloc(cllength + 2, sizeof(char));
  firstdatafile = (char *) calloc(cllength + 2, sizeof(char));
  seconddatafile = (char *) calloc(cllength + 2, sizeof(char));
  outfile = (char *) calloc(cllength + 2, sizeof(char));
  infile = (char *) calloc(cllength + 2, sizeof(char));
  left_abutment = (char *) calloc(cllength + 2, sizeof(char));
  right_abutment = (char *) calloc(cllength + 2, sizeof(char));
  varname = (char *) calloc(cllength + 2, sizeof(char));
  value = (char *) calloc(cllength + 2, sizeof(char));
  inputformat = (char *) calloc(cllength + 2, sizeof(char));
  if(strcmp(argv[1], "-") == 0 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-help") == 0) {
    printf("Usage: %s configuration_file_name [configuration_file_line_length]\n", argv[0]);
    printf("The configuration file name is mandatory and the configuration file line length is optional.\n");
    printf("Specify the configuration_file_line_length as an integer if it exceeds 512.\n");
    printf("Comments in the configuration file begin with # or // and can follow data if separated by one or more spaces.\n");
    printf("Data lines can appear in any order.  Filenames need not be quoted and can include the full path.\n");
    printf("Data lines consist of space-delimited triplets of attribute, \"=\", and value.\n");
    printf("Allowed attributes include datafile, outfile, left_abutment, right_abutment, and marker.\n");
    printf("Generally there will be 10 markers between the abutments.\n");
    return(0);
  }
  nperm = 0;
  CFG = fopen(argv[1], "r");
  while (fgets(cfgline, cllength, CFG)) {
    ierror = 1;
    varname = strtok(cfgline, " ,|\t\n");
    filler = strtok(NULL, " ,|\t\n");
    value = strtok(NULL, " ,|\t\n");
    fflush(NULL);
    if (varname[0] == '#') {
      ierror = 0; //Do nothing; it's a comment.
      continue; //Comments work only if they have three words because of fscanf instead of fgets.
    }
    if (varname[0] == '/' && varname[1] == '/') {
      ierror = 0; //This is also a comment.
      continue;
    }
    if (strcmp(varname, "firstdatafile") == 0) {
      strcpy(firstdatafile, value);
      ierror = 0;
    }
    if (strcmp(varname, "seconddatafile") == 0) {
      strcpy(seconddatafile, value);
      ierror = 0;
    }
    if (strcmp(varname, "outfile") == 0) {
      strcpy(outfile, value);
      ierror = 0;
    }
    if (strcmp(varname, "left_abutment") == 0) {
      strcpy(left_abutment, value);
      ierror = 0;
    }
    if (strcmp(varname, "right_abutment") == 0) {
      strcpy(right_abutment, value);
      ierror = 0;
    }
    if (strcmp(varname, "inputformat") == 0) {
      if (strcmp(value, "mapmaker") == 0 || strcmp(value, "Mapmaker") == 0) {
        strcpy(inputformat, value);
        ierror = 0;
      }
      else { //Deliberately leave ierror unset to zero.
        printf("Only the Mapmaker input format is currently supported.\n");
      }
    }
    if (strcmp(varname, "marker") == 0) {
      nperm++;
      ierror = 0;
    }
    if (ierror != 0) {
      printf("A variable name was misspelled or not supported.\n");
      return(1);
    }
  }
  names = (char **) malloc((3 + nperm) * sizeof(char *));
  for (i = 0; i < nperm + 3; i++)  names[i] = (char *) calloc(cllength, sizeof(char));
  fseek(CFG, 0, 0);
  k = 1;
  while (fgets(cfgline, cllength, CFG)) {
    ierror = 1;
    varname = strtok(cfgline, " ,|\t\n");
    filler = strtok(NULL, " ,|\t\n");
    value = strtok(NULL, " ,|\t\n");
    if (strcmp(varname, "marker") == 0) {
      strcpy(names[k], value);
      k++;
    }
    if (strcmp(varname, "left_abutment") == 0) strcpy(names[0], value);
    if (strcmp(varname, "right_abutment") == 0) strcpy(names[nperm+1], value);
  }
  OUT = fopen(outfile, "w");
  fprintf(OUT, "first data file = %s\nsecond data file = %s\n", firstdatafile, seconddatafile);
  fprintf(OUT, "output file = %s\n", outfile);
  fprintf(OUT, "left abutment = %s\n", names[0]);
  for (i = 1; i < nperm + 1; i++) {fprintf(OUT, "junction marker %d = %s\n", i, names[i]);}
  fprintf(OUT, "right abutment = %s\n", names[nperm+1]);
  chdata = (char **) malloc((3 + nperm) * sizeof(char *));
  for (i = 0; i < nperm + 3; i++) chdata[i] = (char *) calloc(cllength, sizeof(char));
  recombfrac = (double **) malloc((3 + nperm) * sizeof(double *));
  for (i = 0; i < nperm + 3; i++) recombfrac[i] = (double *) calloc(cllength, sizeof(double));
  nMissing = (int **) malloc((3 + nperm) * sizeof(int *));
  for (i = 0; i < nperm + 3; i++) nMissing[i] = (int *) calloc(cllength, sizeof(int));
  plantsTotal = -1;
  DAT = fopen(firstdatafile, "r");
  fgets(line, cllength, DAT);
  fgets(line, cllength, DAT);
  fgets(line, cllength, DAT); //Skip the three header lines.
  while (fgets(line, cllength, DAT)) { 
    //Currently the user must make sure that cllength exceeds the maximum data line length.
    //Otherwise the data line might be truncated.
    //The genotypes are likewise assumed to be concatenated with no space in between mapping individuals.
    //*g63
    //HBAAHAHAAHABHBAHAHBBABAABABHAABBHABAHHAHABABBBAABHAHHHBABBAAHHAAAAHHABHBAABBBHBABAAHAHAHAAAAAHAABAAHHBBHBAAABAAHABBBAAAAHBHAHBAAABABBAAABBBBAHBBABBHBBBABHABAHHABAHHABAABBHABABBBHBBAABHHAHABAAAAHAAABHAAHBBHHBBABAABBBBAAABABAAABHBAABBBBBHHHABBAABBBABHAHBHAHAAABHAHBHBABHHBAABBHAHAABABABAAAAAHHHABABHHHBAAHBBBHAAAAHHBAAAAAHHAABHABAHHHBAHBBHBAABAHHABHBBBAAAAAABBAAABBBBBHBABBHHHABBBBHBHAAABHHAAAAABBHAHBBHBABBHABBBBAHBBABBABBHBABBHHABAHBBAHHBHABABBBHHAHABBHBAAAHAAHAABABAAAAHABBHBBAHABABAHABAABBHABBABBBB
    value = strtok(line, "|\t\n");
    if (value[0] == '*') {
      arb = value + sizeof(char);
      strcpy(currmarker, arb);
      getflag = 0;
      if (strcmp(currmarker, left_abutment) == 0) {
        k = 0;
        getflag = 1;
        //fprintf(OUT, "diag: current marker is left_abutment is %s\n", currmarker);
      }
      for (i = 1; i <= nperm; i++) {
        if (strcmp(currmarker, names[i]) == 0) {
          k = i;
          getflag = 1;
          //fprintf(OUT, "diag: current marker is %s\n", currmarker);
        }
      }
      if (strcmp(currmarker, right_abutment) == 0) {
        k = nperm + 1;
        getflag = 1;
        //fprintf(OUT, "diag: current marker is right_abutment is %s\n", currmarker);
      }
    }
    else if (getflag == 1) {
      if (plantsTotal < 0) plantsTotal = strlen(value);
      else if (plantsTotal != strlen(value)) {
        fprintf(OUT, "The data lines differ in length.  plantsTotal = %d strlen(value) = %lx\n", plantsTotal, strlen(value));
        return(1);
      }
      strcpy(chdata[k], value);
      //fprintf(OUT, "k = %d current marker = %s data = %s\n", k, currmarker, chdata[k]);
    }
  }
  fclose(DAT);
  DAT = fopen(seconddatafile, "r");
  fgets(line, cllength, DAT);
  fgets(line, cllength, DAT);
  fgets(line, cllength, DAT); //Skip the three header lines.
  while (fgets(line, cllength, DAT)) { 
    //Currently the user must make sure that cllength exceeds the maximum data line length.
    //Otherwise the data line might be truncated.
    //The genotypes are likewise assumed to be concatenated with no space in between mapping individuals.
    //*g63
    //HBAAHAHAAHABHBAHAHBBABAABABHAABBHABAHHAHABABBBAABHAHHHBABBAAHHAAAAHHABHBAABBBHBABAAHAHAHAAAAAHAABAAHHBBHBAAABAAHABBBAAAAHBHAHBAAABABBAAABBBBAHBBABBHBBBABHABAHHABAHHABAABBHABABBBHBBAABHHAHABAAAAHAAABHAAHBBHHBBABAABBBBAAABABAAABHBAABBBBBHHHABBAABBBABHAHBHAHAAABHAHBHBABHHBAABBHAHAABABABAAAAAHHHABABHHHBAAHBBBHAAAAHHBAAAAAHHAABHABAHHHBAHBBHBAABAHHABHBBBAAAAAABBAAABBBBBHBABBHHHABBBBHBHAAABHHAAAAABBHAHBBHBABBHABBBBAHBBABBABBHBABBHHABAHBBAHHBHABABBBHHAHABBHBAAAHAAHAABABAAAAHABBHBBAHABABAHABAABBHABBABBBB
    value = strtok(line, "|\t\n");
    if (value[0] == '*') {
      arb = value + sizeof(char);
      strcpy(currmarker, arb);
      getflag = 0;
      if (strcmp(currmarker, left_abutment) == 0) {
        k = 0;
        getflag = 1;
        //fprintf(OUT, "diag: current marker is left_abutment is %s\n", currmarker);
      }
      for (i = 1; i <= nperm; i++) {
        if (strcmp(currmarker, names[i]) == 0) {
          k = i;
          getflag = 1;
          //fprintf(OUT, "diag: current marker is %s\n", currmarker);
        }
      }
      if (strcmp(currmarker, right_abutment) == 0) {
        k = nperm + 1;
        getflag = 1;
        //fprintf(OUT, "diag: current marker is right_abutment is %s\n", currmarker);
      }
    }
    else if (getflag == 1) {
      strcpy(chdata[k], value);
      //fprintf(OUT, "k = %d current marker = %s data = %s\n", k, currmarker, chdata[k]);
    }
  }
  for (i = 0; i < nperm + 2; i++) {
    for (j = i + 1; j < nperm + 2; j++) {
      nMissing[i][j] = 0;
      recombfrac[i][j] = 0;
      for (k = 0; k < plantsTotal; k++) {
        if (chdata[i][k] == '-' || chdata[j][k] == '-') nMissing[i][j]++;
        else {
          if (chdata[i][k] != chdata[j][k]) recombfrac[i][j]++;
          if (chdata[i][k] == 'A' && chdata[j][k] == 'B') recombfrac[i][j]++; //Add in the second breakpoint for homozygous recombinants.
          if (chdata[i][k] == 'B' && chdata[j][k] == 'A') recombfrac[i][j]++;
        }
      }
      recombfrac[i][j] /= 2; //We have two alleles per locus in this case.
    }
  }
  for (i = 0; i < nperm + 2; i++) {
    for (j = i + 1; j < nperm + 2; j++) {
      recombfrac[i][j] /= (plantsTotal - nMissing[i][j]);
      recombfrac[j][i] = recombfrac[i][j]; //Finally fill in the other half of the recombfrac array.
    }
  }
  for (i = 0; i < nperm + 2; i++) {
    fprintf(OUT, "%f", recombfrac[i][0]);
    for (j = 1; j < nperm + 2; j++) fprintf(OUT, " %f", recombfrac[i][j]);
    fprintf(OUT, "\n");
  }
  nf = (int *) malloc((2 + nperm) * sizeof(int));
  if (nperm <= 10) coeffs = (int ***) malloc((1 + nperm) * sizeof(int **));
  else {
    fprintf(OUT, "This program is currently limited to permuting at most 10 loci.\n");
    return(1);
  }
  coeffs[0] = (int **) malloc(sizeof(int *));
  coeffs[1] = (int **) malloc(sizeof(int *));
  coeffs[2] = (int **) malloc(sizeof(int *));
  for (i = 0; i < 11; i++) nf[i] = 0;
  nfac = 2; //The cases for 0, 1, and 2 loci have been done immediately above.
  for (n = 3; n < 11; n++) {
    m = 0;
    nfac *= n;
    nf[n] = nfac;
    coeffs[n] = (int **) malloc((10 + nfac) * sizeof(int *)); //We are trying to eliminate occasional malloc errors that we do not understand.
    for (j = 0; j < nfac; j++) coeffs[n][j] = (int *) calloc(n, sizeof(int));
    for (ia = 0; ia < n; ia++) {
      for (ib = 0; ib < n; ib++) {
        if (ib == ia) continue;
        for (ic = 0; ic < n; ic++) {
          if (ic == ia || ic == ib) continue;
          if (n == 3) { 
            coeffs[3][m][0] = ia;
            coeffs[3][m][1] = ib;
            coeffs[3][m][2] = ic;
            m++;
          }
          else {
            for (id = 0; id < n; id++) {
              if (id == ia || id == ib || id == ic) continue;
              if (n == 4) {
                coeffs[4][m][0] = ia;
                coeffs[4][m][1] = ib;
                coeffs[4][m][2] = ic;
                coeffs[4][m][3] = id;
                m++;
              }
              else {
                for (ie = 0; ie < n; ie++) {
                  if (ie == ia || ie == ib || ie == ic || ie == id) continue;
                  if (n == 5) {
                    coeffs[5][m][0] = ia;
                    coeffs[5][m][1] = ib;
                    coeffs[5][m][2] = ic;
                    coeffs[5][m][3] = id;
                    coeffs[5][m][4] = ie;
                    m++;
                  }
                  else {
                    for (iF = 0; iF < n; iF++) {
                      if (iF == ia || iF == ib || iF == ic || iF == id || iF == ie) continue;
                      if (n == 6) {
                        coeffs[6][m][0] = ia;
                        coeffs[6][m][1] = ib;
                        coeffs[6][m][2] = ic;
                        coeffs[6][m][3] = id;
                        coeffs[6][m][4] = ie;
                        coeffs[6][m][5] = iF;
                        m++;
                      }
                      else {
                        for (ig = 0; ig < n; ig++) {
                          if (ig == ia || ig == ib || ig == ic || ig == id || ig == ie || ig == iF) continue;
                          if (n == 7) {
                            coeffs[7][m][0] = ia;
                            coeffs[7][m][1] = ib;
                            coeffs[7][m][2] = ic;
                            coeffs[7][m][3] = id;
                            coeffs[7][m][4] = ie;
                            coeffs[7][m][5] = iF;
                            coeffs[7][m][6] = ig;
                            m++;
                          }
                          else {
                            for (ih = 0; ih < n; ih++) {
                              if (ih == ia || ih == ib || ih == ic || ih == id || ih == ie || ih == iF || ih == ig) continue;
                              if (n == 8) {
                                coeffs[8][m][0] = ia;
                                coeffs[8][m][1] = ib;
                                coeffs[8][m][2] = ic;
                                coeffs[8][m][3] = id;
                                coeffs[8][m][4] = ie;
                                coeffs[8][m][5] = iF;
                                coeffs[8][m][6] = ig;
                                coeffs[8][m][7] = ih;
                                m++;
                              }
                              else {
                                for (ii = 0; ii < n; ii++) {
                                  if (ii == ia || ii == ib || ii == ic || ii == id || ii == ie || ii == iF || ii == ig || ii == ih) continue;
                                  if (n == 9) {
                                    coeffs[9][m][0] = ia;
                                    coeffs[9][m][1] = ib;
                                    coeffs[9][m][2] = ic;
                                    coeffs[9][m][3] = id;
                                    coeffs[9][m][4] = ie;
                                    coeffs[9][m][5] = iF;
                                    coeffs[9][m][6] = ig;
                                    coeffs[9][m][7] = ih;
                                    coeffs[9][m][8] = ii;
                                    m++;
                                  }
                                  else {
                                    for (ij = 0; ij < n; ij++) {
                                      if (ij == ia || ij == ib || ij == ic || ij == id || ij == ie || ij == iF || ij == ig || ij == ih || ij == ii) continue;
                                      if (n == 10) {
                                        coeffs[10][m][0] = ia;
                                        coeffs[10][m][1] = ib;
                                        coeffs[10][m][2] = ic;
                                        coeffs[10][m][3] = id;
                                        coeffs[10][m][4] = ie;
                                        coeffs[10][m][5] = iF;
                                        coeffs[10][m][6] = ig;
                                        coeffs[10][m][7] = ih;
                                        coeffs[10][m][8] = ii;
                                        coeffs[10][m][9] = ij;
                                        m++;
                                      }
                                      else if (n > 10) {
                                        fprintf(OUT, "This program is designed for a maximum of 10 permuted loci.  You have %d loci.\n", n);
                                        return(1);
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  printf("nperm = %d nf[%d] = %d\n", nperm, nperm, nf[nperm]);
  fflush(NULL);
  for (i = 3; i <= nperm; i++) {
    for (j = 0; j < nf[i]; j++) {
      for (k = 0; k < i; k++) coeffs[i][j][k] += 1;
    }
  }
//  for (j = 0; j < nf[10]; j++) {
//    fprintf(OUT, "%d", coeffs[10][j][0]);
//    for (k = 1; k < 10; k++) fprintf(OUT, " %d", coeffs[10][j][k]);
//    fprintf(OUT, "\n");
//    fflush(NULL);
//  }
  savedk = (int *) calloc(nperm + 1, sizeof(int));
  distances = (double *) calloc(nperm + 1, sizeof(double));
  minmaplength = 100000;
  for (j = 0; j < nf[nperm]; j++) {
    currmaplength = recombfrac[0][coeffs[nperm][j][1]];
    //fprintf(OUT, "beginning currmap length summation: first step = %f\n", currmaplength);
    //fflush(NULL);
    for (k = 0; k < nperm - 1; k++) {
      //fprintf(OUT, "k = %d first coeff = %d ", k, coeffs[nperm][j][k]);
      //fflush(NULL);
      //fprintf(OUT, "second coeff = %d ", coeffs[nperm][j][k+1]);
      //fflush(NULL);
      ia = coeffs[nperm][j][k];
      ib = coeffs[nperm][j][k+1];
      //fprintf(OUT, "recombfrac[%d][%d] = %f\n", coeffs[nperm][j][k], coeffs[nperm][j][k+1], recombfrac[ia][ib]);
      //fflush(NULL);
      currmaplength += recombfrac[ia][ib];
    }
    currmaplength += recombfrac[ib][nperm+1];
    //fprintf(OUT, "currmaplength = %f\n", currmaplength);
    if (currmaplength < minmaplength) {
      minmaplength = currmaplength;
      for (k = 0; k < nperm; k++) savedk[k] = coeffs[nperm][j][k];
    }
  }
  fprintf(OUT, "savedk[0] = %d", savedk[0]);
  for (k = 1; k < nperm; k++) fprintf(OUT, " %d", savedk[k]);
  fprintf(OUT, "\nminimum map length = %f\n", minmaplength);
  fprintf(OUT, "%s", names[0]);
  for (k = 0; k < nperm; k++) {
    ia = savedk[k];
    fprintf(OUT, " %s", names[ia]);
  }
  fprintf(OUT, " %s\n", names[nperm+1]);
  distances[0] = 0;
  ia = savedk[0];
  distances[1] = recombfrac[0][ia];
  for (k = 0; k < nperm - 1; k++) {
    ia = savedk[k];
    ib = savedk[k+1];
    distances[k+2] = distances[k+1] + recombfrac[ia][ib];
  }
  distances[nperm+1] = distances[nperm] + recombfrac[ib][nperm+1];
  fprintf(OUT, "\nVerdict:\n%s %f\n", names[0], distances[0]);
  for (k = 0; k < nperm; k++) {
    ia = savedk[k];
    fprintf(OUT, "%s %f\n", names[ia], distances[k+1]);
  }
  fprintf(OUT, "%s %f\n", names[nperm+1], distances[nperm+1]);
}

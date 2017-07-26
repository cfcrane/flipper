#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define NSTACK 55

void sort2(int, double[], int);
void sort2(int, double[], int[], int);
FILE *OUT;

int main(int argc, char *argv[]) {
  FILE *FAK, *PRD, *CFG, *OUTA, *RNN;
  char rnstates[260], varname[512], value[512], outfile[512], outafile[512], fakfile[512], prdfile[512], rnnfile[512], outputformat[512];
  int i, j, k, kk, m, lm, igen, times, baseNumber, chimax, *chrlengths, inssortlimit, takeflag, nLoci, *sortedPosIndex, *chromosomes;
  int plants, *preDosage, *postDosage, dicount[2], mindist, savedIndexa, savedIndexb, ngenerations, itemp, uniformityflag, ierror;
  int **srtPosIndex, **srtTempIndex, totalones, totalzeroes, totalall, *ones, *zeroes, *tempPosIndex, rmax, todecouple, breakpoints;
  int **printarray;
  long rseed, secseed;
  double r, total, probnothing, chitotal, *bins, *landmarks, *chipos, *sortedpos, *locusPosition, *tempPosition, breakspercell;
  double *sortedPosition, lastTelomere, nextTelomere, lowerBoundary, upperBoundary, lastChiasma, nextChiasma, increment;
  double probmissing, probwrong, missinglimit;
  CFG = fopen(argv[1], "r");
  i = 1;
  j = 0;
  while (i > 0) {
    ierror = 1;
    //i = fscanf(CFG, "%s %*s %s", &varname, &value);
    i = fscanf(CFG, "%s %*s %s", varname, value);
    if (strcmp(varname, "outputformat") == 0) {
      if (strcmp(value, "bygene") == 0 || strcmp(value, "byplant") == 0) {
        strcpy(outputformat, value);
        ierror = 0;
      }
    }
    if (strcmp(varname, "outfile") == 0) {
      strcpy(outfile, value);
      ierror = 0;
    }
    if (strcmp(varname, "outafile") == 0) {
      strcpy(outafile, value);
      ierror = 0;
    }
    if (strcmp(varname, "fakfile") == 0) {
      strcpy(fakfile, value);
      ierror = 0;
    }
    if (strcmp(varname, "prdfile") == 0) {
      strcpy(prdfile, value);
      ierror = 0;
    }
    if (strcmp(varname, "rnnfile") == 0) {
      strcpy(rnnfile, value);
      ierror = 0;
    }
    if (strcmp(varname, "baseNumber") == 0) {
      baseNumber = (int)strtol(value, (char**)NULL, 10);
      chrlengths = (int *) calloc(baseNumber, sizeof(int));
      ierror = 0;
    }
    if (strcmp(varname, "ngenerations") == 0) {
      ngenerations = (int)strtol(value, (char**)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "chimax") == 0) {
      chimax = (int)strtol(value, (char**)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "times") == 0) {
      times = (int)strtol(value, (char**)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "probnothing") == 0) {
      probnothing = strtod(value, (char**)NULL);
      ierror = 0;
    }
    if (strcmp(varname, "nLoci") == 0) {
      nLoci = (int)strtol(value, (char**)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "rseed") == 0) {
      rseed = (int)strtol(value, (char**)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "secseed") == 0) {
      secseed = (int)strtol(value, (char**)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "uniformityflag") == 0) {
      uniformityflag = (int)strtol(value, (char**)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "todecouple") == 0) {
      todecouple = (int)strtol(value, (char**)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "probmissing") == 0) {
      probmissing = strtod(value, (char**)NULL);
      ierror = 0;
    }
    if (strcmp(varname, "probwrong") == 0) {
      probwrong = strtod(value, (char**)NULL);
      ierror = 0;
    }
    if (strcmp(varname, "chrlength") == 0) {
      chrlengths[j] = (int)strtol(value, (char**)NULL, 10);
      j++;
      ierror = 0;
    }
    if (ierror > 0) {
      printf("Variable name %s or value %s was not recognized, check spelling.\n", varname, value);
      return(1);
    }
  }
  if (strcmp(outputformat, "bygene") == 0) {
    printarray = (int **) calloc(2 * (nLoci + 10), sizeof(int *));
    for (i = 0; i < 2 * nLoci; i++) printarray[i] = (int *)calloc(times + 10, sizeof(int));
  }
  missinglimit = probwrong + probmissing;
  printf("outfile = %s outafile = %s fakfile = %s prdfile = %s rnnfile = %s\n", outfile, outafile, fakfile, prdfile, rnnfile);
  printf("baseNumber = %d chimax = %d times = %d probnothing = %f\n", baseNumber, chimax, times, probnothing);
  printf("nLoci = %d rseed = %ld secseed = %ld uniformityflag = %d todecouple = %d\n", nLoci, rseed, secseed, uniformityflag, todecouple);
  for (i = 0; i < baseNumber; i++) printf("chrlengths %d %d\n", i, chrlengths[i]);
  OUT = fopen(outfile, "w");
  OUTA = fopen(outafile, "w");
  FAK = fopen(fakfile, "w");
  PRD = fopen(prdfile, "w");
  RNN = fopen(rnnfile, "w");
  landmarks = (double *) calloc(baseNumber + 10, sizeof(double));
  bins = (double *) calloc(chimax + 10, sizeof(double));
  chipos = (double *) calloc(chimax + 10, sizeof(double));
  locusPosition = (double *) calloc(nLoci + 10, sizeof(double));
  sortedpos = (double *) calloc(chimax + 10, sizeof(double));
  sortedPosition = (double *) calloc(nLoci + 10, sizeof(double));
  sortedPosIndex = (int *) calloc(nLoci + 10, sizeof(int));
  preDosage = (int *) malloc(2 * (nLoci + 10) * sizeof(int));
  postDosage = (int *) malloc(2 * (nLoci + 10) * sizeof(int));
  srtTempIndex = (int **) malloc(3 * sizeof(int));
  for (i = 0; i < 2; i++) srtTempIndex[i] = (int *) malloc((nLoci + 10) * sizeof(int));
  srtPosIndex = (int **) malloc(3 * sizeof(int));
  for (i = 0; i < 2; i++) srtPosIndex[i] = (int *) malloc((nLoci + 10) * sizeof(int));
  ones = (int *) calloc(2 * nLoci + 10, sizeof(int));
  zeroes = (int *) calloc(2 * nLoci + 10, sizeof(int));
  tempPosition = (double *) calloc(nLoci + 10, sizeof(double));
  tempPosIndex = (int *) calloc(nLoci + 10, sizeof(int));
  chromosomes = (int *) malloc(2 * (nLoci + 10) * sizeof(int));
  for (i = 0; i < 2 * nLoci; i++) chromosomes[i] = -1;
  rmax = 1;
  for (i = 0; i < 31; i++) rmax *= 2;
  rmax--;
  printf("rmax = %d\n", rmax);
  initstate(rseed, rnstates, 256);
  inssortlimit = 7;
  total = 0;
  for (i = 0; i < baseNumber; i++) total += chrlengths[i];
  for (i = 1; i <= baseNumber; i++) landmarks[i] = landmarks[i-1] + (double) chrlengths[i-1] / total;
  fprintf(OUT, "baseNumber = %d, chimax = %d, times = %d, probnothing = %f, nLoci = %d\n", baseNumber, chimax, times, probnothing, nLoci);
  fprintf(OUT, "ngenerations = %d, rseed = %ld, secseed = %ld, uniformityflag = %d, todecouple = %d\n", ngenerations, rseed, secseed, uniformityflag, todecouple);
  fprintf(OUT, "landmarks:\n");
  for (i = 0; i <= baseNumber; i++) fprintf(OUT, "%d %f\n", i, landmarks[i]);
  for (i = 0; i < nLoci; i++) tempPosition[i] = (double) (random()) / rmax;
  if (uniformityflag > 0) {
    for (i = 0; i < nLoci; i++) tempPosIndex[i] = i;
    sort2(nLoci, tempPosition, tempPosIndex, inssortlimit);
    increment = 1 / ((double) nLoci + 2);
    for (i = 0; i < nLoci; i++) locusPosition[tempPosIndex[i]] = (i + 1) * increment;
  }
  else {
    for (i = 0; i < nLoci; i++) locusPosition[i] = tempPosition[i];
  }
  for (i = 0; i < nLoci; i++) fprintf(OUT, "DIAG for locusPosition: %d\t%f\n", i, locusPosition[i]);
  for (i = 0; i < nLoci; i++) {
    sortedPosition[i] = locusPosition[i];
    sortedPosIndex[i] = i;
  }
  sort2(nLoci, sortedPosition, sortedPosIndex, inssortlimit);
  for (i = 0; i < nLoci; i++) {
    itemp = sortedPosIndex[i] + nLoci;
    fprintf(OUT, "%d %d %d %f\n", i, sortedPosIndex[i], itemp, sortedPosition[i]);
  }
  for (j = 0; j < 2; j++) {
    lm = 0;
    for (i = 0; i < nLoci; i++) {
      if (sortedPosition[i] >= landmarks[lm]) {
        fprintf(OUTA, "%d %d %f TELO\n", lm, lm + (j + 1) * nLoci, (landmarks[lm] / 2) + j * .5); /*dummy landmark name, no landmarks in sortedPosIndex!*/
        lm++;
      }
      fprintf(OUTA, "%d %d %f 1\n", sortedPosIndex[i], sortedPosIndex[i] + j * nLoci, (sortedPosition[i] / 2) + j * .5);
    }
  }
  fprintf(OUTA, "%d %d %f TELO\n", lm, lm + nLoci, landmarks[lm]);
  mindist = 10000;
  for (i = 0; i < nLoci - 1; i++) {
    if (abs(sortedPosIndex[i] - sortedPosIndex[i+1]) < mindist) {
      mindist = abs(sortedPosIndex[i] - sortedPosIndex[i+1]);
      savedIndexa = sortedPosIndex[i];
      savedIndexb = sortedPosIndex[i+1];
    }
  }
  for (i = 0; i < 2 * nLoci; i++) preDosage[i] = 1;
  for (i = 0; i < 2 * nLoci; i++) {
    fprintf(PRD, "%d ", preDosage[i]);
    if (i % 50 == 49) fprintf(PRD, "\n");
  }
  for (i = 0; i < 2; i++) dicount[i] = 0;
  if (todecouple > 0) initstate(secseed, rnstates, 256);
  totalones = 0;
  totalzeroes = 0;
  totalall = 0;
  for (i = 0; i < 2 * nLoci; i++) {
    for (j = 0; j < baseNumber; j++) {
      if (sortedPosition[i%nLoci] > landmarks[j]) chromosomes[i] = j;
    }
  }
  breakpoints = 0;
  for (plants = 0; plants < times; plants++) {
    if (plants % 1000 == 0) printf("%d plants are done.\n", plants);
    for (i = 0; i < nLoci; i++) srtPosIndex[0][i] = sortedPosIndex[i];
    for (i = 0; i < nLoci; i++) srtPosIndex[1][i] = nLoci + sortedPosIndex[i];
    for (igen = 0; igen < ngenerations; igen++) {
      for (m = 0; m < 2; m++) {
        chitotal = 0;
        kk = 0;
        for (j = 0; j < chimax; j++) chipos[j] = 0;
        for (j = 0; j < chimax; j++) {
          r = (double) random() / rmax;
          if (r > probnothing) {
            chipos[kk] = (double) random() / rmax;
            kk++;
            chitotal++;
          }
        }
        bins[(int)chitotal-1]++;
        sort2((int) chitotal, chipos, inssortlimit);
        lastTelomere = 0;
        nextTelomere = 0;
        lastChiasma = 0;
        nextChiasma = 0;
        upperBoundary = 0;
        lowerBoundary = 0;
        for (i = 0; i < nLoci; i++) {
          if (sortedPosition[i] > nextTelomere) { /*The next loop will always start a new chromosome in this case.*/
            r = (double) random() / rmax;
            if (r < .5) takeflag = 0;
            else takeflag = 1;
          }
          if (sortedPosition[i] > upperBoundary) {
            takeflag = (takeflag + 1) % 2; /*This will flip each interval within a chromosome and still keep new chromosomes random.*/
            /*This method will cause a minor error if loci i and i-1 lie within the same chromosome and are separated by two chiasmata.
            In this case, i and i-1 will not cosegregate.  This should happen only rarely when there are at least tens of loci per chromosome.*/
            for (j = 0; j <= baseNumber; j++) {
              if (landmarks[j] > sortedPosition[i]) {
                nextTelomere = landmarks[j];
                lastTelomere = landmarks[j-1];
                break;
              }
            }
            if (sortedPosition[i] <= chipos[0]) {
              lastChiasma = 0;
              nextChiasma = chipos[0];
            }
            if (sortedPosition[i] > chipos[0] && sortedPosition[i] <= chipos[(int) chitotal-1]) {
              for (j = 0; j < chitotal; j++) {
                if (chipos[j] > sortedPosition[i]) {
                  lastChiasma = chipos[j-1];
                  nextChiasma = chipos[j];
                  break;
                }
              }
            }
            if (sortedPosition[i] > chipos[(int) chitotal-1]) {
              lastChiasma = chipos[(int) chitotal-1];
              nextChiasma = 1;
            }
            if (fabs(sortedPosition[i] - lastChiasma) < fabs(sortedPosition[i] - lastTelomere)) lowerBoundary = lastChiasma;
            else lowerBoundary = lastTelomere;
            if (fabs(sortedPosition[i] - nextChiasma) < fabs(sortedPosition[i] - nextTelomere)) upperBoundary = nextChiasma;
            else upperBoundary = nextTelomere;
          }
          if (sortedPosition[i] >= lowerBoundary && sortedPosition[i] < upperBoundary) {
            srtTempIndex[m][i] = srtPosIndex[takeflag][i];
            dicount[takeflag]++;
          }
        }
      } /*end of m loop*/
      for (i = 0; i < nLoci; i++) {
        srtPosIndex[0][i] = srtTempIndex[0][i];
        srtPosIndex[1][i] = srtTempIndex[1][i];
      }
    } /*end of igen loop*/
    for (i = 0; i < 2 * nLoci; i++) postDosage[i] = 0;
    for (i = 0; i < 2; i++) {
      for (j = 0; j < nLoci; j++) postDosage[srtPosIndex[i][j]] = 1;
    }
    if (strcmp(outputformat, "bygene") == 0) {
      for (i = 0; i < 2 * nLoci; i++) {
        r = (double) random() / rmax;
        if (r < probwrong) {
          if (postDosage[i] == 0) printarray[i][plants] = 1;
          else printarray[i][plants] = 0;
        }
        else {
          if (r < missinglimit) printarray[i][plants] = 2;
          else printarray[i][plants] = postDosage[i];
        }
      }
    }
    if (strcmp(outputformat, "byplant") == 0) {
      fprintf(FAK, ">Plant %d\n", plants);
      for (i = 0; i < 2 * nLoci; i++) {
        r = (double) random() / rmax;
        if (r < probwrong) {
          if (postDosage[i] == 0) fprintf(FAK, "1 ");
          else fprintf(FAK, "0 ");
        }
        else {
          if (r < missinglimit) fprintf(FAK, "- ");
          else fprintf(FAK, "%d ", postDosage[i]);
        }
        if (i % 50 == 49) fprintf(FAK, "\n");
      }
      if (nLoci % 50 != 0) fprintf(FAK, "\n");
    }
    for (i = 1; i < 2 * nLoci; i++) {
      if (chromosomes[i] == chromosomes[i-1]) {
        if (postDosage[sortedPosIndex[i-1]] != postDosage[sortedPosIndex[i]]) breakpoints++;
      }
    }
    for (i = 0; i < 2 * nLoci; i++) {
      totalall++;
      if (postDosage[i] > 0) {
        totalones++;
        ones[i]++;
      }
      else {
        totalzeroes++;
        zeroes[i]++;
      }
    }
  }
  if (strcmp(outputformat, "bygene") == 0) {
    for (i = 0; i < 2 * nLoci; i++) {
      fprintf(FAK, ">Gene %d\n", i);
      for (j = 0; j < times; j++) {
        if (printarray[i][j] == 2) fprintf(FAK, "- ");
        else fprintf(FAK, "%d ", printarray[i][j]);
        if (j % 50 == 49) fprintf(FAK, "\n");
      }
      if (times % 50 != 0) fprintf(FAK, "\n");
    }
  }
  breakspercell = (double) breakpoints / times;
  fprintf(OUT, "breaks per cell = %f\n", breakspercell);
  for (i = 0; i < chimax; i++) bins[i] /= (2 * times);
  fprintf(OUT, "BIN TOTALS:\n");
  for (i = 0; i < chimax; i++) fprintf(OUT, "%d %f\n", i, bins[i]);
  fprintf(OUT, "XTA POSITIONS IN LAST CELL:\n");
  for (i = 0; i < chimax; i++) if (chipos[i] > 0) fprintf (OUT, "%f ", chipos[i]);
  fprintf(OUT, "\n");
  for (i = 0; i < 2 * nLoci; i++) fprintf(OUT, "%d locus = %d zeroes = %d ones = %d\n",i, sortedPosIndex[i], zeroes[i], ones[i]);
  fprintf(OUT, "mindist, indices = %d %d %d\n", mindist, savedIndexa, savedIndexb);
  fprintf(OUT, "dicount for 0, 1 = %d %d\n", dicount[0], dicount[1]);
  fprintf(OUT, "totalzeroes = %d, totalones = %d, totalall = %d, fraction = %f\n", totalzeroes, totalones, totalall, 
    (double) totalzeroes / totalall);
  fprintf(RNN, "%ld\n", random());
  return(0);
}

void sort2(int n, double arr[], int m) {
  int i, ir, j, k, l, *istack, jstack;
  double a, temp;
  ir = n - 1;
  l = 0;
  jstack = 0;
  istack = (int *) malloc(NSTACK * sizeof(int));
  for (;;) {
    if (ir - l < m) {
      for (j = l + 1; j <= ir; j++) {
        a = arr[j];
        for (i = j - 1; i >= 0; i--) {
          if (arr[i] <= a) break;
          arr[i+1] = arr[i];
        }
        arr[i+1] = a;
      }
      if (!jstack) {
        free(istack);
        return;
      }
      ir = istack[jstack];
      l = istack[jstack-1];
      jstack -= 2;
    }
    else {
      k = (l + ir) / 2;
      temp = arr[k];
      arr[k] = arr[l+1];
      arr[l+1] = temp;
      if (arr[l] > arr[ir]) {
        temp = arr[l];
        arr[l] = arr[ir];
        arr[ir] = temp;
      }
      if (arr[l+1] > arr[ir]) {
        temp = arr[l+1];
        arr[l+1] = arr[ir];
        arr[ir] = temp;
      }
      if (arr[l] > arr[l+1]) {
        temp = arr[l];
        arr[l] = arr[l+1];
        arr[l+1] = temp;
      }
      i = l + 1;
      j = ir;
      a = arr[l+1];
      for (;;) {
        do i++; while (arr[i] < a);
        do j--; while (arr[j] > a);
        if (j < i) break;
        temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
      }
      arr[l+1] = arr[j];
      arr[j] = a;
      jstack += 2;
      if (jstack > NSTACK) fprintf(OUT, "NSTACK is too small in sort2.\n");
      if (ir - i + 1 >= j - 1) {
        istack[jstack] = ir;
        istack[jstack-1] = i;
        ir = j - 1;
      }
      else {
        istack[jstack] = j - 1;
        istack[jstack-1] = l;
        l = i;
      }
    }
  }
}

void sort2(int n, double arr[], int brr[], int m) {
  int i, ir, j, k, l, *istack, jstack, b, itemp;
  double a, temp;
  ir = n - 1;
  l = 0;
  jstack = 0;
  istack = (int *) malloc(NSTACK * sizeof(int));
  for (;;) {
    if (ir - l < m) {
      for (j = l + 1; j <= ir; j++) {
        a = arr[j];
        b = brr[j];
        for (i = j - 1; i >= 0; i--) {
          if (arr[i] <= a) break;
          arr[i+1] = arr[i];
          brr[i+1] = brr[i];
        }
        arr[i+1] = a;
        brr[i+1] = b;
      }
      if (!jstack) {
        free(istack);
        return;
      }
      ir = istack[jstack];
      l = istack[jstack-1];
      jstack -= 2;
    }
    else {
      k = (l + ir) / 2;
      temp = arr[k];
      arr[k] = arr[l+1];
      arr[l+1] = temp;
      itemp = brr[k];
      brr[k] = brr[l+1];
      brr[l+1] = itemp;
      if (arr[l] > arr[ir]) {
        temp = arr[l];
        arr[l] = arr[ir];
        arr[ir] = temp;
        itemp = brr[l];
        brr[l] = brr[ir];
        brr[ir] = itemp;
      }
      if (arr[l+1] > arr[ir]) {
        temp = arr[l+1];
        arr[l+1] = arr[ir];
        arr[ir] = temp;
        itemp = brr[l+1];
        brr[l+1] = brr[ir];
        brr[ir] = itemp;
      }
      if (arr[l] > arr[l+1]) {
        temp = arr[l];
        arr[l] = arr[l+1];
        arr[l+1] = temp;
        itemp = brr[l];
        brr[l] = brr[l+1];
        brr[l+1] = itemp;
      }
      i = l + 1;
      j = ir;
      a = arr[l+1];
      b = brr[l+1];
      for (;;) {
        do i++; while (arr[i] < a);
        do j--; while (arr[j] > a);
        if (j < i) break;
        temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
        itemp = brr[i];
        brr[i] = brr[j];
        brr[j] = itemp;
      }
      arr[l+1] = arr[j];
      arr[j] = a;
      brr[l+1] = brr[j];
      brr[j] = b;
      jstack += 2;
      if (jstack > NSTACK) fprintf(OUT, "NSTACK is too small in sort2.\n");
      if (ir - i + 1 >= j - 1) {
        istack[jstack] = ir;
        istack[jstack-1] = i;
        ir = j - 1;
      }
      else {
        istack[jstack] = j - 1;
        istack[jstack-1] = l;
        l = i;
      }
    }
  }
  free(istack);
}

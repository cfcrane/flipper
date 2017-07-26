#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

int main(int argc, char *argv[]) {
  char *datafile, *prdfile, *outfile, *ngrfile, *outafile, *varname, *value, inputformat[512];
  char *line, *token, *lctoken, *queryfile, vars[5][512], **names, *cfgline, *printstring, *inmapfile;
  char *tempname, mmdatatype[20], *markernamesfile, *rfracfile, *filler, **chdata;
  int savedj, imax, queryarray[20], *querycodes, mapqueries[20], mq, *foretestins, *afttestins;
  int i, j, k, l, nLoci, plantsTotal, baseNumber, ploidy, ierror, ***coeffs, mergearray[10], *foretestatt, *afttestatt;
  int *forePointers, *aftPointers, mappedLoci, nLinkageGroups, lowindex, endstatusk, endstatusm, *exemplars;
  int g, m, n, km, nfac, mini, savedv[20], maxnamelength, iv, ivv, outaflag, *membership, highindex, lowcount, highcount;
  int *fracUsedCode, isfull, maxprdlinelength, maxdatlinelength, *preDosage, *postDosage, **nMissing, *plantsUsed;
  int **orthodata, *tempdata, *sectempdata, found, temp, found1, found2, found3, *highlowstatus, lowBs;
  int *lgStarts, *lgEnds, *lgCounts, *sortedlgCounts, *doneLinkageGroups, nf[12], **actmap, totmeas, highAs, highBs, lowAs;
  int ia, ib, ic, id, ie, iF, ig, ih, ii, ij, ip, stemlocusa, stemlocusb, ik, im, *mafts, *mfores, ir;
  int rfileflag, nameflag, claimedloci, claimedplants, As, Bs, missingmeas;
  int flankfore, flankaft, anneallength, foreendm, aftendm, *insertions, savedfore, savedaft, instodo, pused, savedg;
  int gycount, gzcount, commit, combineparents, cllength, *lgnewCounts, *lgnewStarts;
  int *pool, *lgnewEnds, *sortedlgnewCounts, totcount, currcount, printflag;
  int savedjf, savedja, foredist, aftdist, *foreact, *aftact, *foremap, *aftmap;
  int currmaplocus, curractlocus, lastmaplocus, lastactlocus, lg;
  double mappedlength, physlength;
  double segdisthreshold, *segdises, actlength, maxminrf, *minofrfs, **minlgdists;
  double **recombfrac, mapLength, bandwidth, base, top, queriedmaplength, duration, cudist, tdist;
  double frac, minFracFirst, minFracSecond, *lgLengths, minmaplength, localmaplength, attlength, inslength, trial;
  double rfxatt, rfxins, rfxceiling, tripletlength, delta, mindist, trialdista, trialdistb, *lgnewLengths;
  double temprfi, temprfsj, totlength, lengthpercount, *measurements, *sortedmeasurements, meas;
  FILE *CFG, *OUT, *DAT, *PRD, *NGR, *SPE, *QRY, *ACT, *MNA, *RFF, *PUD, *SGD, *DDP, *MEA, *ACM, *INM;
  if(strcmp(argv[1], "-") == 0 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-help") == 0) {
    printf("Usage: %s configuration_file_name [configuration_file_line_length]\n", argv[0]);
    printf("The configuration file name is mandatory and the configuration file line length is optional.\n");
    printf("Specify the configuration file line length only if it needs to exceed 512 characters.\n");
    printf("Comments in the configuration file begin with # or // and can follow data if separated by one or more spaces.\n");
    printf("Data lines can appear in any order.  Filenames need not be quoted and can include the full path.\n");
    printf("Data lines consist of space-delimited triplets of attribute, \"=\", and value.\n");
    printf("Allowed attributes include datafile, prdfile, outfile, ngrfile, outafile,\n");
    printf(" queryfile, markernamesfile, rfracfile, inputformat, nLoci,\n");
    printf(" plantsTotal, maxnamelength, minFracFirst, minFracSecond, ploidy,\n");
    printf(" baseNumber, maxminrf, rfxceiling, segdisthreshold, maxprdlinelength,\n");
    printf(" maxdatlinelength\n");
    printf(" and combineparents.\n");
    return(0);
  }
  printstring = (char *)calloc(200, sizeof(char));
  clock_t timeused, timeused1, timeused2, currtime;
  CFG = fopen(argv[1], "r");
  timeused = clock();
  outaflag = 0; //This controls output of the length of the actual map order in simulations.
  rfileflag = 0; //This controls whether to read in recombination fractions directly.
  combineparents = 0; //This controls whether to map all parental chromosomes separately (0) or together to make baseNumber linkage groups (1).
  maxprdlinelength = 0; //Otherwise this is not getting set with MapMaker input.
  delta = -1;
  cllength = 512;
  if (argc == 3) {
    k = (int)strtol(argv[2], (char**)NULL, 10);
    if (k > cllength) {
      cllength = k;
      printf("Maximum file name length cllength has been increased from 512 to %d\n", cllength);
    }
  }
  cfgline = (char *) calloc(cllength, sizeof(char));
  datafile = (char *) calloc(cllength + 2, sizeof(char));
  prdfile = (char *) calloc(cllength + 2, sizeof(char));
  outfile = (char *) calloc(cllength + 2, sizeof(char));
  ngrfile = (char *) calloc(cllength + 2, sizeof(char));
  outafile = (char *) calloc(cllength + 2, sizeof(char));
  queryfile = (char *) calloc(cllength + 2, sizeof(char));
  markernamesfile = (char *) calloc(cllength + 2, sizeof(char));
  rfracfile = (char *) calloc(cllength + 2, sizeof(char));
  inmapfile = (char *) calloc(cllength + 2, sizeof(char));
  while (fgets(cfgline, cllength, CFG)) {
    fflush(NULL);
    ierror = 1;
    varname = strtok(cfgline, " ,|\t\n");
    filler = strtok(NULL, " ,|\t\n");
    value = strtok(NULL, " ,|\t\n");
    fflush(NULL);
    if (varname[0] == '#') {
      ierror = 0; //Do nothing; it's a comment.
      continue;
    }
    if (varname[0] == '/' && varname[1] == '/') {
      ierror = 0; //This is also a comment.
      continue;
    }
    if (strcmp(varname, "inputformat") == 0) {
      if (strcmp(value, "bygene") == 0 || strcmp(value, "byplant") == 0 || strcmp(value, "mapmaker") == 0) {
        strcpy(inputformat, value);
        ierror = 0;
      }
      else {
        printf("The value of inputformat must be byplant, bygene, or mapmaker.  Please correct the configuration file.\n");
        return(1);
      }
    }
    if (strcmp(varname, "datafile") == 0) {
      strcpy(datafile, value);
      ierror = 0;
    }
    if (strcmp(varname, "rfracfile") == 0) {
      strcpy(rfracfile, value);
      rfileflag = 1;
      ierror = 0;
    }
    if (strcmp(varname, "inmapfile") == 0) {
      strcpy(inmapfile, value);
      ierror = 0;
    }
    if (strcmp(varname, "prdfile") == 0) {
      strcpy(prdfile, value);
      ierror = 0;
    }
    if (strcmp(varname, "outfile") == 0) {
      strcpy(outfile, value);
      ierror = 0;
    }
    if (strcmp(varname, "markernamesfile") == 0) {
      strcpy(markernamesfile, value);
      ierror = 0;
    }
    if (strcmp(varname, "ngrfile") == 0) {
      strcpy(ngrfile, value);
      ierror = 0;
    }
    if (strcmp(varname, "outafile") == 0) {
      strcpy(outafile, value);
      outaflag = 1;
      ierror = 0;
    }
    if (strcmp(varname, "nLoci") == 0) {
      nLoci = (int)strtol(value, (char**)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "plantsTotal") == 0) {
      plantsTotal = (int)strtol(value, (char**)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "maxnamelength") == 0) {
      maxnamelength = (int)strtol(value, (char**)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "minFracFirst") == 0) {
      minFracFirst = strtod(value, (char**)NULL);
      ierror = 0;
    }
    if (strcmp(varname, "minFracSecond") == 0) {
      minFracSecond = strtod(value, (char**)NULL);
      ierror = 0;
    }
    if (strcmp(varname, "baseNumber") == 0) {
      baseNumber = (int)strtol(value, (char**)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "ploidy") == 0) {
      ploidy = (int)strtol(value, (char**)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "maxminrf") == 0) {
      maxminrf = strtod(value, (char **)NULL);
      ierror = 0;
    }
    if (strcmp(varname, "rfxceiling") == 0) {
      rfxceiling = strtod(value, (char **)NULL);
      ierror = 0;
    }
    if (strcmp(varname, "delta") == 0) {
      delta = strtod(value, (char **)NULL);
      ierror = 0;
    }
    if (strcmp(varname, "segdisthreshold") == 0) {
      segdisthreshold = strtod(value, (char **)NULL);
      ierror = 0;
    }
    if (strcmp(varname, "maxprdlinelength") == 0) {
      maxprdlinelength = (int)strtol(value, (char **)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "maxdatlinelength") == 0) {
      maxdatlinelength = (int)strtol(value, (char **)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "combineparents") == 0) {
      combineparents = (int)strtol(value, (char **)NULL, 10);
      ierror = 0;
    }
    if (ierror == 1) {
      printf("The configuration file has a misspelling, an incorrect value, or an extra variable.  Please fix it.\n");
      printf("The offending line contains: %s\n", varname);
      return(1);
    }
  }
  fclose(CFG);
  if (strlen(outfile) > 1) OUT = fopen(outfile, "w");
  else {
    printf("The configuration file did not specify the output file.\n");
    return(1);
  }
  setbuf(OUT, NULL);
  fprintf(OUT, "The configuration file was %s\n", argv[1]);
  if (minFracSecond > minFracFirst) { //This and other parameter-error messages go here.
    fprintf(OUT, "Inconsistent fraction thresholds, first = %f, second = %f\n", minFracFirst, minFracSecond);
    return(1);
  }
  if (strlen(datafile) > 1) DAT = fopen(datafile, "r");
  else {
    printf("The configuration file did not specify the data file.\n");
    return(1);
  }
  if (strlen(prdfile) > 1) {
    if (strcmp(inputformat, "mapmaker") != 0) PRD = fopen(prdfile, "r");
  } 
  else {
    if (strcmp(inputformat, "mapmaker") != 0) {
      printf("The configuration file did not specify the predosage file.\n");
      return(1);
    }
  }
  if (delta < 0) delta = 0.5 / plantsTotal; //The default delta is the bandwidth.
  /*ngrfile is for opening NGR conditionally below.*/
  fprintf(OUT, "ploidy = %d, baseNumber = %d, nLoci = %d, plantsTotal = %d\n", ploidy, baseNumber, nLoci, plantsTotal);
  fprintf(OUT, "maxprdlinelength = %d, maxdatlinelength = %d\n", maxprdlinelength, maxdatlinelength);
  names = (char **) malloc(nLoci * sizeof(char *));
  for (i = 0; i < nLoci; i++) names[i] = (char *) calloc(maxnamelength + 4, sizeof(char));
  tempname = (char *) calloc(maxnamelength, sizeof(char));
  recombfrac = (double **) malloc(nLoci * sizeof(double *));
  for (i = 0; i < nLoci; i++) recombfrac[i] = (double *) calloc(nLoci, sizeof(double));
  nMissing = (int **) calloc(nLoci, sizeof(int *));
  for (i = 0; i < nLoci; i++) nMissing[i] = (int *) calloc(nLoci, sizeof(int));
  preDosage = (int *) calloc(nLoci, sizeof(int));
  postDosage = (int *) calloc(nLoci, sizeof(int));
  plantsUsed = (int *) calloc(nLoci, sizeof(int));
  fracUsedCode = (int *) calloc(nLoci, sizeof(int));
  minofrfs = (double *) calloc(nLoci, sizeof(double));
  segdises = (double *) calloc(nLoci, sizeof(double));
  insertions = (int *) malloc(nLoci * sizeof(int));
  mfores = (int *) malloc(nLoci * sizeof(int));
  mafts = (int *) malloc(nLoci * sizeof(int));
  membership = (int *) malloc(nLoci * sizeof(int));
  foreact = (int *) malloc(nLoci * sizeof(int));
  aftact = (int *) malloc(nLoci * sizeof(int));
  foremap = (int *) malloc(nLoci * sizeof(int));
  aftmap = (int *) malloc(nLoci * sizeof(int));
  if (maxprdlinelength > maxdatlinelength) {
    line = (char *) calloc(maxprdlinelength, sizeof(char));
    token = (char *) calloc(maxprdlinelength, sizeof(char));
    lctoken = (char *) calloc(maxprdlinelength, sizeof(char));
  }
  else {
    line = (char *) calloc(maxdatlinelength, sizeof(char));
    token = (char *) calloc(maxdatlinelength, sizeof(char));
    lctoken = (char *) calloc(maxdatlinelength, sizeof(char));
  }
  timeused1 = clock();
  duration = (double) (timeused1 - timeused) / CLOCKS_PER_SEC;;
  fprintf(OUT, "done with setup, time used = %f\n", duration);
  if (rfileflag == 1) {
    if (strlen(rfracfile) > 1) RFF = fopen(rfracfile, "r"); //This file also must contain the marker names.
    else {
      printf("The configuration file did not specify the name of the file to hold recombination fractions.\n");
      return(1);
    }
    fprintf(OUT, "The rfracfile %s will be read.\n", rfracfile);
    k = -1;
    fgets(line, maxdatlinelength, RFF); //Skip the header.
    while (fgets(line, maxdatlinelength, RFF)) {
      i = 0;
      while (line[i] == ' ' || line[i] == '\t') i++; //Move to the first non-whitespace character in the line.
      if (line[i] == '/' && line[i+1] == '/') continue;
      if (line[i] == '#') continue;
      if (line[i] == '>') {
        j = i; //Emulate a Perl chomp.
        while (line[j] != '\n') j++;
        line[j] = '\0';
        k++;
        strcpy(names[k], line + (i + 1) * sizeof(char)); //Lop off the initial > sign.
        j = 0;
      }
      else {
        token = strtok(line, " ,|\t\n");
        recombfrac[k][j] = strtod(token, (char**)NULL);
        j++;
        while ((token = strtok(NULL, " ,|\t\n"))) {
          recombfrac[k][j] = strtod(token, (char**)NULL);
          j++;
        }
      }
    }
    while (fgets(line, maxdatlinelength, PUD)) {
      i = 0;
      while (line[i] == ' ' || line[i] == '\t') i++; //Move to the first non-whitespace character in the line.
      if (line[i] == '/' && line[i+1] == '/') continue;
      if (line[i] == '#') continue;
      while (line[i] != '\n') i++; //Emulate a Perl chomp.
      line[i] = '\0';
      token = strtok(line, " ,|\t\n");
      j = (int)strtol(token, (char **)NULL, 10);
      token = strtok(NULL, " ,|\t\n");
      plantsUsed[j] = (int)strtol(token, (char **)NULL, 10);
    }
    fclose(PUD);
    while (fgets(line, maxdatlinelength, SGD)) {
      i = 0;
      while (line[i] == ' ' || line[i] == '\t') i++; //Move to the first non-whitespace character in the line.
      if (line[i] == '/' && line[i+1] == '/') continue;
      if (line[i] == '#') continue;
      while (line[i] != '\n') i++; //Emulate a Perl chomp.
      line[i] = '\0';
      token = strtok(line, " ,|\t\n");
      j = (int)strtol(token, (char **)NULL, 10);
      token = strtok(NULL, " ,|\t\n");
      segdises[j] = strtod(token, (char **)NULL);
    }
    printf("Done reading in recombination fractions, plants used, and segregation distortions.\n");
  }
  else {
    printf("We are reading in the original data.\n");
    if (strcmp(inputformat, "mapmaker") != 0) {
      k = 0;
      while (fgets(line, maxprdlinelength, PRD)) {
        token = strtok(line, " ,|\t\n");
        if (token[0] == '/' && token[1] == '/') continue;
        preDosage[k] = (int)strtol(token, (char **)NULL, 10);
        k++;
        while ((token = strtok(NULL, " ,|\t\n"))) {
          if (token[0] == '/' && token[1] == '/') break;
          if (k >= nLoci) {
            fprintf(OUT, "This attempt will exceed the bounds of the preDosage array, k = %d\n", k);
            return(2);
          }
          preDosage[k] = (int)strtol(token, (char **)NULL, 10);
          k++;
        }
      }
      fclose(PRD);
      if (k < nLoci) {
        fprintf(OUT, "There were only %d loci in the preDosage array.  We were expecting %d loci.\n", k, nLoci);
        return(1);
      }
    }
    if (strcmp(inputformat, "byplant") == 0 ) {
      printf("The input format is by plant.\n");
      if (strlen(markernamesfile) > 1) MNA = fopen(markernamesfile, "r");
      else {
        printf("The configuration file did not specify the marker names file.\n");
        return(1);
      }
      j = 0;
      while (fgets(line, maxprdlinelength, MNA)) {
        if (line[0] != '\n') {
          line[maxnamelength-1] = '\0';
          for (i = 0; i < maxnamelength; i++) {
            if (line[i] == '\n') line[i] = '\0';
          }
          strcpy(names[j], line);
          j++;
        }
      }
      fclose(MNA);
      ierror = 1;
      isfull = 0;
      while (fgets(line, maxdatlinelength, DAT)) {
        token = strtok(line, " ,|\t\n");
        if (token[0] == '/' && token[1] == '/') continue;
        if (token[0] == '>') {
          ierror = 0;
          if (isfull == 1) {
            for (i = 0; i < nLoci; i++) {
              for (j = i + 1; j < nLoci; j++) {
                if (postDosage[i] >= 0 && postDosage[j] >= 0) {
                  if (preDosage[i] - postDosage[i] != preDosage[j] - postDosage[j]) {
                    recombfrac[i][j]++;
                  }
                }
                else nMissing[i][j]++;
              }
            }
          }
          k = 0;
          isfull = 0;
          for (i = 0; i < nLoci; i++) postDosage[i] = -1;
        }
        else {
          if (ierror == 1) {
            fprintf(OUT, "Unexpected beginning of the data file.  A fasta-like > sign was expected.  Offending line is %s", line);
            return(1);
          }
          else {
            if (token[0] != '-') {
              postDosage[k] = (int)strtol(token, (char**)NULL, 10);
              plantsUsed[k]++;
            }
            k++;
            while ((token = strtok(NULL, " ,|\t\n"))) {
              if (token[0] == '/' && token[1] == '/') break;
              if (token[0] != '-') {
                postDosage[k] = (int)strtol(token, (char**)NULL, 10);
                plantsUsed[k]++;
              }
              k++;
              if (k == nLoci) isfull = 1;
            }
          }
        }
      }
      //Here we must process the last set of postDosages.
      if (isfull == 1) {
        for (i = 0; i < nLoci; i++) {
          for (j = i + 1; j < nLoci; j++) {
            if (postDosage[i] >= 0 && postDosage[j] >= 0) {
              if (preDosage[i] - postDosage[i] != preDosage[j] - postDosage[j]) recombfrac[i][j]++;
            }
            else nMissing[i][j]++;
          }
        }
      }
    }
    else if (strcmp(inputformat, "bygene") == 0) {
      printf("The input format is by gene.\n");
      //We fill orthodata orthogonally, by columns, so it is orthogonal to the data as input.
      //Once filled, orthodata can be read out just like postDosage to calculate recombination fractions.
      orthodata = (int **) calloc(plantsTotal, sizeof(int *));
      for (i = 0; i < plantsTotal; i++) orthodata[i] = (int *) calloc(nLoci, sizeof(int));
      tempdata = (int *) calloc(plantsTotal, sizeof(int));
      ierror = 1;
      isfull = 0;
      ig = 0; //for loci
      while (fgets(line, maxdatlinelength, DAT)) {
        token = strtok(line, " ,|\t\n");
        if (token[0] == '/' && token[1] == '/') continue;
        if (token[0] == '>') {
          if (isfull == 1) {
            //tempdata contains values by plant to be loaded into orthodata.
            for (i = 0; i < plantsTotal; i++) orthodata[i][ig] = tempdata[i];
            ig++;
          }
          ierror = 0;
          for (i = 0; i < plantsTotal; i++) tempdata[i] = -1;
          ip = 0;
          memset(tempname, '\0', maxnamelength);
          iv = 1;
          while (token[iv] != '\n' && token[iv] != '\0') {
            tempname[iv-1] = token[iv];
            iv++;
          }
          strcpy(names[ig], tempname);
        }
        else {
          if (ierror == 1) {
            fprintf(OUT, "Unexpected beginning of the data file.  A fasta-like > sign was expected.  Offending line is %s", line);
            return(1);
          }
          else {
            if (token[0] != '-') {
              tempdata[ip] = (int)strtol(token, (char**)NULL, 10);
              plantsUsed[ig]++;
            }
            ip++;
            while ((token = strtok(NULL, " ,|\t\n"))) {
              if (token[0] == '/' && token[1] == '/') break;
              if (token[0] != '-') {
                tempdata[ip] = (int)strtol(token, (char**)NULL, 10);
                plantsUsed[ig]++;
              }
              ip++;
              if (ip == plantsTotal) isfull = 1;
            }
          }
        }
      }
      //Now process the last line.
      if (isfull == 1) {
        for (i = 0; i < plantsTotal; i++) orthodata[i][ig] = tempdata[i];
        ig++;
        if (ig != nLoci) {
          fprintf(OUT, "Number of read genes %d under bygene option does not agree with the stated total, %d.\n", ig, nLoci);
          return(1);
        }
      }
      printf("We are done reading data file %s.\n", datafile);
      //Now feed orthodata to recombfrac.
      for (i = 0; i < plantsTotal; i++) {
        printf("We are processing plant %d from orthodata.\n", i);
        for (j = 0; j < nLoci; j++) {
          for (k = j + 1; k < nLoci; k++) {
            if (orthodata[i][j] >= 0 && orthodata[i][k] >= 0) {
              if (preDosage[j] - orthodata[i][j] != preDosage[k] - orthodata[i][k]) recombfrac[j][k]++;
            }
            else nMissing[j][k]++;
          }
        }
      }
      printf("We are done processing orthodata.\n");
      //free(orthodata);
    }
    else if (strcmp(inputformat, "mapmaker") == 0) {
      printf("The input format is mapmaker.\n");
      for (i = 0; i < nLoci; i++) preDosage[i] = 1; //Assume MapMaker input to be diploid, with only presence or absence scored.
      strcpy(mmdatatype, "NULL");
      fprintf(OUT, "We have a MapMaker input file.\n");
      fgets(line, maxdatlinelength, DAT);
      token = strtok(line, " \t\n");
      i = 0;
      while (token[i] != '\0') {
        lctoken[i] = tolower(token[i]);
        i++;
      }
      lctoken[i] = '\0';
      if (strcmp(lctoken, "data") == 0) {
        token = strtok(NULL, " \t\n");
        i = 0;
        while (token[i] != '\0') {
          lctoken[i] = tolower(token[i]);
          i++;
        }
        lctoken[i] = '\0';
        if (strcmp(lctoken, "type") == 0) {
          token = strtok(NULL, " \t\n");
          i = 0;
          while (token[i] != '\0') {
            lctoken[i] = tolower(token[i]);
            i++;
          }
          lctoken[i] = '\0';
          if (strcmp(lctoken, "ri") == 0 || strcmp(lctoken, "ril") == 0 || strcmp(lctoken, "f2") == 0 || strcmp(lctoken, "f3") == 0) {
            token = strtok(NULL, " \t\n");  
            i = 0;
            while (token[i] != '\0') {
              lctoken[i] = tolower(token[i]);
              i++;
            }
            lctoken[i] = '\0';
            if (strcmp(lctoken, "backcross") == 0 || strcmp(lctoken, "intercross") == 0 || strcmp(lctoken, "self") == 0 || strcmp(lctoken, "sib") == 0) {
              fprintf(OUT, "The MapMaker data type is supported; Flipper2 will continue.\n");
              strcpy(mmdatatype, lctoken);
            }
            else {
              fprintf(OUT, "Something is wrong or unsupported with the first line of the MapMaker file.  Flipper is halting.\n");
              fprintf(OUT, "The offending line is %s with token %s\n", line, token);
              return(1);
            }
          }
          else {
            fprintf(OUT, "Something is wrong or unsupported with the first line of the MapMaker file.  Flipper is halting.\n");
            fprintf(OUT, "The offending line is %s with token %s\n", line, token);
            return(1);
          }
        }
        else {
          fprintf(OUT, "Something is wrong or unsupported with the first line of the MapMaker file.  Flipper is halting.\n");
          fprintf(OUT, "The offending line is %s with token %s\n", line, token);
          return(1);
        }
      }
      else {
        fprintf(OUT, "Something is wrong or unsupported with the first line of the MapMaker file.  Flipper is halting.\n");
        fprintf(OUT, "The offending line is %s with token %s\n", line, token);
        return(1);
      }
      fprintf(OUT, "The MapMaker data type is %s\n", mmdatatype);
      fgets(line, maxdatlinelength, DAT);
      token = strtok(line, " \t\n");
      claimedplants = (int)strtol(token, (char**)NULL, 10);
      if (claimedplants != plantsTotal) {
        fprintf(OUT, "Claimed number %d of mapping individuals in MapMaker data does not agree with total %d in the configuration file.\n", claimedplants, plantsTotal);
        return(1);
      }
      token = strtok(NULL, " \t\n");
      claimedloci = (int)strtol(token, (char**)NULL, 10);
      fprintf(OUT, "The MapMaker file claims %d loci, versus %d in the configuration file.\n", claimedloci, nLoci);
      if (combineparents == 1) {
        //Read MapMaker data to produce the haploid number of linkage groups.
        chdata = (char **) calloc(nLoci, sizeof(char *));
        for (i = 0; i < nLoci; i++) chdata[i] = (char *) calloc(plantsTotal, sizeof(char));
        ig = 0;
        while(fgets(line, maxdatlinelength, DAT)) {
          if (line[0] == '\n' || line[0] == '\0' || line[0] == '#') continue; //Skip comments or blank lines.
          if (line[0] == '/' && (line[1] == '/' || line[1] == '\n')) continue;
          if (line[0] == '*') {
            for (ivv = 1; ivv < maxnamelength; ivv++) {
              if (line[ivv] == '\n' || line[ivv] == '\0' || line[ivv] == '#') break;
              if (line[ivv] == '/' && (line[ivv+1] == '/' || line[ivv+1] == '\n')) break;
              tempname[ivv-1] = line[ivv];
            }
            if (strlen(tempname) > 0) strcpy(names[ig], tempname);
            continue;
          }
          if (line[0] == 'A' || line[0] == 'H' || line[0] == 'B' || line[0] == '-') {
            ik = 0;
            for (iv = 0; iv < maxdatlinelength; iv++) {
              if (line[iv] != ' ' && line[iv] != '\t') {
                chdata[ig][ik] = line[iv];
                ik++;
              }
            }
            ig++;
            continue;
          }
          fprintf(OUT, "Bad data line = %s\n", line);
          return(-1);
        }
        fprintf(OUT, "File %s gave %d loci\n", datafile, ig);
        //fprintf(OUT, "%s gave %d loci\nThe data are:\n", datafile, ig);
        //for (i = 0; i < nLoci; i++) {
        //  fprintf(OUT, "%d\t%s\t", i, names[i]);
        //  for (iv = 0; iv < plantsTotal; iv++) fprintf(OUT, "%c", chdata[i][iv]);
        //  fprintf(OUT, "\n");
        //}
        for (i = 0; i < nLoci; i++) {
          for (j = i + 1; j < nLoci; j++) {
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
        for (i = 0; i < nLoci; i++) {
          As = 0;
          Bs = 0;
          highAs = 0;
          highBs = 0;
          lowAs = 0;
          lowBs = 0;
          for (j = 0; j < plantsTotal; j++) {
            if (chdata[i][j] != '-') plantsUsed[i]++;
            if (chdata[i][j] == 'A') {
              As += 2;
            }
            if (chdata[i][j] == 'B') {
              Bs += 2;
            }
            if (chdata[i][j] == 'H') {
              As++;
              Bs++;
            }
          }
          if (strcmp(mmdatatype, "backcross") == 0) {
            if (As > Bs) segdises[i] = 0.5 + 2 * abs(((double) As / (As + Bs)) - 0.75);
            else segdises[i] = 0.5 + 2 * abs(((double) Bs / (As + Bs)) - 0.75);
          }
          else { //This measure does not account for highly heterozygous populations, i.e., with balanced lethal recessives.
            if (As > Bs) segdises[i] = (double) As / (As + Bs);
            else segdises[i] = (double) Bs / (As + Bs);
          }
        }
      }
      else { //Here combineparents == 0.
        //From here on, MapMaker files are processed much like bygene files.  The tempdata array is filled from the file, then fed into orthodata.
        //Read MapMaker data to map the parents separately, totalling the somatic number of linkage groups.
        orthodata = (int **) calloc(plantsTotal, sizeof(int *));
        for (i = 0; i < plantsTotal; i++) orthodata[i] = (int *) calloc(nLoci, sizeof(int));
        tempdata = (int *) calloc(plantsTotal, sizeof(int));
        sectempdata = (int *) calloc(plantsTotal, sizeof(int));
        for (i = 0; i < plantsTotal; i++) {
          tempdata[i] = -1;
          sectempdata[i] = -1;
        }
        for (i = 0; i < plantsTotal; i++) {
          for (j = 0; j < nLoci; j++) orthodata[i][j] = -1;
        }
        isfull = 0;
        if (strcmp(mmdatatype, "backcross") == 0) ig = -1; //for loci
        else ig = -2;
        ip = 0;
        while(fgets(line, maxdatlinelength, DAT)) {
          nameflag = 0;
          for (iv = 0; iv < maxdatlinelength; iv++) {
            if (line[iv] == '\n' || line[iv] == '\0' || line[iv] == '#') break;
            if (line[iv] == '/' && (line[iv+1] == '/' || line[iv+1] == '\n')) break;
            if (line[iv] == '*') {
              if (isfull == 1) {
                //tempdata contains values by plant to be loaded into orthodata.
                if (strcmp(mmdatatype, "backcross") == 0) { 
                  for (i = 0; i < plantsTotal; i++) orthodata[i][ig] = tempdata[i];
                  //Handle segregation distortion.
                  if (As > Bs) { //H now contributes to As and Bs, so their sum is always greater than zero.
                    segdises[ig] = (double) As / (As + Bs);
                  }
                  else {
                    segdises[ig] = (double) Bs / (As + Bs);
                  }
                }
                else {
                  for (i = 0; i < plantsTotal; i++) {
                    orthodata[i][ig] = tempdata[i];
                    orthodata[i][ig+1] = sectempdata[i];
                  }
                  if (As > Bs) { //H now contributes to As and Bs, so their sum is always greater than zero.
                    segdises[ig] = (double) As / (As + Bs);
                    segdises[ig+1] = (double) As / (As + Bs);
                  }
                  else {
                    segdises[ig] = (double) Bs / (As + Bs);
                    segdises[ig+1] = (double) Bs / (As + Bs);
                  }
                }
              }
              else if (ip > 0) {
                fprintf(OUT, "Next MapMaker gene was encountered before tempdata was full, ip = %d\n", ip);
                return(2);
              }
              if (strcmp(mmdatatype, "backcross") == 0) {
                ig++;
                plantsUsed[ig] = 0;
              }
              else {
                ig += 2;
                plantsUsed[ig] = 0;
                plantsUsed[ig+1] = 0;
              }
              nameflag = 1;
              j = 0;
              memset(tempname, '\0', maxnamelength);
              for (ivv = iv + 1; ivv < maxnamelength; ivv++) {
                if (line[ivv] == '\n' || line[ivv] == '\0' || line[ivv] == '#') break;
                if (line[ivv] == '/' && (line[ivv+1] == '/' || line[ivv+1] == '\n')) break;
                tempname[ivv-iv-1] = line[ivv];
              }
              if (strcmp(mmdatatype, "backcross") == 0) {
                if (strlen(tempname) > 0) strcpy(names[ig], tempname);
              }
              else {
                if (strlen(tempname) > 0) {
                  strcpy(names[ig], tempname);
                  strcpy(names[ig+1], tempname);
                  strncat(names[ig], "-1", 2);
                  strncat(names[ig+1], "-2", 2);
                }
              }
              for (i = 0; i < plantsTotal; i++) {
                tempdata[i] = -1;
                sectempdata[i] = -1;
              }
              ip = 0;
              isfull = 0;
              As = 0;
              Bs = 0;
              highAs = 0;
              highBs = 0;
              lowAs = 0;
              lowBs = 0;
              break;
            }
            if (nameflag == 0) {
              if (strcmp(mmdatatype, "backcross") == 0) {
                if (line[iv] == 'h' || line[iv] == 'H') {
                  tempdata[ip] = 1;
                  plantsUsed[ig]++;
                  As++;
                  Bs++;
                  ip++;
                }
                if (line[iv] == 'a' || line[iv] == 'A') {
                  tempdata[ip] = 0;
                  plantsUsed[ig]++;
                  As += 2;
                  ip++;
                }
                if (line[iv] == 'b' || line[iv] == 'B') {
                  tempdata[ip] = 0;
                  plantsUsed[ig]++;
                  Bs += 2;
                  ip++;
                }
                if (line[iv] == '-') ip++;
                if (line[iv] == 'c' || line[iv] == 'C' || line[iv] == 'd' || line[iv] == 'D') {
                  fprintf(OUT, "Flipper2 currently cannot handle MapMaker backcross files scored for dominant markers.\n");
                  fprintf(OUT, "Please reformat your data in Flipper2\'s native bygene format.\n");
                  return(2);
                }
                //Else do nothing. For example, do nothing with a blank space before or after the gene's name.
                if (ip == plantsTotal) isfull = 1;
              }
              else {
                //All inbred populations go here.
                if (line[iv] == 'h' || line[iv] == 'H') {
                  tempdata[ip] = 1;
                  sectempdata[ip] = 1;
                  plantsUsed[ig]++;
                  plantsUsed[ig+1]++;
                  As++;
                  Bs++;
                  ip++;
                }
                if (line[iv] == 'a' || line[iv] == 'A') {
                  tempdata[ip] = 1;
                  sectempdata[ip] = 0;
                  plantsUsed[ig]++;
                  plantsUsed[ig+1]++;
                  As += 2;
                  ip++;
                }
                if (line[iv] == 'b' || line[iv] == 'B') {
                  tempdata[ip] = 0;
                  sectempdata[ip] = 1;
                  plantsUsed[ig]++;
                  plantsUsed[ig+1]++;
                  Bs += 2;
                  ip++;
                }
                if (line[iv] == '-') {
                  tempdata[ip] = -1;
                  sectempdata[ip] = -1;
                  ip++;
                }
                if (line[iv] == 'c' || line[iv] == 'C' || line[iv] == 'd' || line[iv] == 'D') {
                  fprintf(OUT, "Flipper2 currently cannot handle MapMaker inbred files scored for dominant markers.\n");
                  fprintf(OUT, "Please reformat your data in Flipper2\'s native bygene format.\n");
                  return(2);
                }
                //Else do nothing.  
                if (ip == plantsTotal) isfull = 1;
              }
            }
          }
        }
        if (isfull == 1) { //Process the last line here.
          //fprintf(OUT, "ig = %d As = %d Bs = %d tempdata[2] = %d in the last line of the mapmaker input file.\n", ig, As, Bs, tempdata[2]);
          if (strcmp(mmdatatype, "backcross") == 0) {
            if (As > Bs) segdises[i] = 0.5 + 2 * abs(((double) As / (As + Bs)) - 0.75);
            else segdises[i] = 0.5 + 2 * abs(((double) Bs / (As + Bs)) - 0.75);
            for (i = 0; i < plantsTotal; i++) {
              orthodata[i][ig] = tempdata[i];
              tempdata[i] = -1;
            }
            ig++;
          }
          else {
            //fprintf(OUT, "ig = %d to set the last two columns of orthodata.\n", ig);
            if (As > Bs) segdises[ig] = (double) As / (As + Bs);
            else segdises[ig] = (double) Bs / (As + Bs);
            segdises[ig+1] = segdises[ig];
            for (i = 0; i < plantsTotal; i++) {
              orthodata[i][ig] = tempdata[i];
              tempdata[i] = -1;
              orthodata[i][ig+1] = sectempdata[i];
              sectempdata[i] = -1;
            }
            ig += 2;
          }
          if (ig != nLoci) {
            fprintf(OUT, "Number of read MapMaker genes %d does not agree with the stated total, %d.\n", ig, nLoci);
            return(1);
          }
        }
        //Now feed orthodata to recombfrac.
        for (i = 0; i < plantsTotal; i++) {
          for (j = 0; j < nLoci; j++) {
            for (k = j + 1; k < nLoci; k++) {
              if (orthodata[i][j] >= 0 && orthodata[i][k] >= 0) {
                if (preDosage[j] - orthodata[i][j] != preDosage[k] - orthodata[i][k]) recombfrac[j][k]++;
              }
              else nMissing[j][k]++;
            }
          }
        }
        //free(orthodata);
        //for (ka = 0; ka < nLoci; ka++) fprintf(OUT, "names %d\t%s\n", ka, names[ka]);
      } //End of the section to handle orthodata.
    }
    for (i = 0; i < nLoci; i++) {
      for (j = 0; j < i; j++) {
        nMissing[i][j] = nMissing[j][i];
      }
    }
    for (i = 0; i < nLoci; i++) {
      for (j = 0; j < nLoci; j++) {
        if (nMissing[i][j] != nMissing[j][i]) {
          fprintf(OUT, "Reciprocal nMissings differ, miss[%d][%d] = %d, miss[%d][%d] = %d\n", i, j, nMissing[i][j], j, i, nMissing[j][i]);
        }
      }
    }
    k = 0;
    //The following loop calculates recombination fractions on the basis of the plants for which recombination of i and j is available.
    fprintf(OUT, "We are about to calculate recombination fractions.\n");
    for (i = 0; i < nLoci; i++) {
      for (j = i + 1; j < nLoci; j++) {
        recombfrac[i][j] /= (plantsTotal - nMissing[i][j]);
        recombfrac[j][i] = recombfrac[i][j]; //Finally fill in the other half of the recombfrac array.
      }
    }
    for (i = 0; i < nLoci; i++) {
      for (j = 0; j < nLoci; j++) {
        if (recombfrac[i][j] != recombfrac[j][i]) fprintf(OUT, "reciprocals are unequal, ij = %.4f, ji = %.4f\n", recombfrac[i][j], recombfrac[j][i]);
      }
    }
  }
  imax = 0;
  for (i = 0; i < nLoci; i++) {
    frac = (double) plantsUsed[i] / plantsTotal;
    if (frac >= minFracFirst) fracUsedCode[i] = 2;
    else if (frac >= minFracSecond) fracUsedCode[i] = 1;
    else fracUsedCode[i] = 0; //for loci with too many missing data
  }
  timeused2 = clock();
  duration = (double) (timeused2 - timeused1) /CLOCKS_PER_SEC;
  fprintf(OUT, "End of data accumulation, time used = %f\n", duration);
  fflush(OUT);
  for (i = 0; i < nLoci; i++) {
    minofrfs[i] = 1;
    for (j = 0; j < nLoci; j++) {
      if (i == j) continue;
      if (recombfrac[i][j] < minofrfs[i]) minofrfs[i] = recombfrac[i][j];
    }
  }
  j = 0;
  fprintf(OUT, "List of loci excluded by minimum recombination fraction:\n");
  for (i = 0; i < nLoci; i++) {
    if (minofrfs[i] >= maxminrf) {
      fprintf(OUT, "%d\t%.4f\t%s\n", i, minofrfs[i], names[i]);
      j++;
    }
  }
  if (j == 0) fprintf(OUT, "NONE\n");
  fprintf(OUT, "\n");
  j = 0;
  fprintf(OUT, "List of loci excluded by segregation distortion with allelic ratio:\n");
  for (i = 0; i < nLoci; i++) {
    if (segdises[i] > segdisthreshold) {
      fprintf(OUT, "%d\t%.4f\t%s\n", i, segdises[i], names[i]);
      j++;
    }
  }
  if (j == 0) fprintf(OUT, "NONE\n");
  fprintf(OUT, "\n");
  j = 0;
  fprintf(OUT, "List of loci excluded by too many missing data:\n");
  for (i = 0; i < nLoci; i++) {
    if (fracUsedCode[i] <= 1) {
      fprintf(OUT, "%d\t%s\n", i, names[i]);
      j++;
    }
  }
  if (j == 0) fprintf(OUT, "NONE");
  fprintf(OUT, "\n");
  for (i = 0; i < nLoci; i++) {
    foreact[i] = -1;
    aftact[i] = -1;
    foremap[i] = -1;
    aftmap[i] = -1;
  }
  INM = fopen(inmapfile, "r");
  //fgets(line, maxdatlinelength, INM); //This is the header line of general information.
  //fprintf(OUT, "First line of INM is: %s\n", line);
  lastmaplocus = -1;
  lastactlocus = -1;
  while (fgets(line, maxdatlinelength, INM)) {
    if (line[0] == 'F' || line[0] == 'f') {
      if (line[1] == 'o' && line[2] == 'r' && line[3] == ' ') { //For mapped group 0
        if (lastmaplocus >= 0) fprintf(OUT, "linkage group = %d mapped length = %f physical length = %f\n", lg, mappedlength, physlength);
        lastmaplocus = -1;
        lastactlocus = -1;
        mappedlength = 0.0;
        physlength = 0.0; 
        token = strtok(line, " ,|\t\n");
        token = strtok(NULL, " ,|\t\n");
        token = strtok(NULL, " ,|\t\n");
        token = strtok(NULL, " ,|\t\n");
        lg = (int)strtol(token, (char**)NULL, 10);
      }
    }
    else {
      token = strtok(line, " ,|\t\n");
      currmaplocus = (int)strtol(token, (char**)NULL, 10);
      token = strtok(NULL, " ,|\t\n");
      curractlocus = (int)strtol(token, (char**)NULL, 10);
      fprintf(OUT, "diag: currmaplocus = %d curractlocus = %d\n", currmaplocus, curractlocus);
      if (currmaplocus == -1 && curractlocus == -1) fprintf(OUT, "mapped length = %.6f physical length = %.6f\n", mappedlength, physlength);
      else {
        if (lastmaplocus >= 0) mappedlength += recombfrac[lastmaplocus][currmaplocus];
        if (lastactlocus >= 0) physlength += recombfrac[lastactlocus][curractlocus];
        fprintf(OUT, "dial: mappedlength = %f physlength = %f\n", mappedlength, physlength);
      }
      lastmaplocus = currmaplocus;
      lastactlocus = curractlocus;
    }
  }
  if (lastmaplocus >= 0) fprintf(OUT, "linkage group = %d mapped length = %f physical length = %f\n", lg, mappedlength, physlength);
}

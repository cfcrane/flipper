#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <cmath>

void printMap(int, int, int *, int *, double *, double **, FILE *);
void printStats(char *, int, double **, int *, int *, FILE *);
void countEnds(int, int *, int *, int *, int *, FILE *);
int checkForPresence(int, int, int, int*, int*, FILE *);
void printLGLengths(int, int, int *, int *, int *, FILE *);
int checkForCircularLG(int, int *, int *);
int dcomp(const void *, const void *);
int comp (const void *, const void *);
int compdown (const void *, const void *);
void printLinkageGroup(int, int *, int *, double **, FILE *);
void joinpenult(int, int, int *, int *, double **, FILE *);
int getmaxrf(int *, double **, int, double, double, FILE *);
void printtestloci(int, int, int, double **, int *, int *, FILE *);

int main(int argc, char *argv[]) {
  /*This version has a new module to move user-specified markers to their nearest neighbors in the "forced" map.
  New variables include char *flagfile, char **flaggednames, char *flaggedname, int *flagged, int flaglines, and int tryflag.  */
  char *datafile, *prdfile, *outfile, *specialfile, *ngrfile, *outafile, *measfile, *varname, *value, inputformat[512];
  char *line, *token, *lctoken, *queryfile, vars[5][512], **names, *cfgline, *printstring, *flagfile, **flaggednames, *flaggedname;
  char *tempname, mmdatatype[20], *markernamesfile, *rfracfile, *filler, **chdata, *usedcodefile, *segdisfile, *dedupoutfile;
  int savedj, imax, kmax, queryflag, queryarray[20], *querycodes, mapqueries[20], mq, *foretestins, *afttestins, *flaggedindices;
  int i, j, k, l, nLoci, plantsTotal, baseNumber, ploidy, ierror, ***coeffs, mergearray[10], *foretestatt, *afttestatt, nexemplars;
  int *forePointers, *aftPointers, mappedLoci, nLinkageGroups, lowindex, endstatusk, endstatusm, *exemplars, dedupoutflag;
  int g, h, m, n, km, nfac, mini, savedv[20], maxnamelength, iv, ivv, outaflag, *membership, highindex, lowcount, highcount;
  int *fracUsedCode, isfull, maxprdlinelength, maxdatlinelength, *preDosage, *postDosage, **nMissing, *plantsUsed, *donenewLinkageGroups;
  int **orthodata, *tempdata, *sectempdata, found, temp, found1, found2, found3, measflag, *highlowstatus, lowBs, savedi;
  int *lgStarts, *lgEnds, *lgCounts, *sortedlgCounts, *doneLinkageGroups, nf[12], **actmap, totmeas, highAs, highBs, lowAs;
  int ia, ib, ic, id, ie, iF, ig, ih, ii, ij, ip, stemlocusa, stemlocusb, permuteon, ik, im, *mafts, *mfores, ir, kr, insrounds;
  int rfileflag, specialfileflag, nameflag, claimedloci, claimedplants, testlocusa, testlocusb, As, Bs, missingmeas, minflag;
  int flankfore, flankaft, anneallength, foreendm, aftendm, *insertions, savedfore, savedaft, instodo, pused, savedg, excisionflag;
  int gycount, gzcount, commit, rfxflag, combineparents, cllength, mincount, continuation, *lgnewCounts, *lgnewStarts, nexcised;
  int *pool, *lgnewEnds, *sortedlgnewCounts, nfloop, totcount, currcount, printflag, dedupflag, exclusionflag, flaglines, orignlgroups;
  int savedjf, savedja, foredist, aftdist, *minlgothergroup, starta, startb, *nearestneighbors, runtests, tryflag, *flagged;
  double segdisthreshold, *segdises, actlength, maxminrf, *minofrfs, **minlgdists, lowcutoff, highcutoff, *hightolows;
  double **recombfrac, mapLength, bandwidth, base, top, maxrf, queriedmaplength, duration, cudist, tdist, tooclose;
  double frac, minFracFirst, minFracSecond, *lgLengths, minmaplength, localmaplength, attlength, inslength, trial;
  double rfxatt, rfxins, rfxceiling, tripletlength, delta, rfthreshold, mindist, trialdista, trialdistb, *lgnewLengths;
  double secrfthreshold, temprfi, temprfsj, totlength, lengthpercount, *measurements, *sortedmeasurements, meas;
  double *minlgdistany;
  FILE *CFG, *OUT, *DAT, *PRD, *NGR, *SPE, *QRY, *ACT, *MNA, *RFF, *PUD, *SGD, *DDP, *MEA, *FLG;
  if(strcmp(argv[1], "-") == 0 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-help") == 0) {
    printf("Usage: %s configuration_file_name [configuration_file_line_length]\n", argv[0]);
    printf("The configuration file name is mandatory and the configuration file line length is optional.\n");
    printf("Specify the configuration file line length only if it needs to exceed 512 characters.\n");
    printf("Comments in the configuration file begin with # or // and can follow data if separated by one or more spaces.\n");
    printf("Data lines can appear in any order.  Filenames need not be quoted and can include the full path.\n");
    printf("Data lines consist of space-delimited triplets of attribute, \"=\", and value.\n");
    printf("Allowed attributes include datafile, specialfile, prdfile, outfile, ngrfile, outafile,\n");
    printf(" queryfile, markernamesfile, rfracfile, usedcodefile, segdisfile, inputformat, nLoci,\n");
    printf(" plantsTotal, maxnamelength, testlocusa, testlocusb, minFracFirst, minFracSecond, ploidy,\n");
    printf(" baseNumber, maxrf, maxminrf, rfxflag, rfxceiling, segdisthreshold, maxprdlinelength,\n");
    printf(" maxdatlinelength, permuteon, continuation, rfthreshold, dedupflag, dedupoutflag, measfile,\n");
    printf(" combineparents, and flagfile.\n");
    printf("Attributes rfracfile and specialfile are mutually exclusive, since they control if the\n");
    printf(" original data are read or a file of recombination fractions is read. Generally both\n");
    printf(" rfracfile and special file will designate the same file name.\n");
    printf("Attribute continuation specifies which of three maps is output: the full map, the map minus ill-fitting and probably\n");
    printf(" erroneous loci, or the map with ill-fitting loci forced to their optimal positions in the main linkage groups.\n");
    printf("There is the option to run the final adjustment loop more than once by setting nfloop to 2 or higher, but this\n");
    printf(" should be done only if continuation is set to 2.  If nfloop is set, the value of secrfthreshold should be set also.\n");
    return(0);
  }
  printstring = (char *)calloc(200, sizeof(char));
  clock_t timeused, timeused1, timeused2, currtime;
  CFG = fopen(argv[1], "r");
  timeused = clock();
  queryflag = 0; //This controls any action on an input file with a queried marker order.
  outaflag = 0; //This controls output of the length of the actual map order in simulations.
  dedupoutflag = 0; //This controls whether to output a deduplicated data file.
  dedupflag = 0; //This controls whether to map a deduplicated set of markers.
  exclusionflag = 0; //Setting this to 1 will permanently exclude loci that do not recombine with exemplars.
  rfileflag = 0; //This controls whether to read in recombination fractions directly.
  specialfileflag = 0; //This controls whether to write recombination fractions to SPE.
  measflag = 0; //This is 1 if a file of quantitative measurements is used to produce high-low ratios.
  excisionflag = 0; //This controls how to excise misfit markers from the initial map; zero is simple excision if fore and aft rf values are too high.
  combineparents = 0; //This controls whether to map all parental chromosomes separately (0) or together to make baseNumber linkage groups (1).
  maxprdlinelength = 0; //Otherwise this is not getting set with MapMaker input.
  tryflag = 0; //This controls whether to try to move the user-specified markers listed in flagfile.
  insrounds = 1; //This is the number of times the insertion loop runs that puts removed or duplicate markers on the map.
  delta = -1;
  runtests = 0; //Setting this to 1 will output various diagnostic results.
  minflag = 0; //This controls insertion of misfit loci, either with minimal increase in map length (0) or next to nearest neighbor (1).
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
  specialfile = (char *) calloc(cllength + 2, sizeof(char));
  prdfile = (char *) calloc(cllength + 2, sizeof(char));
  outfile = (char *) calloc(cllength + 2, sizeof(char));
  ngrfile = (char *) calloc(cllength + 2, sizeof(char));
  outafile = (char *) calloc(cllength + 2, sizeof(char));
  queryfile = (char *) calloc(cllength + 2, sizeof(char));
  measfile = (char *) calloc(cllength + 2, sizeof(char));
  markernamesfile = (char *) calloc(cllength + 2, sizeof(char));
  rfracfile = (char *) calloc(cllength + 2, sizeof(char));
  usedcodefile = (char *) calloc(cllength + 2, sizeof(char));
  segdisfile = (char *) calloc(cllength + 2, sizeof(char));
  dedupoutfile = (char *) calloc(cllength + 2, sizeof(char)); //Reserve space even if it is not ever filled or used.
  flagfile = (char *) calloc(cllength + 2, sizeof(char)); //Reserve space even if it is not ever filled or used.
  while (fgets(cfgline, cllength, CFG)) {
    fflush(NULL);
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
    if (strcmp(varname, "dedupoutfile") == 0) {
      strcpy(dedupoutfile, value);
      dedupoutflag = 1;
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
    if (strcmp(varname, "specialfile") == 0) {
      strcpy(specialfile, value);
      specialfileflag = 1;
      ierror = 0;
    }
    if (strcmp(varname, "flagfile") == 0) {
      strcpy(flagfile, value);
      tryflag = 1;
      ierror = 0;
    }
    if (strcmp(varname, "usedcodefile") == 0) {
      strcpy(usedcodefile, value);
      ierror = 0;
    }
    if (strcmp(varname, "segdisfile") == 0) {
      strcpy(segdisfile, value);
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
    if (strcmp(varname, "queryfile") == 0) {
      strcpy(queryfile, value);
      queryflag = 1;
      ierror = 0;
    }
    if (strcmp(varname, "measfile") == 0) {
      strcpy(measfile, value);
      measflag = 1;
      ierror = 0;
    }
    if (strcmp(varname, "kmax") == 0) {
      kmax = (int)strtol(value, (char **)NULL, 10);
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
    if (strcmp(varname, "testlocusa") == 0) {
      testlocusa = (int)strtol(value, (char**)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "testlocusb") == 0) {
      testlocusb = (int)strtol(value, (char**)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "runtests") == 0) {
      runtests = (int)strtol(value, (char**)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "insrounds") == 0) {
      insrounds = (int)strtol(value, (char**)NULL, 10);
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
    if (strcmp(varname, "lowcutoff") == 0) {
      lowcutoff = strtod(value, (char **)NULL);
      ierror = 0;
    }
    if (strcmp(varname, "highcutoff") == 0) {
      highcutoff = strtod(value, (char **)NULL);
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
    if (strcmp(varname, "maxrf") == 0) {
      maxrf = strtod(value, (char **)NULL);
      ierror = 0;
    }
    if (strcmp(varname, "maxminrf") == 0) {
      maxminrf = strtod(value, (char **)NULL);
      ierror = 0;
    }
    if (strcmp(varname, "dedupflag") == 0) {
      dedupflag = (int)strtol(value, (char **)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "exclusionflag") == 0) {
      exclusionflag = (int)strtol(value, (char **)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "excisionflag") == 0) {
      excisionflag = (int)strtol(value, (char**)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "minflag") == 0) {
      minflag = (int)strtol(value, (char**)NULL, 10);
      if (minflag != 0 && minflag != 1) ierror = 1;
      else ierror = 0;
    }
    if (strcmp(varname, "rfxflag") == 0) {
      rfxflag = (int)strtol(value, (char **)NULL, 10);
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
    if (strcmp(varname, "permuteon") == 0) {
      permuteon = (int)strtol(value, (char**)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "combineparents") == 0) {
      combineparents = (int)strtol(value, (char **)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "continuation") == 0) {
      continuation = (int)strtol(value, (char**)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "mincount") == 0) {
      mincount = (int)strtol(value, (char**)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "rfthreshold") == 0) {
      rfthreshold = strtod(value, (char **)NULL);
      ierror = 0;
    }
    if (strcmp(varname, "nfloop") == 0) {
      nfloop = (int)strtol(value, (char **)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "secrfthreshold") == 0) {
      secrfthreshold = strtod(value, (char **)NULL);
      ierror = 0;
    }
    if (strcmp(varname, "tooclose") == 0) {
      tooclose = strtod(value, (char **)NULL);
      ierror = 0;
    }
    if (ierror == 1) {
      printf("The configuration file has a misspelling, an incorrect value, or an extra variable.  Please fix it.\n");
      printf("The offending line contains: %s\n", varname);
      return(1);
    }
  }
  fclose(CFG);
  if (dedupoutflag > 0 && dedupflag == 0) {
    printf("To print out a deduplicated data file, dedupflag must also be set to 1 (to map deduplicated markers) or 2 (merely to output them).\n");
    return(1);
  }
  if (rfileflag > 0 && measflag > 0) {
    printf("Flipper currently calculates the measurement frequency ratio as it reads in the genotypes.\n");
    printf("Therefore, it cannot generate the measurement frequency ratio from tabulated recombination fractions.\n");
    printf("If you want to find QTL, set the rfileflag to zero.  Otherwise remove the measurement file line from the configuration file.\n");
    return(1);
  }
  if (strcmp(inputformat, "mapmaker") != 0 && measflag > 0) {
    printf("The QTL finder is currently set up only for MapMaker-formatted segregation data.\n");
    return(1);
  }
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
  if (specialfileflag == 1) {
    if (rfileflag == 1) {
      fprintf(OUT, "You are not allowed to read and write a file of recombination fractions in the same run.\n");
      fprintf(OUT, "Please remove the entry for at least one of these files from the configuration file.\n");
      fprintf(OUT, "The file names are %s for the specialfile and %s for the rfracfile.\n", specialfile, rfracfile);
      return(1);
    }
    if (strlen(specialfile) > 1) SPE = fopen(specialfile, "w");
    else {
      printf("The configuration file did not specify the special file to hold recombination fractions.\n");
      return(1);
    }
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
  if (measflag == 1 && highcutoff < lowcutoff) {
    printf ("The values of highcutoff and lowcutoff appear to have been reversed.  Check their values in the configuration file.\n");
    return(1);
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
  nearestneighbors = (int *) calloc(nLoci, sizeof(int));
  for (i = 0; i < nLoci; i++) nearestneighbors[i] = -1;
  minofrfs = (double *) calloc(nLoci, sizeof(double));
  segdises = (double *) calloc(nLoci, sizeof(double));
  insertions = (int *) malloc(nLoci * sizeof(int));
  mfores = (int *) malloc(nLoci * sizeof(int));
  mafts = (int *) malloc(nLoci * sizeof(int));
  membership = (int *) malloc(nLoci * sizeof(int));
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
  if (measflag == 1) {
     if (combineparents == 0) hightolows = (double *)calloc(2 * nLoci + 20, sizeof(double));
     if (combineparents == 1) hightolows = (double *)calloc(nLoci + 10, sizeof(double));
  }
  timeused1 = clock();
  duration = (double) (timeused1 - timeused) / CLOCKS_PER_SEC;;
  fprintf(OUT, "done with setup, time used = %f\n", duration);
  if (measflag == 1) {
    MEA = fopen(measfile, "r");
    if (MEA == NULL) {
      printf("The measurement file %s could not be found or opened, but the measurement flag is set to 1.\n", measfile);
      return(1);
    }
    measurements = (double*) calloc(plantsTotal + 50, sizeof(double));
    sortedmeasurements = (double*) calloc(plantsTotal + 50, sizeof(double));
    for (i = 0; i < plantsTotal; i++) { //Any missing names are treated as missing data.
      measurements[i] = -1; //code for missing data
      sortedmeasurements[i] = -1;
    }
    //Any header lines must begin with '#' or '//'.
    while (fgets(line, maxdatlinelength, MEA)) {
      if (line[0] == '\n' || line[0] == '\0' || line[0] == '#') continue; //Skip comments or blank lines.
      if (line[0] == '/' && (line[1] == '/' || line[1] == '\n')) continue;
      if (strlen(line) < 3) continue;
      token = strtok(line, " |\t\n,");
      if (token[0] == '/' && token[1] == '/') continue;
      if (token[0] == '#') continue;
      found = 0;
      for (i = 0; i < strlen(token); i++) {
        if (!isdigit(token[i]) && token[i] != '.') found++;
      }
      if (found > 0) fprintf(OUT, "Bad token is %s\n", token);
      if (found == 0) {
        k = (int)strtol(token, (char**)NULL, 10);
        if (k > plantsTotal) {
          printf("A plant number %d exceeds the total population size %d, which is not allowed.\n", k, plantsTotal);
          return(-1);
        }
        token = strtok(NULL, " |\t\n,");
        found = 0;
        for (i = 0; i < strlen(token); i++) {
          if (!isdigit(token[i]) && token[i] != '.') found++;
        }
        if (found == 0) { //Ignore the alternative; these are left by default as -1 and will be picked up as missing data.
          meas = strtod(token, (char **)NULL);
          measurements[k] = meas;
          sortedmeasurements[k] = meas;
        }
      }
    }
    qsort(sortedmeasurements, plantsTotal, sizeof(double), dcomp);
    totmeas = 0;
    missingmeas = 0;
    for (i = 0; i < plantsTotal; i++) {
      //fprintf(OUT, "meas: %.6f %.6f\n", measurements[i], sortedmeasurements[i]);
      if (measurements[i] < 0) missingmeas++;
      else totmeas++;
    }
    lowindex = int(missingmeas + lowcutoff * totmeas);
    highindex = int(missingmeas + highcutoff * totmeas);
    fprintf(OUT, "totmeas = %d missingmeas = %d lowindex = %d highindex = %d\n", totmeas, missingmeas, lowindex, highindex);
    fflush(OUT);
    highlowstatus = (int*) calloc(plantsTotal + 10, sizeof(int));
    lowcount = 0;
    highcount = 0;
    for (i = 0; i < plantsTotal; i++) {
      if (measurements[i] >= 0 && measurements[i] <= sortedmeasurements[lowindex]) {
        highlowstatus[i] = -1;
	lowcount++;
      }
      if (measurements[i] >= sortedmeasurements[highindex]) {
        highlowstatus[i] = 1;
        highcount++;
      }
    }
    //for (i = 0; i < plantsTotal; i++) fprintf(OUT, "Status of plant %d = %d\n", i, highlowstatus[i]);
  }
  if (rfileflag == 1) {
    if (specialfileflag == 1) {
      fprintf(OUT, "You are not allowed to read and write a file of recombination fractions in the same run.\n");
      fprintf(OUT, "Please remove the entry for at least one of these files from the configuration file.\n");
      fprintf(OUT, "The file names are %s for the specialfile and %s for the rfracfile.\n", specialfile, rfracfile);
      return(1);
    }
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
    if (strlen(usedcodefile) > 1) PUD = fopen(usedcodefile, "r");
    else {
      printf("The configuration file did not specify the name of usedcodefile, the file that holds the numbers of plants used.\n");
      return(1);
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
    if (strlen(segdisfile) > 1) SGD = fopen(segdisfile, "r");
    else {
      printf("The configuration file did not specify the name of segdisfile, the file that holds the segregation distortions.\n");
      return(1);
    }
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
              if (line[iv] == '\n' || line[iv] == '\0') break;
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
              if (measflag == 1) {
                if (highlowstatus[j] == -1) lowAs += 2;
                if (highlowstatus[j] == 1) highAs += 2;
              }
            }
            if (chdata[i][j] == 'B') {
              Bs += 2;
              if (measflag == 1) {
                if (highlowstatus[j] == -1) lowBs += 2;
                if (highlowstatus[j] == 1) highBs += 2;
              }
            }
            if (chdata[i][j] == 'H') {
              As++;
              Bs++;
              if (measflag == 1) {
                if (highlowstatus[j] == -1) {
                  lowAs++;
                  lowBs++;
                }
                if (highlowstatus[j] == 1) {
                  highAs++;
                  highBs++;
                }
              }
            }
          }
          if (measflag == 1) {
            if (highAs + highBs > 0 && lowAs > 0) hightolows[i] = (double)highAs * (lowAs + lowBs) / (lowAs * (highAs + highBs));
            else hightolows[i] = 2 / delta;
          }
          if (strcmp(mmdatatype, "backcross") == 0) {
            if (As > Bs) segdises[i] = 0.5 + 2 * std::abs(((double) As / (As + Bs)) - 0.75);
            else segdises[i] = 0.5 + 2 * std::abs(((double) Bs / (As + Bs)) - 0.75);
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
                  if (measflag == 1) { //This is for the QTL finder.
                    if (highAs + highBs > 0 && lowAs > 0) hightolows[ig] = (double)highAs * (lowAs + lowBs) / (lowAs * (highAs + highBs));
                    else hightolows[ig] = 2 / delta;
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
                  if (measflag == 1) {
                    if (highAs + highBs > 0 && lowAs > 0) hightolows[ig] = (double)highAs * (lowAs + lowBs) / (lowAs * (highAs + highBs));
                    else hightolows[ig] = 2 / delta;
                    if (highAs + highBs > 0 && lowBs > 0) hightolows[ig+1] = (double)highBs * (lowAs + lowBs) / (lowBs * (highAs + highBs));
                    else hightolows[ig+1] = 2 / delta;
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
                  if (measflag == 1) {
                    if (highlowstatus[ip] == 1) {
                      highAs++;
                      highBs++;
                    }
                    if (highlowstatus[ip] == -1) {
                      lowAs++;
                      lowBs++;
                    }
                  }
                  ip++;
                }
                if (line[iv] == 'a' || line[iv] == 'A') {
                  tempdata[ip] = 0;
                  plantsUsed[ig]++;
                  As += 2;
                  if (measflag == 1) {
                    if (highlowstatus[ip] == 1) highAs += 2;
                    if (highlowstatus[ip] == -1) lowAs += 2;
                  }
                  ip++;
                }
                if (line[iv] == 'b' || line[iv] == 'B') {
                  tempdata[ip] = 0;
                  plantsUsed[ig]++;
                  Bs += 2;
                  if (measflag == 1) {
                    if (highlowstatus[ip] == 1) highBs += 2;
                    if (highlowstatus[ip] == -1) lowBs += 2;
                  }
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
                  if (measflag == 1) {
                    if (highlowstatus[ip] == 1) {
                      highAs++;
                      highBs++;
                    }
                    if (highlowstatus[ip] == -1) {
                      lowAs++;
                      lowBs++;
                    }
                  }
                  ip++;
                }
                if (line[iv] == 'a' || line[iv] == 'A') {
                  tempdata[ip] = 1;
                  sectempdata[ip] = 0;
                  plantsUsed[ig]++;
                  plantsUsed[ig+1]++;
                  As += 2;
                  if (measflag >= 1) {
                    if (highlowstatus[ip] == 1) highAs += 2;
                    if (highlowstatus[ip] == -1) lowAs += 2;
                  }
                  ip++;
                }
                if (line[iv] == 'b' || line[iv] == 'B') {
                  tempdata[ip] = 0;
                  sectempdata[ip] = 1;
                  plantsUsed[ig]++;
                  plantsUsed[ig+1]++;
                  Bs += 2;
                  if (measflag >= 1) {
                    if (highlowstatus[ip] == 1) highBs += 2;
                    if (highlowstatus[ip] == -1) lowBs += 2;
                  }
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
            if (As > Bs) segdises[i] = 0.5 + 2 * std::abs(((double) As / (As + Bs)) - 0.75);
            else segdises[i] = 0.5 + 2 * std::abs(((double) Bs / (As + Bs)) - 0.75);
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
            if (measflag == 1) {
              if (highAs + highBs > 0 && lowAs > 0) hightolows[ig] = (double)highAs * (lowAs + lowBs) / (lowAs * (highAs + highBs));
              else hightolows[ig] = 2 / delta;
              if (highAs + highBs > 0 && lowBs > 0) hightolows[ig+1] = (double)highBs * (lowAs + lowBs) / (lowBs * (highAs + highBs));
              else hightolows[ig+1] = 2 / delta;
            }
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
    if (specialfileflag == 1) {
      if (strlen(usedcodefile) > 1) PUD = fopen(usedcodefile, "w");
      else {
        printf("The configuration file did not specify the usedcodefile for missing data.\n");
        return(1);
      }
      for (i = 0; i < nLoci; i++) fprintf(PUD, "%d\t%d\n", i, plantsUsed[i]);
      fclose(PUD);
      if (strlen(segdisfile) > 1) SGD = fopen(segdisfile, "w");
      else {
        printf("The configuration file did not specify the segdisfile for segregation distortion.\n");
        return(1);
      }
      for (i = 0; i < nLoci; i++) fprintf(SGD, "%d\t%.4f\n", i, segdises[i]);
      fclose(SGD);
      fprintf(SPE, "Recombination fractions:\n");
      for (i = 0; i < nLoci; i++) {
        fprintf(SPE, ">%s:\n", names[i]);
        for (j = 0; j < nLoci; j++) {
          fprintf(SPE, "%.4f ", recombfrac[i][j]);
          if (j % 20 == 19) fprintf(SPE, "\n");
        }
        if (j % 20 != 19) fprintf(SPE, "\n");
      } 
      fflush(SPE);
    }
    for (i = 0; i < nLoci; i++) {
      for (j = 0; j < nLoci; j++) {
        if (recombfrac[i][j] != recombfrac[j][i]) fprintf(OUT, "reciprocals are unequal, ij = %.4f, ji = %.4f\n", recombfrac[i][j], recombfrac[j][i]);
      }
    }
  }
  imax = 0; //kmax is now set in the configuration file.
  if (queryflag == 1) {
    querycodes = (int *) calloc(nLoci, sizeof(int));
    i = 0;
    if (strlen(queryfile) > 1) QRY = fopen(queryfile, "r");
    else {
      printf("The configuration file did not specify the queryfile, which is needed.\n");
      return(1);
    }
    while (fgets(line, maxdatlinelength, QRY)) {
      fprintf(OUT, "Queried loci = %s", line);
      fflush(OUT);
      token = strtok(line, " ,|\t\n");
      if (token[0] == '/' && token[1] == '/') continue; //Comments are allowed; they begin with '//' or '#'.
      if (token[0] == '#') continue;
      queryarray[i] = (int)strtol(token, (char**)NULL, 10);
      querycodes[queryarray[i]] = 1; 
      i++;
      while ((token = strtok(NULL, " ,|\t\n"))) {
        if (token[0] == '/' && token[1] == '/') break;
        queryarray[i] = (int)strtol(token, (char **)NULL, 10);
        querycodes[queryarray[i]] = 1;
        i++;
        if (i >= 15) break; //We are limiting the number of queried loci to what can be permuted feasibly.
      }
      imax = i;
    }
    fclose(QRY);
  }
  if (imax > kmax) kmax = imax;
  if (kmax <= 10) coeffs = (int ***) malloc((1 + kmax) * sizeof(int **));
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
  if (queryflag == 1) { //This code handles a submitted marker order from the query file.
    queriedmaplength = 0;
    fprintf(OUT, "Queried recombination fractions:\n");
    for (j = 1; j < imax; j++) {
      queriedmaplength += recombfrac[queryarray[j-1]][queryarray[j]];
      fprintf(OUT, "%d %d %f\n", queryarray[j-1], queryarray[j], recombfrac[queryarray[j-1]][queryarray[j]]);
    }
    fprintf(OUT, "Sum of recombination fractions for the queried map: %f\n", queriedmaplength);
    //Now permute the loci, find the shortest local map, and also find the shortest local map with the original ends.
    minmaplength = 2000000000;
    for (i = 0; i < 20; i++) savedv[i] = -1;
    for (i = 0; i < nf[imax]; i++) {
      localmaplength = 0;
      for (j = 1; j < imax; j++) localmaplength += recombfrac[queryarray[coeffs[imax][i][j-1]]][queryarray[coeffs[imax][i][j]]];
      if (localmaplength < minmaplength) {
        minmaplength = localmaplength;
        mini = i;
        for (j = 0; j < imax; j++) savedv[j] = queryarray[coeffs[imax][i][j]];
      }
    }
    fprintf(OUT, "globally minimum local map length = %f\n", minmaplength);
    fprintf(OUT, "Globally minimal marker order is:");
    for (i = 0; i < imax; i++) fprintf(OUT, " %d", savedv[i]);
    fprintf(OUT, "\n");
    fflush(OUT);
  }
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
      if (recombfrac[i][j] < minofrfs[i] && segdises[j] <= segdisthreshold && fracUsedCode[j] > 1) {
        minofrfs[i] = recombfrac[i][j];
        nearestneighbors[i] = j;
      }
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
  //fprintf(OUT, "\nRequested recombfracs:\nrf[65][147] = %.4f\n", recombfrac[65][147]);
  //fprintf(OUT, " rf[210][160] = %.4f rf[160][98] = %.4f rf[98][424] = %.4f\n", recombfrac[210][160], recombfrac[160][98], recombfrac[98][424]);
  bandwidth = 1 / ((double) plantsTotal * ploidy);
  base = -bandwidth / 2;
  top = base + bandwidth;
  forePointers = (int *) malloc(nLoci * sizeof(int));
  aftPointers = (int *) malloc(nLoci * sizeof(int));
  foretestins = (int *) malloc(nLoci * sizeof(int)); //trial arrays for insertions
  afttestins = (int *) malloc(nLoci * sizeof(int));
  foretestatt = (int *) malloc(nLoci * sizeof(int));
  afttestatt = (int *) malloc(nLoci * sizeof(int));
  exemplars = (int *) malloc(nLoci * sizeof(int));
  for (i = 0; i < nLoci; i++) {
    forePointers[i] = -1; /*set fore and aft pointers to an invalid locus*/
    aftPointers[i] = -1;
    foretestatt[i] = -1;
    afttestatt[i] = -1;
    foretestins[i] = -1;
    afttestins[i] = -1;
  }
  if (dedupflag > 0) {
    printf("Deduplicating loci.\n");
    for (k = 0; k < nLoci; k++) {
      if (fracUsedCode[k] <= 1 || minofrfs[k] > maxminrf || segdises[k] > segdisthreshold) continue;
      for (m = k + 1; m < nLoci; m++) {
        if (fracUsedCode[m] <= 1 || minofrfs[m] > maxminrf || segdises[m] > segdisthreshold) continue;
        if (recombfrac[k][m] < top) {
          if (forePointers[k] == -1 || aftPointers[k] == -1) {
            if (forePointers[m] == -1 || aftPointers[m] == -1) {
              //fprintf(OUT, "k = %d m = %d forek = %d aftk = %d forem = %d aftm = %d in dedup\n", k, m, forePointers[k], aftPointers[k], forePointers[m], aftPointers[m]);
              //fflush(OUT);
              found = 0; //Prevent cycles.
              g = k;
              while (forePointers[g] >= 0) g = forePointers[g];
              if (g == m) found = 1;
              while (aftPointers[g] >= 0) { 
                if (aftPointers[g] == m) found = 1;
                g = aftPointers[g];
              }
              if (found == 0) {
                if (aftPointers[k] == -1) {
                  if (forePointers[m] == -1) {
                    aftPointers[k] = m;
                    forePointers[m] = k;
                  }
                  else {
                    if (aftPointers[m] == -1) { //Both are aft ends. Flip m's linkage group.
                      g = m;
                      while (g >= 0) {
                        temp = forePointers[g];
                        forePointers[g] = aftPointers[g];
                        aftPointers[g] = temp;
                        g = temp;
                      }
                      aftPointers[k] = m;
                      forePointers[m] = k;
                    }
                  }
                }
                else {
                  if (forePointers[k] == -1) {
                    if (aftPointers[m] == -1) {
                      forePointers[k] = m;
                      aftPointers[m] = k;
                    }
                  }
                  else {
                    if (forePointers[m] == -1) { //Both are fore ends. Flip m's linkage group.
                    g = m;
                      while (g >= 0) {
                        temp = aftPointers[g];
                        aftPointers[g] = forePointers[g];
                        forePointers[g] = temp;
                        g = temp;
                      }
                      forePointers[k] = m;
                      aftPointers[m] = k;
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
  for (i = 0; i < nLoci; i++) exemplars[i] = 0;
  nexemplars = 0;
  if (dedupflag > 0) {
    for (i = 0; i < nLoci; i++) {
      if (fracUsedCode[i] <= 1 || minofrfs[i] > maxminrf || segdises[i] > segdisthreshold) continue;
      if (forePointers[i] == -1 && aftPointers[i] == -1) {exemplars[i] = 1; nexemplars++;}
    }
    for (i = 0; i < nLoci; i++) {
      if (forePointers[i] == -1 && aftPointers[i] >= 0) { //A linkage group is not supposed to contain excluded loci, see about 55 lines up.
        g = i;
        pused = plantsUsed[g];
        savedg = g;
        while (g >= 0) {
          if (plantsUsed[g] > pused) {
            pused = plantsUsed[g];
            savedg = g;
          }
          g = aftPointers[g];
        }
        exemplars[savedg] = 1;
        nexemplars++;
      }
    } //Now the exemplars == 1 form a deduplicated subset.
    fprintf(OUT, "There are %d exemplar markers.\n", nexemplars);
    if (dedupoutflag == 1) { //print out the MapMaker entries for the deduplicated subset.
      if (strlen(dedupoutfile) > 1) DDP = fopen(dedupoutfile, "w");
      rewind(DAT);
      if (strcmp(inputformat, "mapmaker") == 0) {
        //*9
        //HBBABHHHHAAHHHHHHAHHH
        fgets(line, maxdatlinelength, DAT);
        fprintf(DDP, "%s", line);
        fgets(line, maxdatlinelength, DAT);
        fprintf(DDP, "%s", line);
        fgets(line, maxdatlinelength, DAT);
        fprintf(DDP, "%s", line);
        while(fgets(line, maxdatlinelength, DAT)) {
          if (line[0] == '\n' || line[0] == '\0' || line[0] == '#') {
            fprintf(DDP, "%s", line);
            continue;
          }
          if (line[0] == '/' && (line[1] == '/' || line[1] == '\n')) {
            fprintf(DDP, "%s", line);
            continue;
          }
          if (line[0] == '*') {
            for (i = 0; i < maxnamelength; i++) tempname[i] = '\0';
            for (iv = 1; iv < maxdatlinelength; iv++) {
              if (line[iv] == '\n' || line[iv] == '\0' || line[iv] == '#') break;
              if (line[iv] == '/' && (line[iv+1] == '/' || line[iv+1] == '\n')) break;
              tempname[iv-1] = line[iv];
            }
            if (strlen(tempname) > 0) {
              printflag = 0;
              if (combineparents == 0) strncat(tempname, "-1", 2);
              for (i = 0; i < nLoci; i++) {
                if (strcmp(names[i], tempname) == 0) {
                  if (exemplars[i] == 1) {
                    printflag = 1;
                    fprintf(DDP, "%s", line);
                    break;
                  }
                }
              }
            }
          }
          else {
            if (printflag == 1) fprintf(DDP, "%s", line);
          }
        }
      }
    }
  }
  if (dedupflag == 2) return(0); //Exit here if this flag merely indicates to output a deduplicated data file.
  //Erase all the zero-length linkage groups.  We will build with only the exemplars if dedupflag == 1, else with all eligible loci.
  for (i = 0; i < nLoci; i++) { //Erase all the zero-length linkage groups.  We will build with only the exemplars if dedupflag == 1, else with all.
    forePointers[i] = -1; /*set fore and aft pointers to an invalid locus*/
    aftPointers[i] = -1;
    foretestatt[i] = -1;
    afttestatt[i] = -1;
    foretestins[i] = -1;
    afttestins[i] = -1;
  }
  printf("Starting to order loci.\n");
  while (top < maxrf) {
    //fprintf(OUT, "top = %f forep[%d] = %d aftp[%d] = %d forep[%d] = %d aftp[%d] = %d\n", top, testlocusa, forePointers[testlocusa], testlocusa, aftPointers[testlocusa], testlocusb, forePointers[testlocusb], testlocusb, aftPointers[testlocusb]);
    //if (forePointers[testlocusa] >= 0) fprintf(OUT, "rf[%d][%d] = %f\n", forePointers[testlocusa], testlocusa, recombfrac[forePointers[testlocusa]][testlocusa]);
    //if (aftPointers[testlocusa] >= 0) fprintf(OUT, "rf[%d][%d] = %f\n", testlocusa, aftPointers[testlocusa], recombfrac[testlocusa][aftPointers[testlocusa]]);
    //if (forePointers[testlocusb] >= 0) fprintf(OUT, "rf[%d][%d] = %f\n", forePointers[testlocusb], testlocusb, recombfrac[forePointers[testlocusb]][testlocusb]);
    //if (aftPointers[testlocusb] >= 0) fprintf(OUT, "rf[%d][%d] = %f\n", testlocusb, aftPointers[testlocusb], recombfrac[testlocusb][aftPointers[testlocusb]]);
    for (k = 0; k < nLoci; k++) {
      if (fracUsedCode[k] <= 1 || minofrfs[k] > maxminrf || segdises[k] > segdisthreshold) continue;
      if (dedupflag == 1 && exemplars[k] == 0) continue;
      for (m = k + 1; m < nLoci; m++) {
        if (fracUsedCode[m] <= 1 || minofrfs[m] > maxminrf || segdises[m] > segdisthreshold) continue;
        if (dedupflag == 1 && exemplars[m] == 0) continue;
        if (forePointers[k] >= 0 && aftPointers[k] >= 0) {
          if (forePointers[m] >= 0 && aftPointers[m] >= 0) continue;
        } //Join the loci only if there is a least one free end among them.
        found = 0;
        g = k;
        while (forePointers[g] >= 0) g = forePointers[g];
        if (g == m) found = 1;
        while (aftPointers[g] >= 0) { //This should prevent the formation of all cycles.
          if (aftPointers[g] == m) found = 1;
          g = aftPointers[g];
        }
        if (found == 0 && recombfrac[k][m] >= base && recombfrac[k][m] < top) {
          if (forePointers[k] == -1 && aftPointers[k] == -1 && forePointers[m] == -1 && aftPointers[m] == -1) {
            aftPointers[k] = m; //Uncontested join of isolated loci.
            forePointers[m] = k;
          }
          else {
            endstatusk = 0;
            //If either k or m does not have a free end, its end status will be zero, and it will go to the insertion block below.
            if (forePointers[k] == -1) endstatusk = 1;
            else if (forePointers[forePointers[k]] == -1) endstatusk = 2;
            if (aftPointers[k] == -1) endstatusk = 3;
            else if (aftPointers[aftPointers[k]] == -1) endstatusk = 4;
            endstatusm = 0;
            if (forePointers[m] == -1) endstatusm = 1;
            else if (forePointers[forePointers[m]] == -1) endstatusm = 2;
            if (aftPointers[m] == -1) endstatusm = 3;
            else if (aftPointers[aftPointers[m]] == -1) endstatusm = 4;
            if (endstatusk > 0 && endstatusm > 0) { //Perform an annealed-end join.
              for (i = 0; i < 10; i++) mergearray[i] = -1;
              flankfore = -1;
              flankaft = -1;
              //We have the standard four cases of fore-fore, fore-aft, aft-fore, and aft-aft.
              if (endstatusk < 3 && endstatusm > 2) { //Join aftm to forek.
                if (endstatusm == 3) g = m;
                else if (endstatusm == 4) g = aftPointers[m];
                i = 0;
                while (i < 4 && g >= 0) {
                  mergearray[i] = g;
                  i++;
                  g = forePointers[g];
                }
                flankfore = g;
                if (endstatusk == 1) g = k;
                else if (endstatusk == 2) g = forePointers[k];
                while (i < 8 and g >= 0) { //Not balanced, biased toward group that contains k if m's group has fewer than four loci.
                  mergearray[i] = g;
                  i++;
                  g = aftPointers[g];
                }
                flankaft = g;
                anneallength = i; //Join has been set up.  One code block will do all joins below.
              }
              else if (endstatusk > 2 && endstatusm < 3) { //Join aftk to forem.
                if (endstatusk == 3) g = k;
                else if (endstatusk == 4) g = aftPointers[k];
                i = 0;
                while (i < 4 && g >= 0) {
                  mergearray[i] = g;
                  i++;
                  g = forePointers[g];
                }
                flankfore = g;
                if (endstatusm == 1) g = m;
                else if (endstatusm == 2) g = forePointers[m];
                while (i < 8 && g >= 0) { //Not necessarily balanced
                  mergearray[i] = g;
                  i++;
                  g = aftPointers[g];
                }
                flankaft = g;
                anneallength = i; //Join has been set up.
              }
              else if (endstatusk < 3 && endstatusm < 3) { //Flip group that contains k, then join aftk to forem.
                //Both loci are currently at the fore end of their linkage groups.
                //The strategy is to fill mergearray, set flanking loci, and then flip k's linkage group.
                if (endstatusk == 1) g = k;
                else if (endstatusk == 2) g = forePointers[k];
                i = 0;
                while (i < 4 && g >= 0) {
                  mergearray[i] = g;
                  i++;
                  g = aftPointers[g];
                } //At the end of this loop, g will be the flanking locus or -1.
                flankfore = g;
                if (endstatusm == 1) g = m;
                else if (endstatusm == 2) g = forePointers[m];
                while (i < 8 && g >= 0) { //Not necessarily balanced
                  mergearray[i] = g;
                  i++;
                  g = aftPointers[g];
                }
                flankaft = g;
                anneallength = i;
                if (endstatusk == 1) g = k; //Now flip k's linkage group.
                else if (endstatusk == 2) g = forePointers[k];
                while (g >= 0) {
                  temp = aftPointers[g];
                  aftPointers[g] = forePointers[g];
                  forePointers[g] = temp;
                  g = temp;
                } //The join has been set up.
              }
              else if (endstatusk > 2 && endstatusm > 2) {
                //The strategy is to fill mergearray, set flanking loci, and then flip m's linkage group.
                if (endstatusk == 3) g = k;
                else if (endstatusk == 4) g = aftPointers[k];
                i = 0;
                while (i < 4 && g >= 0) {
                  mergearray[i] = g;
                  i++;
                  g = forePointers[g];
                } //At the end of this loop, g will be the flanking locus or -1.
                flankfore = g;
                if (endstatusm == 3) g = m;
                else if (endstatusm == 4) g = aftPointers[m];
                while (i < 8 && g >= 0) {
                  mergearray[i] = g;
                  i++;
                  g = forePointers[g];
                } //At the end of this loop, g will be the flanking locus or -1.
                flankaft = g;
                anneallength = i;
                if (endstatusm == 3) g = m;
                else if (endstatusm == 4) g = aftPointers[m];
                while (g >= 0) {
                  temp = forePointers[g];
                  forePointers[g] = aftPointers[g];
                  aftPointers[g] = temp;
                  g = temp;
                } //The join has been set up.
              }
              //All four cases have been set up as flankfore-mergearray-flankaft.
              minmaplength = 2000000000;
              for (i = 0; i < 20; i++) savedv[i] = -1;
              for (i = 0; i < nf[anneallength]; i++) {
                if (flankfore >= 0) localmaplength = recombfrac[flankfore][mergearray[coeffs[anneallength][i][0]]];
                else localmaplength = 0;
                for (j = 1; j < anneallength; j++) {
                  localmaplength += recombfrac[mergearray[coeffs[anneallength][i][j-1]]][mergearray[coeffs[anneallength][i][j]]];
                }
                if (flankaft >= 0) localmaplength += recombfrac[mergearray[coeffs[anneallength][i][anneallength-1]]][flankaft];
                if (localmaplength < minmaplength) {
                  minmaplength = localmaplength;
                  mini = i;
                  for (j = 0; j < anneallength; j++) savedv[j] = mergearray[coeffs[anneallength][i][j]];
                }
              }
              if (flankfore >= 0) rfxatt = recombfrac[flankfore][savedv[0]];
              else rfxatt = 0;
              for (i = 0; i < anneallength - 1; i++) {
                if (recombfrac[savedv[i]][savedv[i+1]] > rfxatt) rfxatt = recombfrac[savedv[i]][savedv[i+1]];
              }
              if (flankaft >= 0) {
                if (recombfrac[savedv[anneallength-1]][flankaft] > rfxatt) rfxatt = recombfrac[savedv[anneallength-1]][flankaft];
              }
              commit = 0;
              if (rfxflag == 0) commit = 1;
              else if (rfxflag == 1 || rfxflag == 4) {
                if (rfxatt <= maxrf) commit = 1;
              }
              else if (rfxflag == 2 || rfxflag == 5) {
                if (rfxatt <= top) commit = 1;
              }
              else if (rfxflag == 3 || rfxflag == 6) {
                if (rfxatt <= rfxceiling) commit = 1;
              }
              else if (rfxflag == 4 || rfxflag == 8) {
                if (rfxatt <= top + delta) commit = 1;
              }
              //fprintf(OUT, "Value of commit for end-join = %d rfxflag = %d\n", commit, rfxflag);
              if (commit == 1) {
                if (flankfore >= 0) aftPointers[flankfore] = savedv[0];
                forePointers[savedv[0]] = flankfore;
                for (i = 0; i < anneallength - 1; i++) {
                  aftPointers[savedv[i]] = savedv[i+1];
                  forePointers[savedv[i+1]] = savedv[i];
                }
                aftPointers[savedv[anneallength-1]] = flankaft;
                if (flankaft >= 0) forePointers[flankaft] = savedv[anneallength-1];
              }
            }
            else { //Set up trial annealed end-join and whole mergers.
              if (forePointers[k] == -1 || aftPointers[k] == -1) { //ik has the free end(s), im has no free ends.
                ik = k;
                im = m;
              } 
              else if (forePointers[m] == -1 || aftPointers[m] == -1) {
                ik = m;
                im = k;
              }
              else fprintf(OUT, "Line 1177 has been reached, which should not happen.\n");
              foreendm = im;
              while (forePointers[foreendm] >= 0) foreendm = forePointers[foreendm];
              aftendm = im;
              while (aftPointers[aftendm] >= 0) aftendm = aftPointers[aftendm];
              for (i = 0; i < nLoci; i++) {
                foretestatt[i] = forePointers[i];
                afttestatt[i] = aftPointers[i];
              }
              for (i = 0; i < 10; i++) mergearray[i] = -1;
              flankfore = -1;
              flankaft = -1;
              if (recombfrac[ik][foreendm] < recombfrac[ik][aftendm]) { //Anneal-join ik to the fore end of the group that contains m.
                if (afttestatt[ik] >= 0) { //Flip ik's group.
                  if (foretestatt[ik] >= 0) fprintf(OUT, "Line 1190, ik has no free ends!\n");
                  g = ik;
                  while (g >= 0) {
                    temp = afttestatt[g];
                    afttestatt[g] = foretestatt[g];
                    foretestatt[g] = temp;
                    g = temp;
                  }
                } //Now ik must be at the aft end of its linkage group.
                i = 0;
                g = ik;
                while (i < 4 && g >= 0) {
                  mergearray[i] = g;
                  i++;
                  g = foretestatt[g];
                } //At the end of this loop, g will be the flanking locus or -1.
                flankfore = g;
                g = foreendm;
                while (i < 8 && g >= 0) {
                  mergearray[i] = g;
                  i++;
                  g = afttestatt[g];
                }
                anneallength = i;
                flankaft = g; //The join has been set up.
              }
              else { //Anneal-join ik to the aft end of the group that contains m.
                if (foretestatt[ik] >= 0) {//Flip ik's group.
                  if (afttestatt[ik] >= 0) fprintf(OUT, "Line 1218, ik has no free ends!\n");
                  g = ik;
                  while (g >= 0) {
                    temp = foretestatt[g];
                    foretestatt[g] = afttestatt[g];
                    afttestatt[g] = temp;
                    g = temp;
                  }
                } //Now ik must be at the fore end of its linkage group.
                i = 0;
                g = aftendm;
                while (i < 4 && g >= 0) {
                  mergearray[i] = g;
                  i++;
                  g = foretestatt[g];
                }
                flankfore = g;
                g = ik;
                while (i < 8 && g >= 0) {
                  mergearray[i] = g;
                  i++;
                  g = afttestatt[g];
                }
                anneallength = i;
                flankaft = g; //The join has been set up.
              }
              for (i = 0; i < 20; i++) savedv[i] = -1;
              minmaplength = 2000000000;
              for (i = 0; i < nf[anneallength]; i++) {
                if (flankfore >= 0) localmaplength = recombfrac[flankfore][mergearray[coeffs[anneallength][i][0]]];
                else localmaplength = 0;
                for (j = 1; j < anneallength; j++) {
                  localmaplength += recombfrac[mergearray[coeffs[anneallength][i][j-1]]][mergearray[coeffs[anneallength][i][j]]];
                }
                if (flankaft >= 0) localmaplength += recombfrac[coeffs[anneallength][i][anneallength-1]][flankaft];
                if (localmaplength < minmaplength) {
                  minmaplength = localmaplength;
                  mini = i;
                  for (j = 0; j < anneallength; j++) savedv[j] = mergearray[coeffs[anneallength][i][j]];
                }
              }
              if (flankfore >= 0) afttestatt[flankfore] = savedv[0];
              foretestatt[savedv[0]] = flankfore;
              for (i = 0; i < anneallength - 1; i++) {
                afttestatt[savedv[i]] = savedv[i+1];
                foretestatt[savedv[i+1]] = savedv[i];
              }
              afttestatt[savedv[anneallength-1]] = flankaft;
              if (flankaft >= 0) foretestatt[flankaft] = savedv[anneallength-1]; //The annealed-end trial is done.
              //Now try inserting the group that contains ik into the group that contains m.  Again, ik should be a free end.
              for (i = 0; i < nLoci; i++) {
                foretestins[i] = forePointers[i];
                afttestins[i] = aftPointers[i];
                insertions[i] = -1;
                mfores[i] = -1;
                mafts[i] = -1;
              }
              g = ik;
              while (foretestins[g] >= 0) g = foretestins[g];
              i = 0;
              while (g >= 0) { //Dissolve ik's linkage group.
                insertions[i] = g;
                temp = afttestins[g];
                foretestins[g] = -1;
                afttestins[g] = -1;
                g = temp;
                i++;
              }
              instodo = i;
              for (i = 0; i < instodo; i++) { //Find closest neighbors in group im for each locus liberated from ik's group.
                tripletlength = recombfrac[insertions[i]][foreendm] + recombfrac[foreendm][afttestins[foreendm]];
                //fprintf(OUT, "initial tripletlength = %.4f = %.4f + %.4f for %d %d %d\n", tripletlength, recombfrac[insertions[i]][foreendm], recombfrac[foreendm][afttestins[foreendm]], insertions[i], foreendm, afttestins[foreendm]);
                mfores[i] = -1;
                mafts[i] = foreendm;
                g = foreendm;
                while (afttestins[g] >= 0) { 
                  trial = recombfrac[g][insertions[i]] + recombfrac[insertions[i]][afttestins[g]];
                  if (trial < tripletlength) {
                    tripletlength = trial;
                    mfores[i] = g;
                    mafts[i] = afttestins[g];
                  }
                  g = afttestins[g];
                }
                trial = recombfrac[foretestins[g]][g] + recombfrac[g][insertions[i]];
                if (trial < tripletlength) {
                  tripletlength = trial;
                  mfores[i] = g;
                  mafts[i] = -1;
                }
              }
              for (i = 0; i < instodo; i++) {
                savedfore = -1;
                savedaft = -1;
                if (mfores[i] < 0) {
                  if (mafts[i] >= 0) {//The new locus goes before mafts[i].
                    g = mafts[i];
                    while (foretestins[g] >= 0) g = foretestins[g];
                    tripletlength = recombfrac[insertions[i]][g] + recombfrac[g][afttestins[g]];
                    savedfore = -1;
                    savedaft = g;
                    while (g != mafts[i]) {
                      trial = recombfrac[g][insertions[i]] + recombfrac[insertions[i]][afttestins[g]];
                      if (trial < tripletlength) {
                        tripletlength = trial;
                        savedfore = g;
                        savedaft = afttestins[g];
                      }
                      g = afttestins[g];
                    } //Upon finishing this loop, g is mafts[i], so we already have the trial length from foretestins[mafts[i]] to mafts[i].
                    if (savedfore >= 0) afttestins[savedfore] = insertions[i];
                    foretestins[insertions[i]] = savedfore;
                    afttestins[insertions[i]] = savedaft;
                    foretestins[savedaft] = insertions[i];
                  }
                  else {
                    fprintf(OUT, "Insertion was attempted to a single locus, which should not happen.\n");
                    return(5);
                  }
                }
                else {
                  if (mafts[i] >= 0) {
                    g = mfores[i];
                    tripletlength = recombfrac[g][insertions[i]] + recombfrac[afttestins[g]][insertions[i]];
                    savedfore = g;
                    savedaft = afttestins[g];
                    g = afttestins[g];
                    while (g != mafts[i]) {
                      trial = recombfrac[g][insertions[i]] + recombfrac[insertions[i]][afttestins[g]];
                      if (trial < tripletlength) {
                        tripletlength = trial;
                        savedfore = g;
                        savedaft = afttestins[g];
                      }
                      g = afttestins[g];
                    } //When done, g == mafts[i].
                    afttestins[savedfore] = insertions[i];
                    foretestins[insertions[i]] = savedfore;
                    afttestins[insertions[i]] = savedaft;
                    foretestins[savedaft] = insertions[i];
                  }
                  else { //The new locus goes behind mfores[i].
                    g = mfores[i];
                    tripletlength = recombfrac[foretestins[g]][g] + recombfrac[g][insertions[i]];
                    savedfore = g;
                    savedaft = afttestins[g]; //which can be -1
                    while (afttestins[g] >= 0) {
                      trial = recombfrac[g][insertions[i]] + recombfrac[insertions[i]][afttestins[g]];
                      if (trial < tripletlength) {
                        tripletlength = trial;
                        savedfore = g;
                        savedaft = afttestins[g];
                      }
                      g = afttestins[g];
                    }
                    trial = recombfrac[foretestins[g]][g] + recombfrac[g][insertions[i]];
                    if (trial < tripletlength) {
                      tripletlength = trial;
                      savedfore = g;
                      savedaft = -1;
                    }
                    afttestins[savedfore] = insertions[i];
                    foretestins[insertions[i]] = savedfore;
                    afttestins[insertions[i]] = savedaft;
                    if (savedaft >= 0) foretestins[savedaft] = insertions[i];
                  }
                }
              }
              g = ik;
              while (foretestatt[g] >= 0) g = foretestatt[g];
              found = 0; //This use of found differs from a cycle check.
              attlength = 0;
              rfxatt = 0;
              while (g >= 0) {
                if (g == im) {
                  found = 1;
                }
                if (afttestatt[g] >= 0) {
                   if (recombfrac[g][afttestatt[g]] > rfxatt) rfxatt = recombfrac[g][afttestatt[g]];
                   //fprintf(OUT, "attlength %.4f + recombfrac[%d][%d] %.4f\n", attlength, g, afttestatt[g], recombfrac[g][afttestatt[g]]);
                   attlength += recombfrac[g][afttestatt[g]];
                }
                g = afttestatt[g];
              }
              if (found == 0) fprintf(OUT, "The end-attached group at line 1321 contains %d but not %d\n", ik, im);
              g = ik;
              while (foretestins[g] >= 0) g = foretestins[g];
              found = 0;
              inslength = 0;
              rfxins = 0;
              while (g >= 0) {
                if (g == im) {
                  found = 1;
                }
                if (afttestins[g] >= 0) {
                  if (recombfrac[g][afttestins[g]] > rfxins) rfxins = recombfrac[g][afttestins[g]];
                  inslength += recombfrac[g][afttestins[g]];
                }
                g = afttestins[g];
              }
              if (found == 0) fprintf(OUT, "The post-insertion group at line 1333 contains %d but not %d\n", ik, im);
              commit = 0;
              if (rfxflag == 0) commit = 1;
              else if (rfxflag == 1 && rfxatt <= maxrf && rfxins <= maxrf) commit = 1;
              else if (rfxflag == 2 && rfxatt <= top && rfxins <= top) commit = 1;
              else if (rfxflag == 3 && rfxatt <= rfxceiling && rfxins <= rfxceiling) commit = 1;
              else if (rfxflag == 4 && rfxatt <= top + delta && rfxins <= top + delta) commit = 1; 
              else if (rfxflag == 5) {
                if (rfxatt <= maxrf || rfxins <= maxrf) commit = 1;
              }
              else if (rfxflag == 6) {
                if (rfxatt <= top || rfxins <= top) commit = 1;
              }
              else if (rfxflag == 7) {
                if (rfxatt <= rfxceiling || rfxins <= rfxceiling) commit = 1;
              }
              else if (rfxflag == 8) {
                if (rfxatt <= rfxins) {
                  if (rfxatt <= top + delta) commit = 1;
                }
                else {
                  if (rfxins <= top + delta) commit = 1;
                }
              }
              if (commit == 1) {
                if (inslength < attlength) {//Commit the insertion by copying foretestins and afttestins.
                  g = ik;
                  while (foretestins[g] >= 0) g = foretestins[g];
                  while (g >= 0) {
                    forePointers[g] = foretestins[g];
                    aftPointers[g] = afttestins[g];
                    g = afttestins[g];
                  }
                }
                else { //Commit the end-attachment by copying foretestatt and afttestatt.
                  g = ik;
                  while (foretestatt[g] >= 0) g = foretestatt[g];
                  while (g >= 0) {
                    forePointers[g] = foretestatt[g];
                    aftPointers[g] = afttestatt[g];
                    g = afttestatt[g];
                  }
                }
              }
            }
          }
          gycount = 0;
          gzcount = 0;
          for (g = 0; g < nLoci; g++) {
            if (forePointers[g] < 0) gycount++;
            if (aftPointers[g] < 0) gzcount++;
          }
          if (gycount != gzcount) return(-4);
          g = k;
          while (forePointers[g] >= 0) g = forePointers[g];
          while (aftPointers[g] >= 0) {
            if (recombfrac[g][aftPointers[g]] > maxrf) {
              fprintf(OUT, "Line 1607: Bad join for k = %d m = %d rf = %.4f maxrf = %.4f top = %.4f\n", k, m, recombfrac[g][aftPointers[g]], maxrf, top);
            }
            g = aftPointers[g];
          }
          g = m;
          while (forePointers[g] >= 0) g = forePointers[g];
          while (aftPointers[g] >= 0) {
            if (recombfrac[g][aftPointers[g]] > maxrf) {
              fprintf(OUT, "Line 1615: Bad join for k = %d m = %d rf = %.4f maxrf = %.4f top = %.4f\n", k, m, recombfrac[g][aftPointers[g]], maxrf, top);
            }
            g = aftPointers[g];
          }
        }
      }
    }
    base += bandwidth;
    top += bandwidth;
  } //End of the band loop.
  if (permuteon == 1) { //Permute a sliding window of kmax loci.
    for (k = 0; k < nLoci; k++) {
      if (forePointers[k] == -1) {
        stemlocusa = -1;
        ia = 0;
        g = k;
        while (ia < kmax) {
          if (g < 0) break;
          mergearray[ia] = g;
          g = aftPointers[g];
          ia++;
        }
        stemlocusb = g;
        while (mergearray[kmax-1] >= 0) { //Permute
          minmaplength = 2000000000;
          for (i = 0; i < nf[kmax]; i++) {
            localmaplength = 0;
            for (j = 1; j < kmax; j++) {
              localmaplength += recombfrac[mergearray[coeffs[kmax][i][j-1]]][mergearray[coeffs[kmax][i][j]]];
            }
            if (stemlocusa >= 0) localmaplength += recombfrac[mergearray[coeffs[kmax][i][0]]][stemlocusa];
            if (stemlocusb >= 0) localmaplength += recombfrac[mergearray[coeffs[kmax][i][kmax-1]]][stemlocusb];
            if (localmaplength < minmaplength) {
              minmaplength = localmaplength;
              mini = i;
              for (j = 0; j < kmax; j++) savedv[j] = mergearray[coeffs[kmax][i][j]];
            }
          }
          //Join
          if (stemlocusa >= 0) aftPointers[stemlocusa] = savedv[0];
          forePointers[savedv[0]] = stemlocusa;
          for (i = 1; i < kmax; i++) {
            aftPointers[savedv[i-1]] = savedv[i];
            forePointers[savedv[i]] = savedv[i-1];
          }
          aftPointers[savedv[kmax-1]] = stemlocusb;
          if (stemlocusb >= 0) forePointers[stemlocusb] = savedv[kmax-1];
          //Shift
          stemlocusa = savedv[0];
          for (i = 0; i < kmax - 1; i++) mergearray[i] = savedv[i+1];
          mergearray[kmax-1] = stemlocusb;
          stemlocusb = aftPointers[stemlocusb];
        }
      }
    }
  }
  if (runtests > 0) {
    fprintf(OUT, "Environment of test loci in descending map:\n");
    i = 5;
    printtestloci(testlocusa, testlocusb, i, recombfrac, forePointers, aftPointers, OUT);
  }
  mapLength = 0;
  mappedLoci = 0;
  nLinkageGroups = 0;
  for (i = 0; i < nLoci; i++) {
    if (forePointers[i] < 0 && aftPointers[i] >= 0) {
      nLinkageGroups++;
      g = i;
      mappedLoci++;
      while (aftPointers[g] >= 0) {
        mappedLoci++;
        mapLength += recombfrac[g][aftPointers[g]];
        g = aftPointers[g];
      }
    }
  }
  fprintf(OUT, "\nList of unmapped loci:\n");
  j = 0;
  for (i = 0; i < nLoci; i++) {
    if (forePointers[i] < 0 && aftPointers[i] < 0) {
      fprintf(OUT, " %d", i);
      j++;
      if (j % 20 == 0) fprintf(OUT, "\n");
    }
  }
  fprintf(OUT, "\n");
  if (j == 0) fprintf(OUT, "None\n");
  fprintf(OUT, "\nTotal mapped loci = %d map length = %.4f in %d linkage groups\n", mappedLoci, mapLength, nLinkageGroups);
  lgStarts = (int *) calloc(nLinkageGroups + 1, sizeof(int));
  lgEnds = (int *) calloc(nLinkageGroups + 1, sizeof(int));
  lgCounts = (int *) calloc(nLinkageGroups + 1, sizeof(int));
  lgLengths = (double *) calloc(nLinkageGroups + 1, sizeof(double));
  sortedlgCounts = (int *) calloc(nLinkageGroups + 1, sizeof(int));
  doneLinkageGroups = (int *) calloc(nLinkageGroups + 1, sizeof(int));
  k = 0;
  for (i = 0; i < nLoci; i++) {
    if (forePointers[i] < 0 && aftPointers[i] >= 0) {
      lgStarts[k] = i;
      g = i;
      lgCounts[k] = 1;
      while (aftPointers[g] >= 0) { //This has to be this way to collect recombination fractions.
        lgCounts[k]++;
        lgLengths[k] += recombfrac[g][aftPointers[g]];
        lgEnds[k] = aftPointers[g];
        g = aftPointers[g];
      }
      k++;
    }
  }
  for (i = 0; i < nLinkageGroups; i++) {
    sortedlgCounts[i] = lgCounts[i];
    doneLinkageGroups[i] = 0;
  }
  qsort(sortedlgCounts, nLinkageGroups, sizeof(int), compdown);
  for (i = 0; i < nLinkageGroups; i++) {
    for (j = 0; j < nLinkageGroups; j++) {
      if (lgCounts[j] == sortedlgCounts[i] && doneLinkageGroups[j] == 0) {
        fprintf(OUT, "linkage group %d\tstart = %d\tend = %d\tcount = %d\tlength = %f\n", i, lgStarts[j], lgEnds[j], lgCounts[j], lgLengths[j]);
        doneLinkageGroups[j] = 1;
        break;
      }
    }
  }
  minlgdistany = (double *) malloc(nLinkageGroups * sizeof(double));
  minlgothergroup = (int *) malloc(nLinkageGroups * sizeof(int));
  minlgdists = (double **) malloc(nLinkageGroups * sizeof(double *));
  for (i = 0; i < nLinkageGroups; i++) minlgdists[i] = (double *) malloc(nLinkageGroups * sizeof(double));
  orignlgroups = nLinkageGroups;
  for (i = 0; i < nLinkageGroups; i++) doneLinkageGroups[i] = 0; //Reload to print out the map itself.
  for (i = 0; i < nLoci; i++) membership[i] = -1;
  for (i = 0; i < nLinkageGroups; i++) {
    for (j = 0; j < nLinkageGroups; j++) minlgdists[i][j] = 1;
  }
  fprintf(OUT, "\nThere are %d linkage groups in descending order of locus count:\n", nLinkageGroups);
  mq = 0;
  for (i = 0; i < nLinkageGroups; i++) {
    for (j = 0; j < nLinkageGroups; j++) {
      if (lgCounts[j] == sortedlgCounts[i] && doneLinkageGroups[j] == 0) {
        fprintf(OUT, "linkage group %d\tstart = %d\tend \t%d\tcount = %d\tlength = %f\n", i, lgStarts[j], lgEnds[j], lgCounts[j], lgLengths[j]);
        g = lgStarts[j];
        cudist = 0.0;
        while (aftPointers[g] >= 0) {
          if (queryflag == 1 && querycodes[g] == 1) {
            mapqueries[mq] = g;
            mq++;
          }
          fprintf(OUT, "%s\t%d\t%d\t%f\t%f\t%f\t%s\t%f", names[g], g, aftPointers[g], recombfrac[g][aftPointers[g]], cudist, minofrfs[g], names[nearestneighbors[g]], segdises[g]);
          if (measflag == 1) fprintf(OUT, "\t%6.4f", hightolows[g]);
          fprintf(OUT, "\n");
          membership[g] = i;
          cudist += recombfrac[g][aftPointers[g]];
          g = aftPointers[g];
        }
        if (queryflag == 1 && querycodes[g] == 1) {
          mapqueries[mq] = g;
          mq++;
        }
        fprintf(OUT, "%s\t%d\t%d\t-\t%f\t%f\t%s\t%f", names[g], g, aftPointers[g], cudist, minofrfs[g], names[nearestneighbors[g]], segdises[g]);
        if (measflag == 1) fprintf(OUT, "\t%f", hightolows[g]);
        fprintf(OUT, "\n");
        membership[g] = i;
        doneLinkageGroups[j] = 1;
        break;
      }
    }
  }
  if (queryflag == 1) {
    fprintf(OUT, "mq = %d\n", mq);
    queriedmaplength = 0;
    fprintf(OUT, "Recombination fractions for queried markers in the preceding map:\n");
    for (i = 1; i < mq; i++) {
      queriedmaplength += recombfrac[mapqueries[i-1]][mapqueries[i]];
      fprintf(OUT, "%d %d %f\n", mapqueries[i-1], mapqueries[i], recombfrac[mapqueries[i-1]][mapqueries[i]]);
    }
    fprintf(OUT, "Sum of recombination fractions for these queried markers: %f\n", queriedmaplength);
  }
  //Special block to test proximity of two linkage groups.
  if (runtests > 0) {
    fprintf(OUT, "\nMutual recombination fractions of linkage groups containing loci %d and %d:\n", testlocusa, testlocusb);
    starta = testlocusa;
    while (forePointers[starta] >= 0) starta = forePointers[starta];
    startb = testlocusb;
    while (forePointers[startb] >= 0) startb = forePointers[startb];
    g = starta;
    while (g >= 0) {
      fprintf(OUT, "%d:", g);
      h = startb;
      while (h >= 0) {
        fprintf(OUT, " %f", recombfrac[g][h]);
        h = aftPointers[h];
      }
      fprintf(OUT, "\n");
      g = aftPointers[g];
    }
    fprintf(OUT, "\n");
  }
  for (i = 0; i < nLoci; i++) {
    if (membership[i] >= 0) {
      for (j = 0; j < nLoci; j++) {
        if (membership[j] >= 0) {
          if (recombfrac[i][j] < minlgdists[membership[i]][membership[j]]) {
            minlgdists[membership[i]][membership[j]] = recombfrac[i][j];
            minlgdists[membership[j]][membership[i]] = recombfrac[i][j]; //Keep all members synchronized.
          }
        }
      }
    }
  }
  for (i = 0; i < nLinkageGroups; i++) {
    minlgdistany[i] = 100.0;
    minlgothergroup[i] = -1;
  }
  fprintf(OUT, "\nMinimum recombination fractions between linkage groups:\n");
  for (i = 0; i < nLinkageGroups; i++) {
    fprintf(OUT, "%d:", i);
    for (j = 0; j < nLinkageGroups; j++) {
      fprintf(OUT, "\t%.4f", minlgdists[i][j]);
      if (j % 10 == 9) fprintf(OUT, "\n");
      if (i != j && minlgdists[i][j] < minlgdistany[i]) {
        minlgdistany[i] = minlgdists[i][j];
        minlgothergroup[i] = j;
      }
    }
    fprintf(OUT, "\n");
  }
  fprintf(OUT, "\n");
  for (i = 0; i < nLinkageGroups; i++) fprintf(OUT, "minlgdistany[%d] = %f minlgothergroup[%d] = %d\n", i, minlgdistany[i], i, minlgothergroup[i]);
  fprintf(OUT, "\n");
  for (i = 0; i < nLinkageGroups; i++) {
    for (j = 0; j < nLinkageGroups; j++) {
      if (i != j && minlgdists[i][j] < maxrf) {
        fprintf(OUT, "Possible freeze-out involves linkage groups %d and %d at minimum recombination %.4f\n", i, j, minlgdists[i][j]);
      }
    }
  }
  if (outaflag == 1) {
    if (strlen(outafile) > 1) ACT = fopen(outafile, "r");
    else {
      printf("The configuration file did not specify outafile, which is needed here.\n");
      return(1);
    }
    actmap = (int**) malloc(2 * baseNumber * sizeof(int*));
    for (i = 0; i < 2 * baseNumber; i++) {
      actmap[i] = (int*) malloc((nLoci + 1) * sizeof(int));
      for (j = 0; j < nLoci; j++) actmap[i][j] = -1;
    }
    k = -1;
    while (fgets(line, maxdatlinelength, ACT)) {
      i = 0;
      token = strtok(line, " ,|\t\n");
      if (token[0] == '/' && token[1] == '/') continue; //Comments are allowed; they begin with '//' or '#'.
      if (token[0] == '#') continue;
      strcpy(vars[i], token);
      while ((token = strtok(NULL, " ,|\t\n"))) {
        i++;
        strcpy(vars[i], token);
      }
      if (strcmp(vars[3], "TELO") == 0) {
        k++;
        j = 0;
      }
      else {
        actmap[k][j] = (int)strtol(vars[1], (char**)NULL, 10);
        j++;
      }
    }
    actlength = 0;
    for (i = 0; i < 2 * baseNumber; i++) {
      for (j = 1; j < nLoci; j++) {
        if (actmap[i][j] == -1) break;
        actlength += recombfrac[actmap[i][j-1]][actmap[i][j]];
      }
    }
    fprintf(OUT, "Map length with actual marker order and the recombination data is %f\n", actlength);
  }
  if (continuation == 0) {
    currtime = clock();
    duration = (double) (currtime - timeused2) / CLOCKS_PER_SEC;
    fprintf(OUT, "done with program, time used = %f\n", duration);
  }
  else {
    if (runtests > 0) {
      fprintf(OUT, "DIAG: before stripping misfit loci, value of rfthreshold = %f\n", rfthreshold);
      fprintf(OUT, "Environment of test loci just before stripping misfit loci:\n");
      i = 5;
      printtestloci(testlocusa, testlocusb, i, recombfrac, forePointers, aftPointers, OUT);
    }
    for (l = 0; l < nfloop; l++) {
      for (i = 0; i < nLoci; i++) { //Splice out any loci that recombine too much with their neighbors.
        nexcised = 0;
        if (forePointers[i] == -1) {
          if (aftPointers[i] >= 0) {
            if (recombfrac[i][aftPointers[i]] > rfthreshold) {
              forePointers[aftPointers[i]] = -1;
              aftPointers[i] = -1;
              nexcised++;
            }
          }
        }
        else if (forePointers[i] >= 0) {
          if (aftPointers[i] == -1) {
            if (recombfrac[i][forePointers[i]] > rfthreshold) {
              aftPointers[forePointers[i]] = -1;
              forePointers[i] = -1;
              nexcised++;
            }
          }
          else {
            if (excisionflag == 0) {
              if (recombfrac[i][forePointers[i]] > rfthreshold && recombfrac[i][aftPointers[i]] > rfthreshold) {
                aftPointers[forePointers[i]] = aftPointers[i];
                forePointers[aftPointers[i]] = forePointers[i];
                aftPointers[i] = -1;
                forePointers[i] = -1;
                nexcised++;
              }
            }
            else {
              if (recombfrac[i][forePointers[i]] + recombfrac[i][aftPointers[i]] - recombfrac[forePointers[i]][aftPointers[i]] > rfthreshold) {
                aftPointers[forePointers[i]] = aftPointers[i];
                forePointers[aftPointers[i]] = forePointers[i];
                aftPointers[i] = -1;
                forePointers[i] = -1;
                nexcised++;
              }
              else if (aftPointers[aftPointers[i]] >= 0 && recombfrac[i][forePointers[i]] + recombfrac[aftPointers[i]][aftPointers[aftPointers[i]]] - recombfrac[forePointers[i]][aftPointers[aftPointers[i]]] > rfthreshold) {
                aftPointers[forePointers[i]] = aftPointers[aftPointers[i]];
                forePointers[aftPointers[aftPointers[i]]] = forePointers[i];
                forePointers[aftPointers[i]] = -1;
                aftPointers[aftPointers[i]] = -1;
                forePointers[i] = -1;
                aftPointers[i] = -1;
                nexcised += 2;
              }
            }
          }
        }
      }
      fprintf(OUT, "Number of markers excised = %d\n", nexcised);
      if (runtests > 0) {
        fprintf(OUT, "Environment of test loci just after stripping misfit loci:\n");
        i = 5;
        printtestloci(testlocusa, testlocusb, i, recombfrac, forePointers, aftPointers, OUT);
      }
      //Cleave linkage groups at gaps that exceed maxrf.  Such gaps can form upon excision of misfit markers.
      //Cleave before dismembering too-short groups, since cleavage might produce short groups.
      for (i = 0; i < nLoci; i++) {
        if (aftPointers[i] >= 0) {
          if (recombfrac[i][aftPointers[i]] >= maxrf) {
            forePointers[aftPointers[i]] = -1;
            aftPointers[i] = -1;
          }
        }
      }
      for (i = 0; i < nLoci; i++) {
        if (forePointers[i] == -1 && aftPointers[i] >= 0) {
          g = i;
          currcount = 0;
          while (g >= 0) {
            currcount++;
            g = aftPointers[g];
          }
          //Dismember the too-short and presumed extraneous linkage group if it has close markers in another linkage group.
          //The hard-coded 6 is needed to prevent a segmentation fault when filling mergearray, since a two-locus linkage
          //group is now permitted.
          if (currcount < 6 || (currcount < mincount && minlgdistany[membership[i]] < tooclose)) {
            fprintf(OUT, "Linkage group %d starting with locus %d is being dismembered with currcount = %d and mincount = %d.\n", membership[i], i, currcount, mincount);
            g = i;
            while (g >= 0) {
              j = aftPointers[g];
              forePointers[g] = -1;
              aftPointers[g] = -1;
              g = j;
            }
          }
        }
      }
      strcpy(printstring, "after dismembering short linkage groups");
      printStats(printstring, nLoci, recombfrac, forePointers, aftPointers, OUT);
      nLinkageGroups = 0; 
      //fprintf(OUT, "Unordered map after removing high-recombining loci:\n");
      for (i = 0; i < nLoci; i++) {
        if (forePointers[i] == -1 && aftPointers[i] >= 0) {
          nLinkageGroups++;
//          fprintf(OUT, "Linkage group beginning with locus %d:\n", i);
//          cudist = 0.0;
//          g = i;
//          while (aftPointers[g] >= 0) {
//            fprintf(OUT, "%s\t%d\t%d\t%f\t%f\t%f\t%f", names[g], g, aftPointers[g], recombfrac[g][aftPointers[g]], cudist, minofrfs[g], segdises[g]);
//            if (measflag == 1) fprintf(OUT, "\t%6.4f", hightolows[g]);
//            fprintf(OUT, "\n");
//            cudist += recombfrac[g][aftPointers[g]];
//            g = aftPointers[g];
//          }
//          fprintf(OUT, "%s\t%d\t%d\t-\t%f\t%f\t%f", names[g], g, aftPointers[g], cudist, minofrfs[g], segdises[g]);
//          if (measflag == 1) fprintf(OUT, "\t%6.4f", hightolows[g]);
//          fprintf(OUT, "\n");
        }
      }
      pool = (int *) calloc(2 * nLinkageGroups + 1, sizeof(int));
      for (i = 0; i < 2 * nLinkageGroups + 1; i++) pool[i] = -1;
      k = 0;
      for (i = 0; i < nLoci; i++) {
        found = 0;
        if (forePointers[i] == -1 && aftPointers[i] >= 0) found = 1;
        if (forePointers[i] >= 0 && aftPointers[i] == -1) found = 2;
        if (found >= 1) {
          pool[k] = i;
          if (found == 1) {fprintf(OUT, "pool[%d] = %d, which is a fore end with fP = %d and aP = %d\n", k, pool[k], forePointers[pool[k]], aftPointers[pool[k]]);} 
          if (found == 2) {fprintf(OUT, "pool[%d] = %d, which is an aft end with fP = %d and aP = %d\n", k, pool[k], forePointers[pool[k]], aftPointers[pool[k]]);}
          k++;
        }
      }
      km = k;
      fprintf(OUT, "Pool has been filled. There are %d ends and %d linkage groups.\n\n", k, nLinkageGroups);
      if (runtests > 0) {
        fprintf(OUT, "Environment of test loci after filling the pool:\n");
        i = 5;
        printtestloci(testlocusa, testlocusb, i, recombfrac, forePointers, aftPointers, OUT);
      }
      base = -bandwidth / 2;
      top = base + bandwidth;
      while (top < maxrf) {
        for (k = 0; k < km; k++) { //Iterate over the pool.
          for (m = k + 1; m < km; m++) { //m must be greater than k.
            //fprintf(OUT, "DIAG at line 2497: km = %d k = %d m = %d pool[k] = %d pool[m] = %d", km, k, m, pool[k], pool[m]);
            //fprintf(OUT, " rf[%d][%d] = %f\n", pool[k], pool[m], recombfrac[pool[k]][pool[m]]);
            //fflush(OUT);
            if (recombfrac[pool[k]][pool[m]] >= base && recombfrac[pool[k]][pool[m]] < top) { 
              //Test if pool[k] and pool[m] are still free ends.
              if ((forePointers[pool[k]] == -1 && aftPointers[pool[k]] >= 0) || (forePointers[pool[k]] >= 0 && aftPointers[pool[k]] == -1)) {
                if ((forePointers[pool[m]] == -1 && aftPointers[pool[m]] >= 0) || (forePointers[pool[m]] >= 0 && aftPointers[pool[m]] == -1)) {
                  found = 0; //Prevent cycles.
                  g = pool[k];
                  while (forePointers[g] >= 0) {
                    if (g == pool[m]) {found = 1; break;}
                    g = forePointers[g];
                  }
                  while (g >= 0) {
                    if (g == pool[m]) {found = 1; break;}
                    g = aftPointers[g];
                  }
                  //fprintf(OUT, "DIAG at line 2513, found = %d for pool[k] = %d and pool[m] = %d\n", found, pool[k], pool[m]);
                  //fprintf(OUT, "DIAG at line 2514, fP[%d] = %d aP[%d] = %d fP[%d] = %d aP[%d] = %d\n", pool[k], forePointers[pool[k]], pool[k], aftPointers[pool[k]], pool[m], forePointers[pool[m]], pool[m], aftPointers[pool[m]]);
                  if (found == 0) { //We can anneal-join pool[k] and pool[m].
                    if (aftPointers[pool[k]] == -1 && forePointers[pool[k]] >= 0) {
                      if (forePointers[pool[m]] == -1 && aftPointers[pool[m]] >= 0) {
                        ik = pool[k];
                        im = pool[m];
                      }
                      else if (forePointers[pool[m]] >= 0 && aftPointers[pool[m]] == -1) { //Flip pool[m]'s group.
                        g = pool[m];
                        while (g >= 0) {
                          temp = forePointers[g];
                          forePointers[g] = aftPointers[g];
                          aftPointers[g] = temp;
                          g = temp;
                        }
                        ik = pool[k];
                        im = pool[m];
                      }
                      else {
                        fprintf(OUT, "Locus %d is not an end in attempted join at line 2261.\n", pool[m]);
                        continue; //Don't join pool[k] with pool[m].
                      }
                    }
                    else if (aftPointers[pool[k]] >= 0 && forePointers[pool[k]] == -1) {
                      if (aftPointers[pool[m]] == -1 && forePointers[pool[m]] >= 0) {
                        ik = pool[m];
                        im = pool[k];
                      }
                      else if (forePointers[pool[m]] == -1 && aftPointers[pool[m]] >= 0) { //Flip pool[k]'s group.
                        g = pool[k];
                        while (g >= 0) {
                          temp = aftPointers[g];
                          aftPointers[g] = forePointers[g];
                          forePointers[g] = temp;
                          g = temp;
                        }
                        ik = pool[k];
                        im = pool[m];
                      }
                      else {
                        fprintf(OUT, "Locus %d is not an end in attempted join at line 2291.\n", pool[m]);
                        continue; //Don't join pool[k] with pool[m].
                      }
                    }
                    fprintf(OUT, "Join proposed of aft %d with fore %d with base = %.6f and top = %.6f\n", ik, im, base, top);
                    fflush(OUT);
                    //We assume that the minimum accepted linkage-group length is at least 5.
                    //Join ik with im using end-annealing.
                    for (i = 0; i < 10; i++) mergearray[i] = -1;
                    g = ik;
                    for (i = 0; i < 5; i++) {
                      mergearray[i] = g;
                      g = forePointers[g];
                    }
                    flankfore = g;
                    g = im;
                    for (i = 5; i < 10; i++) {
                      mergearray[i] = g;
                      g = aftPointers[g];
                    }
                    flankaft = g;
                    minmaplength = 2000000000;
                    for (i = 0; i < nf[10]; i++) {
                      if (flankfore >= 0) localmaplength = recombfrac[flankfore][mergearray[coeffs[10][i][0]]];
                      else localmaplength = 0;
                      for (j = 1; j < anneallength; j++) {
                        localmaplength += recombfrac[mergearray[coeffs[10][i][j-1]]][mergearray[coeffs[10][i][j]]];
                      }
                      if (flankaft >= 0) localmaplength += recombfrac[mergearray[coeffs[10][i][anneallength-1]]][flankaft];
                      if (localmaplength < minmaplength) {
                        minmaplength = localmaplength;
                        mini = i;
                        for (j = 0; j < 10; j++) savedv[j] = mergearray[coeffs[10][i][j]];
                      }
                    }
                    if (flankfore >= 0) aftPointers[flankfore] = savedv[0];
                    forePointers[savedv[0]] = flankfore;
                    for (i = 0; i < 9; i++) {
                      aftPointers[savedv[i]] = savedv[i+1];
                      forePointers[savedv[i+1]] = savedv[i];
                    }
                    aftPointers[savedv[9]] = flankaft;
                    if (flankaft >= 0) forePointers[flankaft] = savedv[9];
                  }
                }
              }
            }
          }
        }
        base += bandwidth;
        top += bandwidth;
      }
      if (runtests > 0) {
        fprintf(OUT, "Environment of test loci after cleaning and joining ends:\n");
        k = 5;
        printtestloci(testlocusa, testlocusb, k, recombfrac, forePointers, aftPointers, OUT);
      }
      fprintf(OUT, "Done with cleaned-end joins loop.\n"); fflush(NULL);
      if (permuteon == 2) { //Permute a sliding window of kmax loci beginning at the start of each linkage group.
        printf("Locally permuting markers.\n");
        for (k = 0; k < nLoci; k++) {
          if (forePointers[k] == -1 && aftPointers[k] >= 0) {
            stemlocusa = -1;
            ia = 0;
            g = k;
            while (ia < kmax) {
              if (g < 0) break;
              mergearray[ia] = g;
              g = aftPointers[g];
              ia++;
            }
            stemlocusb = g;
            while (mergearray[kmax-1] >= 0) { 
              //Permute
              minmaplength = 2000000000;
              for (i = 0; i < nf[kmax]; i++) {
                localmaplength = 0;
                for (j = 1; j < kmax; j++) {
                  localmaplength += recombfrac[mergearray[coeffs[kmax][i][j-1]]][mergearray[coeffs[kmax][i][j]]];
                }
                if (stemlocusa >= 0) localmaplength += recombfrac[mergearray[coeffs[kmax][i][0]]][stemlocusa];
                if (stemlocusb >= 0) localmaplength += recombfrac[mergearray[coeffs[kmax][i][kmax-1]]][stemlocusb];
                if (localmaplength < minmaplength) {
                  minmaplength = localmaplength;
                  mini = i;
                  for (j = 0; j < kmax; j++) savedv[j] = mergearray[coeffs[kmax][i][j]];
                }
              }
              //Join
              if (stemlocusa >= 0) aftPointers[stemlocusa] = savedv[0];
              forePointers[savedv[0]] = stemlocusa;
              for (i = 1; i < kmax; i++) {
                aftPointers[savedv[i-1]] = savedv[i];
                forePointers[savedv[i]] = savedv[i-1];
              }
              aftPointers[savedv[kmax-1]] = stemlocusb;
              if (stemlocusb >= 0) forePointers[stemlocusb] = savedv[kmax-1];
              //Shift
              stemlocusa = savedv[0];
              for (i = 0; i < kmax - 1; i++) mergearray[i] = savedv[i+1];
              mergearray[kmax-1] = stemlocusb;
              if (stemlocusb >= 0) stemlocusb = aftPointers[stemlocusb];
            }
          }
        }
        printf("Done with %d-locus permutation of cleaned map.\n", kmax);
        if (runtests > 0) {
          fprintf(OUT, "Environment of test loci after local permutation:\n");
          k = 5;
          printtestloci(testlocusa, testlocusb, k, recombfrac, forePointers, aftPointers, OUT);
        }
      }
      for (i = 0; i < nLinkageGroups + 1; i++) { //Recount the linkage groups and redetermine their ends.
        lgStarts[i] = -1;
        lgEnds[i] = -1;
        lgLengths[i] = 0;
        lgCounts[i] = 0;
        sortedlgCounts[i] = 0;
      }
      k = 0;
      for (i = 0; i < nLoci; i++) {
        if (forePointers[i] < 0 && aftPointers[i] >= 0) {
          lgStarts[k] = i;
          g = i;
          lgCounts[k] = 1;
          sortedlgCounts[k] = 1;
          while (aftPointers[g] >= 0) { //This has to be this way to collect recombination fractions.
            lgCounts[k]++;
            sortedlgCounts[k]++;
            lgLengths[k] += recombfrac[g][aftPointers[g]];
            lgEnds[k] = aftPointers[g];
            g = aftPointers[g];
          }
          k++;
        }
      }
      nLinkageGroups = k;
      qsort(sortedlgCounts, nLinkageGroups, sizeof(int), compdown);
      if (continuation >= 1) { //Output the map as it stands.
        if (runtests > 0) {
          fprintf(OUT, "Environment of test loci in misfit map:\n");
          i = 5;
          printtestloci(testlocusa, testlocusb, i, recombfrac, forePointers, aftPointers, OUT);
        }
        fprintf(OUT, "\nThere are %d linkage groups after removing misfit loci:\n", nLinkageGroups); //Mutation alert: j is re-used immediately below.
        for (i = 0; i < nLinkageGroups; i++) doneLinkageGroups[i] = 0;
        for (i = 0; i < nLinkageGroups; i++) { //We are using values of sortedlgCounts from above.
          for (j = 0; j < nLinkageGroups; j++) {
            if (lgCounts[j] == sortedlgCounts[i] && doneLinkageGroups[j] == 0) {
              fprintf(OUT, "linkage group %d\tstart = %d\tend \t%d\tcount = %d\tlength = %f\n", i, lgStarts[j], lgEnds[j], lgCounts[j], lgLengths[j]);
              cudist = 0.0;
              g = lgStarts[j];
              while (aftPointers[g] >= 0) {
                fprintf(OUT, "%s\t%d\t%d\t%f\t%f\t%f\t%s\t%f", names[g], g, aftPointers[g], recombfrac[g][aftPointers[g]], cudist, minofrfs[g], names[nearestneighbors[g]], segdises[g]);
                if (measflag == 1) fprintf(OUT, "\t%f", hightolows[g]);
                fprintf(OUT, "\n");
                cudist += recombfrac[g][aftPointers[g]];
                g = aftPointers[g];
              }
              fprintf(OUT, "%s\t%d\t%d\t-\t%f\t%f\t%s\t%f", names[g], g, aftPointers[g], cudist, minofrfs[g], names[nearestneighbors[g]], segdises[g]);
              if (measflag == 1) fprintf(OUT, "\t%f", hightolows[g]);
              fprintf(OUT, "\n");
              doneLinkageGroups[j] = 1;
              break; //This is essential for the consecutive numbering of the linkage groups.
            }
          }
        }
      }
      fprintf(OUT, "\nList of unmapped loci:\n");
      j = 0;
      for (i = 0; i < nLoci; i++) {
        if (forePointers[i] < 0 && aftPointers[i] < 0) {
          fprintf(OUT, " %d", i);
          j++;
          if (j%20 == 0) fprintf(OUT, "\n");
        }
      }
      fprintf(OUT, "\n");
      if (j == 0) fprintf(OUT, "None\n");
      if (continuation == 1) {
        currtime = clock();
        duration = (double) (currtime - timeused2) / CLOCKS_PER_SEC;
        fprintf(OUT, "done with program, time used = %f\n", duration);
        return(1);
      }
      else { //Force all unbound loci next to their nearest neighbor that is within a linkage group.
        printf("Inserting unbound markers onto the map.\n");
        for (kr = 0; kr < insrounds; kr++) {
          for (i = 0; i < nLoci; i++) {
            if (forePointers[i] < 0 && aftPointers[i] < 0) {
              //Exclude loci that exceed segregation distortion or missing data thresholds.
              if (fracUsedCode[i] <= 1 || minofrfs[i] > maxminrf || segdises[i] > segdisthreshold) continue;
              if (exclusionflag == 1 && exemplars[i] == 0) continue;
              mindist = 1000000000;
              savedj = -1;
              for (j = 0; j < nLoci; j++) {
                if (forePointers[j] < 0 && aftPointers[j] < 0) continue; //Insert into an existing linkage group.
                //Handle zero-recombining loci separately, without local permutation.
                //All this looks overly complicated to insert a zero-recombining locus, but missing data can affect recombination
                //fraction with other nearby loci.
                if (recombfrac[i][j] < 0.00000001) { //Avoid test failure due to round-off error.
                  mindist = recombfrac[i][j];
                  if (forePointers[j] == -1) {
                    trialdista = recombfrac[i][j] + recombfrac[j][aftPointers[j]];
                    trialdistb = recombfrac[i][j] + recombfrac[i][aftPointers[j]];
                    if (trialdista < trialdistb) {
                      forePointers[j] = i;
                      aftPointers[i] = j;
                    }
                    else {
                      forePointers[aftPointers[j]] = i;
                      aftPointers[i] = aftPointers[j];
                      aftPointers[j] = i;
                      forePointers[i] =j;
                    }
                  }
                  else if (aftPointers[j] == -1) {
                    trialdista = recombfrac[i][j] + recombfrac[j][forePointers[j]];
                    trialdistb = recombfrac[i][j] + recombfrac[i][forePointers[j]];
                    if (trialdista < trialdistb) {
                      aftPointers[j] = i;
                      forePointers[i] = j;
                    }
                    else {
                      aftPointers[forePointers[j]] = i;
                      forePointers[i] = forePointers[j];
                      forePointers[j] = i;
                      aftPointers[i] = j;
                    }
                  }
                  else { //Insert into the interior of the linkage group.
                    trialdista = recombfrac[forePointers[j]][i] + recombfrac[i][j] + recombfrac[j][aftPointers[j]];
                    trialdistb = recombfrac[forePointers[j]][j] + recombfrac[i][j] + recombfrac[i][aftPointers[j]];
                    if (trialdista < trialdistb) {
                      aftPointers[forePointers[j]] = i;
                      forePointers[i] = forePointers[j];
                      forePointers[j] = i;
                      aftPointers[i] = j;
                    }
                    else {
                      forePointers[aftPointers[j]] = i;
                      aftPointers[i] = aftPointers[j];
                      aftPointers[j] = i;
                      forePointers[i] = j;
                    }
                  }
                  break;
                }
                //The following statements refer to loci that recombine with all other loci at least once.
                if (minflag == 0) { //Minimize added map length.
                  if (forePointers[j] >= 0 && aftPointers[j] >= 0) {
                    tdist = recombfrac[forePointers[j]][i] + recombfrac[i][j];
                    if (tdist < mindist) {
                      mindist = tdist;
                      savedjf = forePointers[j];
                      savedja = j;
                    }
                    tdist = recombfrac[i][j] + recombfrac[i][aftPointers[j]];
                    if (tdist < mindist) {
                      mindist = tdist;
                      savedjf = j;
                      savedja = aftPointers[j];
                    }
                  }
                  else if (forePointers[j] >= 0 && aftPointers[j] == -1) { //j is an aft end.
                    tdist = 2 * recombfrac[i][j]; //Always compare sum of two recombination fractions involving i to avoid biased placement to end.
                    if (tdist < mindist) {
                      mindist = tdist;
                      savedjf = j;
                      savedja = -1;
                    }
                    tdist = recombfrac[i][j] + recombfrac[i][forePointers[j]];
                    if (tdist < mindist) {
                      mindist = tdist;
                      savedjf = forePointers[j];
                      savedja = j;
                    }
                  }
                  else if (forePointers[j] == -1 && aftPointers[j] >= 0) { //j is a fore end.
                    tdist = 2 * recombfrac[i][j];
                    if (tdist < mindist) {
                      mindist = tdist;
                      savedjf = -1;
                      savedja = j;
                    }
                    tdist = recombfrac[i][j] + recombfrac[i][aftPointers[j]];
                    if (tdist < mindist) {
                      mindist = tdist;
                      savedjf = j;
                      savedja = aftPointers[j];
                    }
                  }
                }
                else { //Insert i next to nearest neighbor j.
                  if (recombfrac[i][j] < mindist) {
                    mindist = recombfrac[i][j];
                    if (forePointers[j] >= 0) {
                      if (aftPointers[j] == -1) { //j is an aft end.
                        if (recombfrac[forePointers[j]][j] <= recombfrac[i][forePointers[j]]) { //Marker i goes on the end.
                          savedjf = j;
                          savedja = -1;
                        }
                        else { //Marker i is penultimate.
                          savedjf = forePointers[j];
                          savedja = j;
                        }
                      }
                      else { //j lies within a linkage group.
                        //Take care of the violating case made possible by genotyping errors.
                        if (recombfrac[i][forePointers[j]] < recombfrac[j][forePointers[j]] && recombfrac[i][aftPointers[j]] < recombfrac[j][aftPointers[j]]) {
                          if (recombfrac[i][forePointers[j]] < recombfrac[i][aftPointers[j]]) {
                            savedjf = forePointers[j];
                            savedja = j;
                          }
                          else {
                            savedjf = j;
                            savedja = aftPointers[j];
                          }
                        }
                        else { //Here i follows the normal rules of proximity.
                          if (recombfrac[i][forePointers[j]] < recombfrac[j][forePointers[j]]) {
                            savedjf = forePointers[j];
                            savedja = j;
                          }
                          else {
                            savedjf = j;
                            savedja = aftPointers[j];
                          }
                        }
                        /*This is the way before 10-16-2016.
                        if (recombfrac[i][forePointers[j]] < recombfrac[i][aftPointers[j]]) {
                          savedjf = forePointers[j];
                          savedja = j;
                        }
                        else {
                          savedjf = j;
                          savedja = aftPointers[j];
                        }*/
                      }
                    }
                    else { //j is a fore end.
                      if (recombfrac[i][aftPointers[j]] < recombfrac[j][aftPointers[j]]) {
                        savedjf = j;
                        savedja = aftPointers[j];
                      }
                      else {
                        savedjf = -1;
                        savedja = j;
                      }
                    }
                  }
                }
              }
              if (mindist > maxrf) continue; //Keep ill-fitting markers off the map.
              if (mindist < 0.00000001) continue; //Zero-recombining markers were inserted uncontestedly above.
              foredist = 0; //The following lines assume that the target linkage group is at least six markers long.
              if (savedjf >= 0) {
                g = savedjf;
                foredist = 0;
                while (forePointers[g] >= 0 && foredist < 11) {
                  foredist++;
                  g = forePointers[g];
                }
              }
              aftdist = 0;
              if (savedja >= 0) {
                g = savedja;
                while (aftPointers[g] >= 0 && aftdist < 11) {
                  aftdist++;
                  g = aftPointers[g];
                }
              }
              if (foredist < 5) {
                for (ir = 0; ir < 10; ir++) mergearray[ir] = -1;
                mergearray[0] = i;
                if (savedjf >= 0) g = savedjf;
                else g = savedja;
                while (forePointers[g] >= 0) g = forePointers[g];
                for (ir = 1; ir <= 5; ir++) {
                  mergearray[ir] = g;
                  g = aftPointers[g];
                }
                flankaft = g;
                minmaplength = 2000000000;
                for (ia = 0; ia < nf[6]; ia++) {
                  localmaplength = 0;
                  for (ib = 1; ib < 6; ib++) {
                    localmaplength += recombfrac[mergearray[coeffs[6][ia][ib-1]]][mergearray[coeffs[6][ia][ib]]];
                  }
                  localmaplength += recombfrac[mergearray[coeffs[6][ia][5]]][flankaft];
                  if (localmaplength < minmaplength) {
                    minmaplength = localmaplength;
                    for (ib = 0; ib < 6; ib++) savedv[ib] = mergearray[coeffs[6][ia][ib]];
                  }
                }
                forePointers[savedv[0]] = -1;
                for (ia = 1; ia < 6; ia++) {
                  aftPointers[savedv[ia-1]] = savedv[ia];
                  forePointers[savedv[ia]] = savedv[ia-1];
                }
                aftPointers[savedv[5]] = flankaft;
                forePointers[flankaft] = savedv[5];
              }
              else if (aftdist < 5) {
                for (ir = 0; ir < 10; ir++) mergearray[ir] = -1;
                mergearray[0] = i;
                if (savedja >= 0) g = savedja;
                else g = savedjf;
                while (aftPointers[g] >= 0) g = aftPointers[g];
                for (ir = 1; ir <= 5; ir++) {
                  mergearray[ir] = g;
                  g = forePointers[g];
                }
                flankfore = g;
                minmaplength = 2000000000;
                for (ia = 0; ia < nf[6]; ia++) {
                  localmaplength = 0;
                  for (ib = 1; ib < 6; ib++) {
                    localmaplength += recombfrac[mergearray[coeffs[6][ia][ib-1]]][mergearray[coeffs[6][ia][ib]]];
                  }
                  localmaplength += recombfrac[mergearray[coeffs[6][ia][0]]][flankfore];
                  if (localmaplength < minmaplength) {
                    minmaplength = localmaplength;
                    for (ib = 0; ib < 6; ib++) savedv[ib] = mergearray[coeffs[6][ia][ib]];
                  }
                }
                aftPointers[flankfore] = savedv[0];
                forePointers[savedv[0]] = flankfore;
                for (ia = 1; ia < 6; ia++) {
                  aftPointers[savedv[ia-1]] = savedv[ia];
                  forePointers[savedv[ia]] = savedv[ia-1];
                }
                aftPointers[savedv[5]] = -1;
              }
              else { /*savedjf and savedja lie deeper within the linkage group. Insert i with local permutation.*/
                for (ir = 0; ir < 10; ir++) mergearray[ir] = -1;
                mergearray[0] = i;
                g = savedjf;
                for (ir = 1; ir < 6; ir++) {
                  mergearray[ir] = g;
                  g = forePointers[g];
                }
                flankfore = g;
                g = savedja;
                for (ir = 6; ir < 10; ir++) {
                  mergearray[ir] = g;
                  g = aftPointers[g];
                }
                flankaft = g;
                minmaplength = 2000000000;
                for (ia = 0; ia < nf[10]; ia++) {
                  localmaplength = 0;
                  for (ib = 1; ib < 10; ib++) {
                    localmaplength += recombfrac[mergearray[coeffs[10][ia][ib-1]]][mergearray[coeffs[10][ia][ib]]];
                  }
                  if (flankfore >= 0) localmaplength += recombfrac[mergearray[coeffs[10][ia][0]]][flankfore];
                  if (flankaft >= 0) localmaplength += recombfrac[mergearray[coeffs[10][ia][9]]][flankaft];
                  if (localmaplength < minmaplength) {
                    minmaplength = localmaplength;
                    for (ib = 0; ib < 10; ib++) savedv[ib] = mergearray[coeffs[10][ia][ib]];
                  }
                }
                aftPointers[flankfore] = savedv[0];
                forePointers[savedv[0]] = flankfore;
                for (ia = 1; ia < 10; ia++) {
                  aftPointers[savedv[ia-1]] = savedv[ia];
                  forePointers[savedv[ia]] = savedv[ia-1];
                }
                aftPointers[savedv[9]] = flankaft;
                forePointers[flankaft] = savedv[9]; //lgStarts should be unaffected.
              }
            }
          }
        }
        //Cleave any joins that exceed maxrf after insertion with local permutation.
        for (i = 0; i < nLoci; i++) {
          if (aftPointers[i] >= 0 && recombfrac[i][aftPointers[i]] > maxrf) {
            forePointers[aftPointers[i]] = -1;
            aftPointers[i] = -1;
          }
        }
        //Recalculate lengths and counts of linkage groups.  The rankings of counts could have changed from the original map.
        nLinkageGroups = 0; //Now update nLinkageGroups and allocate all the new arrays.
        for (i = 0; i < nLoci; i++) {
          if (forePointers[i] == -1 && aftPointers[i] >= 0) nLinkageGroups++;
        }
        lgnewLengths = (double *) calloc(nLinkageGroups, sizeof(double));
        lgnewCounts = (int *) calloc(nLinkageGroups, sizeof(int));
        lgnewStarts = (int *) calloc(nLinkageGroups, sizeof(int));
        lgnewEnds = (int *) calloc(nLinkageGroups, sizeof(int));
        sortedlgnewCounts = (int *) calloc(nLinkageGroups, sizeof(int));
        donenewLinkageGroups = (int *) calloc(nLinkageGroups + 1, sizeof(int));
        for (i = 0; i < nLinkageGroups; i++) {
          lgnewStarts[i] = -1;
          lgnewEnds[i] = -1;
        }
        k = 0;
        for (i = 0; i < nLoci; i++) {
          if (forePointers[i] == -1 && aftPointers[i] >= 0) {
            lgnewStarts[k] = i;
            g = i;
            while (aftPointers[g] >= 0) { //This has to be this way to collect recombination fractions.
              lgnewCounts[k]++;
              lgnewLengths[k] += recombfrac[g][aftPointers[g]];
              lgnewEnds[k] = aftPointers[g];
              g = aftPointers[g];
            }
            lgnewCounts[k]++; //to include the last locus, whose aftPointer is -1.
            k++;
          }
        }
        for (i = 0; i < nLinkageGroups; i++) {
          sortedlgnewCounts[i] = lgnewCounts[i];
          donenewLinkageGroups[i] = 0;
        }
        qsort(sortedlgnewCounts, nLinkageGroups, sizeof(int), compdown);
        totlength = 0;
        totcount = 0;
        ic = 0;
        fprintf(OUT, "\nSummary of final map:\n");
        for (i = 0; i < nLinkageGroups; i++) {
          for (j = 0; j < nLinkageGroups; j++) {
            if (lgnewCounts[j] == sortedlgnewCounts[i] && donenewLinkageGroups[j] == 0) {
              fprintf(OUT, "linkage group %d\tstart = %d\tend = %d\tcount = %d\tlength = %f\n", ic, lgnewStarts[j], lgnewEnds[j], lgnewCounts[j], lgnewLengths[j]);
              donenewLinkageGroups[j] = 1;
              totlength += lgnewLengths[j];
              totcount += lgnewCounts[j];
              ic++;
            }
          }
        }
        if (totcount > 0) lengthpercount = totlength / totcount;
        else lengthpercount = 0;
        fprintf(OUT, "Total length of revised map = %.6f count = %d mean = %.6f\n", totlength, totcount, lengthpercount);
        for (i = 0; i < nLinkageGroups; i++) donenewLinkageGroups[i] = 0; //Reload to print out the map itself.
        for (i = 0; i < nLoci; i++) membership[i] = -1; //Reload memberships.
        fprintf(OUT, "\nFinal map with forced placement of frozen-out loci:\n");
        for (i = 0; i < nLinkageGroups; i++) {
          for (j = 0; j < nLinkageGroups; j++) {
            if (lgnewCounts[j] == sortedlgnewCounts[i] && donenewLinkageGroups[j] == 0) {
              fprintf(OUT, "linkage group %d\tstart = %d\tend \t%d\tcount = %d\tlength = %f\n", i, lgnewStarts[j], lgnewEnds[j], lgnewCounts[j], lgnewLengths[j]);
              cudist = 0.0;
              g = lgnewStarts[j];
              while (aftPointers[g] >= 0) {
                fprintf(OUT, "%s\t%d\t%d\t%f\t%f\t%f\t%s\t%f", names[g], g, aftPointers[g], recombfrac[g][aftPointers[g]], cudist, minofrfs[g], names[nearestneighbors[g]], segdises[g]);
                if (measflag == 1) fprintf(OUT, "\t%f", hightolows[g]);
                fprintf(OUT, "\n");
                membership[g] = i;
                cudist += recombfrac[g][aftPointers[g]];
                g = aftPointers[g];
              }
              if (queryflag == 1 && querycodes[g] == 1) {
                mapqueries[mq] = g;
                mq++;
              }
              fprintf(OUT, "%s\t%d\t%d\t-\t%f\t%f\t%s\t%f", names[g], g, aftPointers[g], cudist, minofrfs[g], names[nearestneighbors[g]], segdises[g]);
              if (measflag == 1) fprintf(OUT, "\t%f", hightolows[g]);
              fprintf(OUT, "\n");
              membership[g] = i;
              donenewLinkageGroups[j] = 1;
              break;
            }
          }
        }
        fprintf(OUT, "\n"); //blank line for plotorders.pl to read the output file
      }
    }
    if (tryflag == 1) {
      FLG = fopen(flagfile, "r");
      flaggedname = (char *) calloc(maxnamelength + 4, sizeof(char));
      flaglines = 0;
      while (fgets(line, maxnamelength, FLG)) flaglines++;
      flaggednames = (char **) malloc((flaglines + 1) * sizeof(char *));
      flaggedindices = (int *) malloc((flaglines + 1) * sizeof(int));
      for (i = 0; i <= flaglines; i++) {
        flaggednames[i] = (char *) calloc(maxnamelength + 4, sizeof(char));
        flaggedindices[i] = -1;
      }
      fseek(FLG, 0, 0);
      ig = 0;
      while (fgets(line, maxnamelength, FLG)) { //Read in the flagged markers.
        flaggedname = strtok(line, " ,|\t\n"); //Emulate a Perl chomp.
        strcpy(flaggednames[ig], flaggedname);
        for (i = 0; i < nLoci; i++) {
          if (strcmp(flaggedname, names[i]) == 0) {
            flaggedindices[ig] = i;
            break;
          }
        }
        ig++;
      }
      fprintf(OUT, "\n");
      for (i = 0; i < flaglines; i++) fprintf(OUT, "flagged name = %s its index = %d\n", flaggednames[i], flaggedindices[i]);
      for (i = 0; i < flaglines; i++) { //Cut the flagged markers out of the map.
        if (forePointers[flaggedindices[i]] >= 0) {
          if (aftPointers[flaggedindices[i]] >= 0) {
            aftPointers[forePointers[flaggedindices[i]]] = aftPointers[flaggedindices[i]];
            forePointers[aftPointers[flaggedindices[i]]] = forePointers[flaggedindices[i]];
            forePointers[flaggedindices[i]] = -1;
            aftPointers[flaggedindices[i]] = -1;
          }
          else {
            aftPointers[forePointers[flaggedindices[i]]] = -1;
            forePointers[flaggedindices[i]] = -1;
          }
        }
        else {
          if (aftPointers[flaggedindices[i]] >= 0) {
            forePointers[aftPointers[flaggedindices[i]]] = -1;
            aftPointers[flaggedindices[i]] = -1;
          }
          else fprintf(OUT, "Marker %d has already been removed from the map.\n", flaggedindices[i]);
        }
      }
      for (i = 0; i < flaglines; i++) fprintf(OUT, "After removing flagged marker %s from the map, foreP = %d marker = %d aftP = %d\n", names[flaggedindices[i]], forePointers[flaggedindices[i]], flaggedindices[i], aftPointers[flaggedindices[i]]);
      for (k = 0; k < nLoci; k++) {
        if (recombfrac[k][31742] < 0.18) fprintf(OUT, "diag at line 3163: k = %d names[k] = %s recombfrac[%d][31742] = %f\n", k, names[k], k, recombfrac[k][31742]);

      }
      for (k = 0; k < flaglines; k++) {
        mindist = 1.0;
        for (i = 0; i < flaglines; i++) { //Find the flagged and non-flagged loci that are closest together in the map.
          if (forePointers[flaggedindices[i]] >= 0 || aftPointers[flaggedindices[i]] >= 0) continue;
          for (j = 0; j < nLoci; j++) { //
            if (j == flaggedindices[i]) continue; //Do not join to self.  Actually, the next line should take care of this.
            if (forePointers[j] == -1 && aftPointers[j] == -1) continue; //Do not join to an unbound locus.
            if (fracUsedCode[j] <= 1 || minofrfs[j] > maxminrf || segdises[j] > segdisthreshold) continue;
            if (forePointers[j] >= 0) {
              tdist = recombfrac[forePointers[j]][flaggedindices[i]] + recombfrac[flaggedindices[i]][j];
              if (tdist < mindist) { //Do this regardless of whether j is internal or an aft end.
                mindist = tdist;
                savedjf = forePointers[j];
                savedja = j;
                savedi = flaggedindices[i];
              }
              if (aftPointers[j] == -1) { //j is an aft end; test after it.
                tdist = recombfrac[flaggedindices[i]][j] + recombfrac[forePointers[j]][j];
                if (tdist < mindist) {
                  mindist = tdist;
                  savedjf = j;
                  savedja = -1;
                  savedi = flaggedindices[i];
                }
              }
            }
            else if (aftPointers[j] >= 0) { //j is a fore end
              tdist = recombfrac[flaggedindices[i]][j] + recombfrac[j][aftPointers[j]];
              if (tdist < mindist) {
                mindist = tdist;
                savedjf = -1;
                savedja = j;
                savedi = flaggedindices[i];
              }
            }
          }
        }
        fprintf(OUT, "About to insert marker %s savedi = %d savedjf = %d savedja = %d\n", names[savedi], savedi, savedjf, savedja);
        foredist = 0; //The following lines assume that the target linkage group is at least six markers long.
        if (savedjf >= 0) {
          g = savedjf;
          foredist = 0;
          while (forePointers[g] >= 0 && foredist < 11) {
            foredist++;
            g = forePointers[g];
          }
        }
        aftdist = 0;
        if (savedja >= 0) {
          g = savedja;
          while (aftPointers[g] >= 0 && aftdist < 11) {
            aftdist++;
            g = aftPointers[g];
          }
        }
        if (foredist < 5) {
          for (ir = 0; ir < 10; ir++) mergearray[ir] = -1;
          mergearray[0] = savedi;
          if (savedjf >= 0) g = savedjf;
          else g = savedja;
          while (forePointers[g] >= 0) g = forePointers[g];
          for (ir = 1; ir <= 5; ir++) {
            mergearray[ir] = g;
            g = aftPointers[g];
          }
          flankaft = g;
          minmaplength = 2000000000;
          for (ia = 0; ia < nf[6]; ia++) {
            localmaplength = 0;
            for (ib = 1; ib < 6; ib++) {
              localmaplength += recombfrac[mergearray[coeffs[6][ia][ib-1]]][mergearray[coeffs[6][ia][ib]]];
            }
            localmaplength += recombfrac[mergearray[coeffs[6][ia][5]]][flankaft];
            if (localmaplength < minmaplength) {
              minmaplength = localmaplength;
              for (ib = 0; ib < 6; ib++) savedv[ib] = mergearray[coeffs[6][ia][ib]];
            }
          }
          forePointers[savedv[0]] = -1;
          for (ia = 1; ia < 6; ia++) {
            aftPointers[savedv[ia-1]] = savedv[ia];
            forePointers[savedv[ia]] = savedv[ia-1];
          }
          aftPointers[savedv[5]] = flankaft;
          forePointers[flankaft] = savedv[5];
        }
        else if (aftdist < 5) {
          for (ir = 0; ir < 10; ir++) mergearray[ir] = -1;
          mergearray[0] = savedi;
          if (savedja >= 0) g = savedja;
          else g = savedjf;
          while (aftPointers[g] >= 0) g = aftPointers[g];
          for (ir = 1; ir <= 5; ir++) {
            mergearray[ir] = g;
            g = forePointers[g];
          }
          flankfore = g;
          minmaplength = 2000000000;
          for (ia = 0; ia < nf[6]; ia++) {
            localmaplength = 0;
            for (ib = 1; ib < 6; ib++) {
              localmaplength += recombfrac[mergearray[coeffs[6][ia][ib-1]]][mergearray[coeffs[6][ia][ib]]];
            }
            localmaplength += recombfrac[mergearray[coeffs[6][ia][0]]][flankfore];
            if (localmaplength < minmaplength) {
              minmaplength = localmaplength;
              for (ib = 0; ib < 6; ib++) savedv[ib] = mergearray[coeffs[6][ia][ib]];
            }
          }
          aftPointers[flankfore] = savedv[0];
          forePointers[savedv[0]] = flankfore;
          for (ia = 1; ia < 6; ia++) {
            aftPointers[savedv[ia-1]] = savedv[ia];
            forePointers[savedv[ia]] = savedv[ia-1];
          }
          aftPointers[savedv[5]] = -1;
        }
        else { /*savedjf and savedja lie deeper within the linkage group. Insert i with local permutation.*/
          for (ir = 0; ir < 10; ir++) mergearray[ir] = -1;
          mergearray[0] = savedi;
          g = savedjf;
          for (ir = 1; ir < 6; ir++) {
            mergearray[ir] = g;
            g = forePointers[g];
          }
          flankfore = g;
          g = savedja;
          for (ir = 6; ir < 10; ir++) {
            mergearray[ir] = g;
            g = aftPointers[g];
          }
          flankaft = g;
          minmaplength = 2000000000;
          for (ia = 0; ia < nf[10]; ia++) {
            localmaplength = 0;
            for (ib = 1; ib < 10; ib++) {
              localmaplength += recombfrac[mergearray[coeffs[10][ia][ib-1]]][mergearray[coeffs[10][ia][ib]]];
            }
            if (flankfore >= 0) localmaplength += recombfrac[mergearray[coeffs[10][ia][0]]][flankfore];
            if (flankaft >= 0) localmaplength += recombfrac[mergearray[coeffs[10][ia][9]]][flankaft];
            if (localmaplength < minmaplength) {
              minmaplength = localmaplength;
              for (ib = 0; ib < 10; ib++) savedv[ib] = mergearray[coeffs[10][ia][ib]];
            }
          }
          aftPointers[flankfore] = savedv[0];
          forePointers[savedv[0]] = flankfore;
          for (ia = 1; ia < 10; ia++) {
            aftPointers[savedv[ia-1]] = savedv[ia];
            forePointers[savedv[ia]] = savedv[ia-1];
          }
          aftPointers[savedv[9]] = flankaft;
          forePointers[flankaft] = savedv[9]; //lgStarts should be unaffected.
        }
      }
      //Now recalculate the map length and output the map again.  See lines 3029 through 3081.
      nLinkageGroups = 0; //Now update nLinkageGroups and reinitialize all the new arrays.
      for (i = 0; i < nLoci; i++) {
        if (forePointers[i] == -1 && aftPointers[i] >= 0) nLinkageGroups++;
      }
      for (i = 0; i < nLinkageGroups; i++) {
        lgnewLengths[i] = 0;
        lgnewCounts[i] = 0;
        lgnewStarts[i] = -1;
        lgnewEnds[i] = -1;
        sortedlgnewCounts[i] = 0;
      }
      k = 0;
      for (i = 0; i < nLoci; i++) {
        if (forePointers[i] == -1 && aftPointers[i] >= 0) {
          lgnewStarts[k] = i;
          g = i;
          while (aftPointers[g] >= 0) { //This has to be this way to collect recombination fractions.
            lgnewCounts[k]++;
            lgnewLengths[k] += recombfrac[g][aftPointers[g]];
            lgnewEnds[k] = aftPointers[g];
            g = aftPointers[g];
          }
          lgnewCounts[k]++; //to include the last locus, whose aftPointer is -1.
          k++;
        }
      }
      for (i = 0; i < nLinkageGroups; i++) {
        sortedlgnewCounts[i] = lgnewCounts[i];
        donenewLinkageGroups[i] = 0;
      }
      qsort(sortedlgnewCounts, nLinkageGroups, sizeof(int), compdown);
      totlength = 0;
      totcount = 0;
      ic = 0;
      fprintf(OUT, "\nSummary of map after moving flagged loci:\n");
      for (i = 0; i < nLinkageGroups; i++) {
        for (j = 0; j < nLinkageGroups; j++) {
          if (lgnewCounts[j] == sortedlgnewCounts[i] && donenewLinkageGroups[j] == 0) {
            fprintf(OUT, "linkage group %d\tstart = %d\tend = %d\tcount = %d\tlength = %f\n", ic, lgnewStarts[j], lgnewEnds[j], lgnewCounts[j], lgnewLengths[j]);
            donenewLinkageGroups[j] = 1;
            totlength += lgnewLengths[j];
            totcount += lgnewCounts[j];
            ic++;
          }
        }
      }
      if (totcount > 0) lengthpercount = totlength / totcount;
      else lengthpercount = 0;
      fprintf(OUT, "Total length of post-flagged map = %.6f count = %d mean = %.6f\n", totlength, totcount, lengthpercount);
      for (i = 0; i < nLinkageGroups; i++) donenewLinkageGroups[i] = 0; //Reload to print out the map itself.
      for (i = 0; i < nLoci; i++) membership[i] = -1; //Reload memberships.
      fprintf(OUT, "\nMap after moving flagged, user-selected loci:\n");
      for (i = 0; i < nLinkageGroups; i++) {
        for (j = 0; j < nLinkageGroups; j++) {
          if (lgnewCounts[j] == sortedlgnewCounts[i] && donenewLinkageGroups[j] == 0) {
            fprintf(OUT, "linkage group %d\tstart = %d\tend \t%d\tcount = %d\tlength = %f\n", i, lgnewStarts[j], lgnewEnds[j], lgnewCounts[j], lgnewLengths[j]);
            cudist = 0.0;
            g = lgnewStarts[j];
            while (aftPointers[g] >= 0) {
              fprintf(OUT, "%s\t%d\t%d\t%f\t%f\t%f\t%s\t%f", names[g], g, aftPointers[g], recombfrac[g][aftPointers[g]], cudist, minofrfs[g], names[nearestneighbors[g]], segdises[g]);
              if (measflag == 1) fprintf(OUT, "\t%f", hightolows[g]);
              fprintf(OUT, "\n");
              membership[g] = i;
              cudist += recombfrac[g][aftPointers[g]];
              g = aftPointers[g];
            }
            if (queryflag == 1 && querycodes[g] == 1) {
              mapqueries[mq] = g;
              mq++;
            }
            fprintf(OUT, "%s\t%d\t%d\t-\t%f\t%f\t%s\t%f", names[g], g, aftPointers[g], cudist, minofrfs[g], names[nearestneighbors[g]], segdises[g]);
            if (measflag == 1) fprintf(OUT, "\t%f", hightolows[g]);
            fprintf(OUT, "\n");
            membership[g] = i;
            donenewLinkageGroups[j] = 1;
            break;
          }
        }
      }
      fprintf(OUT, "\n"); //blank line for plotorders.pl to read the output file
    } //End of tryflag loop.
  }
  //These free statements become necessary if the content of main is ever called in a subroutine.  As the program
  //stands now, all memory is freed without these free statements when the program finishes.
  free(forePointers); free(aftPointers);
  free(coeffs[0]); free(coeffs[1]); free(coeffs[2]);
  for (i = 3; i <= 10; i++) {
    for (j = 0; j < nf[i]; j++) free(coeffs[i][j]);
    free(coeffs[i]);
  }
  free(coeffs);
  for (i = 0; i < nLoci; i++) {
    free(recombfrac[i]); free(nMissing[i]); free(names[i]);
  }
  free(recombfrac); free(nMissing); free(tempname); free(names);
  if (continuation >= 1) {free(lgnewStarts); free(lgnewEnds); free(lgnewCounts); free(lgnewLengths); free(sortedlgCounts);}
  free(doneLinkageGroups); free(membership); free(donenewLinkageGroups);
  for (i = 0; i < orignlgroups; i++) free(minlgdists[i]);
  free(minlgdists); free(printstring); free(cfgline); free(datafile); free(specialfile); free(prdfile); free(outfile);
  free(ngrfile); free(outafile); free(queryfile); free(measfile); free(markernamesfile); free(rfracfile); free(usedcodefile);
  free(segdisfile); 
  free(dedupoutfile); 
  free(line);
  //free(token);
  free(lctoken); 
  free(preDosage); 
  free(postDosage);
  free(plantsUsed); free(insertions); free(mfores); free(mafts);
  if (measflag == 1) free(hightolows);
  if (continuation > 0) free(pool);
  currtime = clock();
  duration = (double) (currtime - timeused2) / CLOCKS_PER_SEC;
  fprintf(OUT, "done with program, time used = %f\n", duration);
  return(0);
}
 
int checkForPresence(int nLoci, int k, int m, int *foreP, int *aftP, FILE *OUT) {
  int found, found1, found2, found3, iv, g;
  /*Check if k and m are in the same linkage group from both directions for both groups.*/
  found = 0; //Check that the two loci are not already in the same linkage group; prevent cycles.
  found1 = 0;
  found2 = 0; //These three founds are temporary, to make sure that all joins are consistent with one another.
  found3 = 0;
  g = k;
  iv = 0;
  while (foreP[g] >= 0) {
    if (iv > nLoci) {
      fprintf(OUT, "iv = %d, in checkForPresence line 2649, having completed a cycle in the map at locus %d.\n", iv, g);
      return(4);
    }
    iv++;
    g = foreP[g];
  }
  iv = 0;
  fprintf(OUT, "subdiag a: g = %d\n", g);
  while (aftP[g] >= 0) {
    if (iv > nLoci) {
      fprintf(OUT, "iv = %d, in checkForPresence line 2659, having completed a cycle in the map at locus %d.\n", iv, g);
      return(4);
    }
    if (g == m) found = 1;
    iv++;
    g = aftP[g];
  }
  if (g == m) found = 1;
  iv = 0;
  fprintf(OUT, "subdiag b: g = %d\n", g);
  while (g >= 0) {
    if (iv > nLoci) {
      fprintf(OUT, "checkForPresence line 2671, iv = %d, having completed a cycle in the map at locus %d.\n", iv, g);
      return(4);
    }
    if (g == m) found1 = 1;
    iv++;
    g = foreP[g];
  }
  iv = 0;
  g = m;
  while (foreP[g] >= 0) {
    if (iv > nLoci) {
      fprintf(OUT, "in checkForPresence, line 2682, iv = %d, having completed a cycle in the map at locus %d.\n", iv, g);
      return(4);
    }
    iv++;
    g = foreP[g];
  }
  fprintf(OUT, "subdiag c: g = %d\n", g);
  iv = 0;
  while (aftP[g] >= 0) {
    if (iv > nLoci) {
      fprintf(OUT, "in checkForPresence, line 2692, iv = %d, having completed a cycle in the map at locus %d.\n", iv, g);
      return(4);
    }
    if (g == k) found2 = 1;
    iv++;
    g = aftP[g];
  }
  if (g == k) found2 = 1;
  iv = 0;
  fprintf(OUT, "subdiag d: g = %d\n", g);
  while (g >= 0) {
    if (iv > nLoci) {
      fprintf(OUT, "in checkForPresence, line 2704, iv = %d, having completed a cycle in the map at locus %d.\n", iv, g);
      return(4);
    }
    if (g == k) found3 = 1;
    iv++;
    g = foreP[g];
  }
  if (found != found1 || found != found2 || found != found3) {
    fprintf(OUT, "forePointers and aftPointers arrays are inconsistent for k = %d and m = %d: found = %d found1 = %d found2 = %d found3 = %d\n", k, m, found, found1, found2, found3);
    return(1);
  }
  if (found > 0 || found1 > 0 || found2 > 0 || found3 > 0) return(1);
  else return(0);
}

int dcomp(const void *r, const void *s) {
  if (*(double*)r > *(double*)s) return(1);
  else if (*(double*)r == *(double*)s) return(0);
  else if (*(double*)r < *(double*)s) return(-1);
  else return(-8);
}

int comp(const void *r, const void *s) {
  if (*((int*)r) > *((int*)s)) return(1);
  else if (*((int*)r) == *((int*)s)) return(0);
  else if (*((int*)r) < *((int*)s)) return(-1);
  else return(-8);
}

int compdown(const void *r, const void *s) {
  if (*((int*)r) > *((int*)s)) return(-1);
  else if (*((int*)r) == *((int*)s)) return(0);
  else if (*((int*)r) < *((int*)s)) return(1);
  else return(8); //This impossible condition silences a compiler warning.
}

int getmaxrf(int *foreP, double **rf, int n, double mxrf, double tp, FILE *OUT) {
  int i, svi, svj;
  double xrf;
  xrf = 0;
  for (i = 0; i < n; i++) {
    if (foreP[i] >= 0) {
      if (rf[i][foreP[i]] > xrf) {
        xrf = rf[i][foreP[i]];
        svi = i;
        svj = foreP[i];
      }
    }
  }
  if (xrf > mxrf || xrf > 3 * tp) {
    fprintf(OUT, "Maximum rf %.4f was exceeded for loci %d and %d with rf = %.4f\n", mxrf, svj, svi, xrf);
    return(1);
  }
  else return(0);
}

void printStats(char* printstring, int nLoci, double **rf, int *foreP, int *aftP, FILE *OUT) {
  int i, g, starts, ends, count;
  double maplength;
  starts = 0; ends = 0; count = 0;
  maplength = 0.0;
  for (i = 0; i < nLoci; i++) {
    if (foreP[i] == -1 && aftP[i] >= 0) {
      starts++;
      g = i;
      while (g >= 0) {
        count++;
        g = aftP[g];
      }
      g = i;
      while (aftP[g] >= 0) {
        maplength += rf[g][aftP[g]];
        g = aftP[g];
      }
    }
    if (foreP[i] >= 0 && aftP[i] == -1) ends++;
  }
  fprintf(OUT, "Printed statistics for %s starts = %d ends = %d map contains %d loci map length = %.6f\n", printstring, starts, ends, count, maplength);
}

void printtestloci(int m, int n, int window, double **rf, int *foreP, int *aftP, FILE *OUT) {
  int i, g, found;
  g = m;
  for (i = 0; i < window; i++) {
    if (foreP[g] >= 0) g = foreP[g];
    else break;
  }
  fprintf(OUT, "DIAG for test locus %d on the folowing lines:\n", m);
  found = 0;
  for (i = 0; i < 2 * window; i++) {
    if (g == n) found = 1;
    fprintf(OUT, "%d\t%d\t%.6f\n", g, aftP[g], rf[g][aftP[g]]);
    if (aftP[g] >= 0) g = aftP[g];
    else break;
  }
  if (found == 0) {
    g = n;
    for (i = 0; i < window; i++) {
      if (foreP[g] >= 0) g = foreP[g];
      else break;
    }
    fprintf(OUT, "DIAG for test locus %d on the folowing lines:\n", n);
    for (i = 0; i < 2 * window; i++) {
      fprintf(OUT, "%d\t%d\t%.6f\n", g, aftP[g], rf[g][aftP[g]]);
      if (aftP[g] >= 0) g = aftP[g];
      else break;
    }
  }
  fprintf(OUT, "\n");
}

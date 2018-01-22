#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <cmath>
#include <errno.h>
#include <regex.h>

/*This program splits its input ABH genotypes file into separate linkage groups on the basis of a preliminary subset map.*/
/* function prototypes:
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
*/

int main(int argc, char *argv[]) {
  char *line, *token, *cfgline, *cfgfile, *inmapfile, *datafile, *filler, *value, *varname, *keyword, **refnames, *datname, **refchdata;
  char *outfilename, *outfilestem, *header1, *header2, *header3, *savedname, *pos;
  int i, j, k, lg, rval, cllength, ierror, maxmaplinelength, maxdatlinelength, regresultkey, regresultlg, getflag, currlength;
  int refnLoci, *reflgs, nlgs, maxrefnamelength, maxoutfilenamelength, plantsTotal, mindiff, currdiff, savedgroup, maxdatnamelength;
  double minrf;
  FILE *CFG, *INM, *DAT, **OUD;
  clock_t timeused, timeused1, timeused2, currtime;
  regex_t patternkey, patternlg, patternast;
  if(strcmp(argv[1], "-") == 0 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-help") == 0) {
    printf("Usage: %s configuration_file_name [configuration_file_line_length]\n", argv[0]);
    printf("The configuration file name is mandatory and the configuration file line length is optional.\n");
    printf("Specify the configuration file line length only if it needs to exceed 512 characters.\n");
    printf("Comments in the configuration file begin with # or // and can follow data if separated by one or more spaces.\n");
    printf("Data lines can appear in any order.  Filenames need not be quoted and can include the full path.\n");
    printf("Data lines consist of space-delimited triplets of attribute, \"=\", and value.\n");
    printf("We need to finish this block.\n");
  }
  cllength = 512;
  if (argc == 3) {
    k = (int)strtol(argv[2], (char**)NULL, 10);
    if (k > cllength) {
      cllength = k;
      printf("Maximum file name length cllength has been increased from 512 to %d\n", cllength);
    }
  }
  maxdatlinelength = 0;
  maxmaplinelength = 0;
  maxoutfilenamelength = 0;
  cfgline = (char *) calloc(cllength, sizeof(char));
  datafile = (char *) calloc(cllength + 2, sizeof(char));
  inmapfile = (char *) calloc(cllength + 2, sizeof(char));
  filler = (char *) calloc(cllength, sizeof(char));
  varname = (char *) calloc(cllength, sizeof(char));
  keyword = (char *) calloc(cllength, sizeof(char));
  CFG = fopen(argv[1], "r");
  printf("The configuration file has been opened successfully.\n");
  fflush(NULL);
  while (fgets(cfgline, cllength, CFG)) {
    ierror = 1;
    varname = strtok(cfgline, " ,|\t\n");
    filler = strtok(NULL, " ,|\t\n");
    value = strtok(NULL, " ,|\t\n");
    if (varname[0] == '#') {
      ierror = 0; //Do nothing; it's a comment.
      continue;
    }
    if (varname[0] == '/' && varname[1] == '/') {
      ierror = 0; //This is also a comment.
      continue;
    }
    if (strcmp(varname, "datafile") == 0) {
      strcpy(datafile, value);
      ierror = 0;
    }
    if (strcmp(varname, "inmapfile") == 0) {
      strcpy(inmapfile, value);
      ierror = 0;
    }
    if (strcmp(varname, "keyword") == 0) {
      strcpy(keyword, value);
      ierror = 0;
    }
    if (strcmp(varname, "maxmaplinelength") == 0) {
      maxmaplinelength = (int)strtol(value, (char **)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "maxdatlinelength") == 0) {
      maxdatlinelength = (int)strtol(value, (char **)NULL, 10);
      ierror = 0;
    }
    if (strcmp(varname, "maxoutfilenamelength") == 0) {
      maxoutfilenamelength = (int)strtol(value, (char **)NULL, 10);
      outfilestem = (char *) calloc(maxoutfilenamelength, sizeof(char));
      outfilename = (char *) calloc(maxoutfilenamelength, sizeof(char));
      ierror = 0;
    }
    if (strcmp(varname, "outfilestem") == 0) {
      if (maxoutfilenamelength == 0) {
        printf("The line for maxoutfilenamelength must precede the line for outfilestem in file %s\n", argv[1]);
        return(1);
      }
      strcpy(outfilestem, value);
      ierror = 0;
    }
  }
  printf("The configuration file has been read.\n");
  fflush(NULL);
  fclose(CFG);
  if (strlen(outfilestem) > maxoutfilenamelength - 7) {
    printf("The name length of the proposed output file name %s will exceed the maximum allowed by maxoutfilenamelength = %d\n", outfilestem, maxoutfilenamelength);
    return(1);
  }
  if (maxdatlinelength == 0) {
    printf("No value was given for maxdatlinelength.\n");
    return(1);
  }
  if (maxmaplinelength == 0) {
    printf("No value was given for maxmaplinelength.\n");
    return(1);
  }
  if (maxmaplinelength > maxdatlinelength) {
    line = (char *) calloc(maxmaplinelength, sizeof(char));
    token = (char *) calloc(maxmaplinelength, sizeof(char));
  }
  else {
    line = (char *) calloc(maxdatlinelength, sizeof(char));
    token = (char *) calloc(maxdatlinelength, sizeof(char));
  }
  header1 = (char *) calloc(maxdatlinelength, sizeof(char));
  header2 = (char *) calloc(maxdatlinelength, sizeof(char));
  header3 = (char *) calloc(maxdatlinelength, sizeof(char));
  if (strlen(inmapfile) > 0) INM = fopen(inmapfile, "r");
  else {
    printf("The configuration file did not specify the inmapfile.\n");
    return(1);
  }
  printf("Inmap file %s has been opened or errno = %d.\n", inmapfile, errno);
  fflush(NULL);
  if (INM==NULL) {
    printf("INM is null after all.\n");
    fflush(NULL);
  }
  
  rval = regcomp(&patternkey, keyword, 0);
  if (rval != 0) {
    printf("Compilation of %s as a regular expression failed.\n", keyword);
    return(1);
  }
  rval = regcomp(&patternlg, "linkage group", 0);
  if (rval != 0) {
    printf("Compilation of \"linkage group\" as a regular expression failed.\n");
    return(1);
  }
  rval = regcomp(&patternast, "*", 0);
  if (rval != 0) {
    printf("Compilation of leading asterisk as a regular expression failed.\n");
    return(1);
  }
  getflag = 0;
  nlgs = 0;
  refnLoci = 0;
  maxrefnamelength = 0;
  while (fgets(line, maxmaplinelength, INM)) {
    if (strlen(line) <= 1) getflag = 0;
    if (getflag == 1) {
      if (regexec(&patternlg, line, 0, NULL, 0)) {
        refnLoci++;
        token = strtok(line, " |\t\n,");
        if (strlen(token) > maxrefnamelength) maxrefnamelength = strlen(token);
        //printf("nlgs = %d refnLoci = %d token = %s\n", nlgs, refnLoci, token);
      }
      else nlgs++;
    }
    /*int regexec (regex_t *compiled, char *string, size_t nmatch, regmatch_t matchptr [], int eflags)*/
    if (!regexec(&patternkey, line, 0, NULL, 0)) getflag = 1;
  }
  printf ("Ending value of getflag is %d nlgs = %d refnLoci = %d maxrefnamelength = %d\n", getflag, nlgs, refnLoci, maxrefnamelength);
  fflush(NULL);
  datname = (char *) calloc(maxrefnamelength * 2, sizeof(char)); //Consider changing this limit.
  refnames = (char **) calloc(refnLoci + 5, sizeof(char *));
  for (i = 0; i < refnLoci + 5; i++) refnames[i] = (char *) calloc(maxrefnamelength + 5, sizeof(char));
  savedname = (char *) calloc(maxrefnamelength + 5, sizeof(char));
  reflgs = (int *) calloc(refnLoci + 5, sizeof(int));
  for (i = 0; i < refnLoci + 5; i++) reflgs[i] = -1;
  fseek(INM, 0, SEEK_SET);
  getflag = 0;
  lg = -1;
  j = 0;
  printf("Ready to begin the second pass through the small map.\n");
  fflush(NULL);
  while (fgets(line, maxmaplinelength, INM)) {
    if (strlen(line) <= 1) getflag = 0;
    if (getflag == 1) {
      if (regexec(&patternlg, line, 0, NULL, 0)) {
        token = strtok(line, " |\t\n,");
        strcpy(refnames[j], token);
        reflgs[j] = lg;
        j++;
      }
      else lg++;
    }
    if (!regexec(&patternkey, line, 0, NULL, 0)) getflag = 1;
  }
  printf("INM has been read twice.\n");
  fflush(NULL);
  /*for (i = 0; i < refnLoci; i++) printf("line 213 of program: refnames[%d] = %s reflgs[%d] = %d\n", i, refnames[i], i, reflgs[i]);*/
  DAT = fopen(datafile, "r");
  printf("The main datafile %s has been opened\n", datafile);
  fflush(NULL);
  fgets(line, maxdatlinelength, DAT); /*There should be three header lines.*/
  fgets(line, maxdatlinelength, DAT);
  fgets(line, maxdatlinelength, DAT);
  if (strlen(line) > 1) {
    printf("The third line of %s is not blank but instead is:\n%s", datafile, line);
    return(1);
  }
  fgets(line, maxdatlinelength, DAT);
  fgets(line, maxdatlinelength, DAT);
  if (regexec(&patternast, line, 0, NULL, 0)) plantsTotal = strlen(line) - 1; //Account for the clinging newline.
  else {
    printf("We are seriously out of register with our data lines.  This line is %s", line);
    return(1);
  }
  maxdatnamelength = 0;
  while (fgets(line, maxdatlinelength, DAT)) {
    if (regexec(&patternast, line, 0, NULL, 0)) {
      if (strlen(line) - 1 != plantsTotal) {
        printf("Data line of deviant length %lu is %s", strlen(line) - 1, line);
        return(1);
      }
    }
    else {
      if (strlen(line) > maxdatnamelength) maxdatnamelength = strlen(line);
    }
  }
  printf("The value of plantsTotal in %s is %d and maxdatnamelength = %d\n", datafile, plantsTotal, maxdatnamelength);
  fflush(NULL);
  refchdata = (char **) calloc(refnLoci + 5, sizeof(char *));
  for (i = 0; i < refnLoci + 5; i++) refchdata[i] = (char *) calloc(plantsTotal + 5, sizeof(char));
  fseek(DAT, 0, SEEK_SET);
  fgets(line, maxdatlinelength, DAT); /*There should be three header lines.*/
  fgets(line, maxdatlinelength, DAT);
  fgets(line, maxdatlinelength, DAT);
  printf("Data headers have been read.\n");
  fflush(NULL);
  while (fgets(line, maxdatlinelength, DAT)) {
    //Single deflines are assumed to alternate with single data lines.
    //*5
    //BABBABHBAAABABABABHHAAAHAAAAABHAHHAHHBHBBAAABAHBAAABAA-HAAAHHABBHA-BHB-BHABBAHHHABBABA-ABBAHBBAAHHAABBHBAAAABBAAHAHHAHBABAHBHA-HHHBHAAAHABBAABBAHB--HBAHAAAAHBBABBABABHAAA-AAABBBBBH-AABBBBBBAAHA-BBAABBBHABABAABBABHBBAAAABBHHABHAHHAAHBHHAAHABBBABBBBAAAABABBABBBBHBHBHBHBHHBBBBHHBBHHABBAABABABBBHHBBAHABBABABABBBBBBABBABBBHBBAB-BABHAHBHHBBBBBHAABB-HABAAABBBBBHBABBAAHAAHHBBHBABHHAABHHABABHHBBHABH-BHBABABBHAAABABBBABBBABBBHHBBBAHAABAHAAABABHHBBHAHABHBBBB-ABB-BAAAAAAHBBHAHB-HBHABBHBAABAHBH-B-HBBBABHAHAH
    if (!regexec(&patternast, line, 0, NULL, 0)) {
      token = strtok(line, " |\t\n,");
      strcpy(datname, token + 1);
      if ((pos = strchr(datname, '\n')) != NULL) *pos = '\0';
      //printf("datafile leading token = %s datname = %s that is all\n", token, datname);
      for (k = 0; k < refnLoci; k++) {
        if (strncmp(datname, refnames[k], maxdatnamelength) == 0) { //fill in refchdata[k]
          //printf("Matching name found, datname = %s refnames[%d] = %s linkage group = %d\n", datname, k, refnames[k], reflgs[k]);
          fgets(line, maxdatlinelength, DAT);
          if ((pos = strchr(line, '\n')) != NULL) *pos = '\0'; //Take off the newline.
          for (i = 0; i < (int)strlen(line); i++) refchdata[k][i] = toupper(line[i]);
          refchdata[k][(int)strlen(line)] = '\0';
          //printf("Data for marker %s = refnames[%d] are %s\n", refnames[k], k, refchdata[k]);
          //printf("Data for marker %s = refnames[%d] are ", refnames[k], k);
          //for (i = 0; i < (int)strlen(line); i++) printf("%c", refchdata[k][i]);
          //printf("\n");
        }
      }
    }
    //else do nothing
  }
  //for (i = 0; i < refnLoci; i++) printf("line 273 of program: refnames[%d] = %s reflgs[%d] = %d refchdata[%d] = %s\n", i, refnames[i], i, reflgs[i], i, refchdata[i]);
  OUD = (FILE **) malloc(nlgs * sizeof(FILE *));
  for (i = 0; i < nlgs; i++) {
    OUD[i] = (FILE *) malloc(sizeof(FILE *));
    sprintf(outfilename, "%s%.3d%s", outfilestem, i, ".txt");
    /*printf("outfile name is %s\n", outfilename);*/
    OUD[i] = fopen(outfilename, "w");
  }
  fseek(DAT, 0, SEEK_SET);
  fgets(header1, maxdatlinelength, DAT); //This will be used intact.
  fgets(header2, maxdatlinelength, DAT); //This will be split to insert correct marker count once it is known
  fgets(header3, maxdatlinelength, DAT); //This should be a blank line, '^\n'.
  for (i = 0; i < nlgs; i++) fprintf(OUD[i], "%s%s%s", header1, header2, header3); //Correct nLoci later.
  fflush(NULL);
  k = 0;
  while (fgets(line, maxdatlinelength, DAT)) {
    if (!regexec(&patternast, line, 0, NULL, 0)) {
      token = strtok(line, " |\t\n,");
      strcpy(datname, token + 1);
    }
    else {
      //data are in line; compare to refchdata
      /*if (strncmp(datname, "g2501", 5) == 0) printf("data line for g2501 is %s", line); //line has a newline?*/
      mindiff = plantsTotal; //This would happen if every character differed between datum and reference.
      savedgroup = -1;
      if ((pos = strchr(line, '\n')) != NULL) *pos = '\0'; //Take off the newline.
      for (i = 0; i < strlen(line); i++) line[i] = toupper(line[i]);
      for (i = 0; i < refnLoci; i++) {
        currdiff = 0;
        for (j = 0; j < plantsTotal; j++) {
          if (line[j] != refchdata[i][j]) currdiff++;
        }
        if (currdiff < mindiff) {
          mindiff = currdiff;
          savedgroup = reflgs[i];
          strcpy(savedname, refnames[i]);
        }
        /*if (strncmp(datname, "g2501", 5) == 0) {
          printf("currdiff = %d mindiff = %d refnames[%d] = %s reflgs[%d] = %d refchdata[%d] = %s\n", currdiff, mindiff, i, refnames[i], i, reflgs[i], i, refchdata[i]);
        }*/
      }
      if (savedgroup < 0) printf("marker %s was not matched to any group\n", datname);
      fprintf(OUD[savedgroup], "*%s\n%s\n", datname, line);
      k++;
    }
  }
}

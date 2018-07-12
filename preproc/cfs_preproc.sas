/**/
/*
/* Preprocessing of CFS data
/* This file basically loads data from different xls-sources
/* and compiles a large database that is supposed to be the 
/* masterfile for any subsequent analyses 
/*
/* PH, 09/2016
/*
/**/

/**/
/* import section
/**/
/* load scr data from previous analyses */
proc import datafile="/folders/myshortcuts/mssm/tmp/cfs/preproc/cfs_scrfiles.xlsx" 
out=scr
dbms=xlsx replace;
*sheet="Data";
run;

/* load eprime data */
proc import datafile="/folders/myshortcuts/mssm/tmp/cfs/preproc/cfs_eprimedata_new.xlsx" 
out=eprime
dbms=xlsx replace;
*sheet="Data";
run;

/* load additional subject information */
proc import datafile="/folders/myshortcuts/mssm/tmp/cfs/preproc/cfs_data.xlsx" 
out=cfs
dbms=xlsx replace;
sheet='DataFinal';
*range="Data$A1:AQ110";
run;

/* load dcm estimates (using eSNA as outcome measure) */
proc import datafile="/folders/myshortcuts/mssm/tmp/cfs/preproc/cfs_dcmout3.csv" 
dbms=csv
  out=dcm3 replace;
  getnames=no;
run;

/* load dcm estimates (using aSNA as outcome measure) */
proc import datafile="/folders/myshortcuts/mssm/tmp/cfs/preproc/cfs_dcmout4.csv" 
dbms=csv
  out=dcm4 replace;
  getnames=no;
run;

/* load additional random sample of SCR peak2peak data scored by Jingchu */
proc import datafile="/folders/myshortcuts/mssm/tmp/cfs/preproc/cfs_randompick.xlsx" 
dbms=xlsx
  out=rand
  replace;
run;
/**/
/* end import
/**/


/**/
/* data preprocessing
/**/

/* rename dcm data */
data dcm;
set dcm3(rename=(Var1=id Var2=trial Var3=dcm1 Var4=dcm2 Var5=dcm3));
/* set dcm4(rename=(Var1=id Var2=trial Var3=dcm1 Var4=dcm2 Var5=dcm3 Var6=dcm4)); */
run;

/* sort eprime, scr, and cfs for merging */
proc sort data=eprime; by id;
proc sort data=scr; by id;
proc sort data=cfs; by id;

/* drop practice trials in eprime */
data eprime;
set eprime; 
if Block in (1 2 3 4 5 38) then delete;
run;


/* add a trial and stimulus variable in eprime */
proc sort data=eprime; by id block;
data eprime;
set eprime;
length stimulus $4;
by id;
if first.id then 
do;
	trial = 1;
end;
else 
do;
	trial + 1;
end;
if trial < 17 and stim eq 'CSminus' then stimulus = 'CS--';
if trial < 17 and stim eq 'CSplus' then stimulus = 'CS++';
if trial > 16 and stim eq 'CSminus' then stimulus = 'CS+-';
if trial > 16 and stim eq 'CSplus' then stimulus = 'CS-+';
run;


/* clean up cfs */
data cfs;
set cfs;
if id eq . then delete;
run;

/* merge cfs and eprime */
data final;
merge cfs eprime;
by id;
run;

/* merge final dataset with dcm by id and trial */
proc sort data=dcm; by id trial;
proc sort data=final; by id trial;
data final;
merge final dcm;
by id trial;
run;


/* add some variables to final dataset */
data final;
set final;

/* assign remaining ids */
if id in (949 7746 7838 7847 9362 9452) then group = 1;
if id in (3017 3982 5049 5508 6793 7121 7382 7908 8147 9076 9247 9525) then group = 0; 
if id = 8147999 then delete;

/*
visrespn = input(visresp, 8.);
confrespn = input(confresp, 8.);
visrtn = input(visrt, 8.);
confrtn = input(confrt, 8.);
length stage $16;
if trial < 17 then stage='Acquisition';
if trial > 16 then stage='Reversal';
*/

/* add information about former exclusion and regrouping status */
prevexcl=0;
prevgroup=group;
regrouped=0;
bleedlevel=0;
if id in (9697 891 2846 4554 6054 1493 1875 1758 9917 3942 1952 5115 1491 999) then prevexcl=1;
if id in (3653 8977 8961 3982 8147 4690 4398 9090 3553 5586 9256 3184 3937 2012) then do; prevgroup=1; regrouped=1; end;
if id in (204 524 654 2670 3324 3464 3695 3869 5049 6793 6946 7003 7245 7328 7908 8072 8334 9178
2387 6104 4116 949 9452 9659 7747 7847 2147 7714 9525 999 8688 7271 5765 6534 3763 8021) then prevexcl=1;

/* introduce bleedlevel to mark clear and probable bleed-throughs */
if group = 0 and visresp = 2 and confresp = 3 then bleedlevel=2;
if group = 0 and visresp = 2 and confresp = 2 then bleedlevel=1;
if group = 1 then bleedlevel=0;

length spider $1;
if stimulus eq 'CS++' then spider = 'A';
if stimulus eq 'CS--' then spider = 'B';
if stimulus eq 'CS+-' then spider = 'A';
if stimulus eq 'CS-+' then spider = 'B';
length groupc $12;
if group = 0 then groupc = 'Unaware';
if group = 1 then groupc = 'Aware';

/* complete missing information (as extracted from eprime files) */
if id in (204 891 949 2845 2846 3553 4554 5508 5765 6793 7791 8746 9525) then order = 'A';
if id in (2012 3019 3982 3586 4996 5049 5445 7063 7121 7838 7847 8961 9452) then order = 'B';
if order eq 'A' and stim eq '' and trial in (1 3 5 8 9 11 13 16 18 20 21 24 26 28 29 32) then spider='B';
if order eq 'A' and stim eq '' and trial in (2 4 6 7 10 12 14 15 17 19 22 23 25 27 30 31) then spider='A';
if order eq 'B' and stim eq '' and trial in (1 3 5 6 8 10 12 15 18 20 22 23 26 27 29 31) then spider='A';
if order eq 'B' and stim eq '' and trial in (2 4 7 9 11 13 14 16 17 19 21 24 25 28 30 32) then spider='B';

run;


data finalprobs;
set final(keep=id bleedlevel);
run;


/* sort by bleedlevel (descending) so that each subject is 
assigned its maximal problem level */
proc sort data=finalprobs; by id descending bleedlevel;
data finalprobs;
set finalprobs;
by id;
if not first.id then delete;
bleedstatus=bleedlevel;
run;

proc sort data=final; by id trial;
proc sort data=finalprobs; by id;
data final;
merge final finalprobs(drop=bleedlevel);
by id;
run;

/* add a variable "ctrial" to allow for trial-wise comparisons */
proc sort data=final; by id spider;
data final;
set final;
by id spider;
if first.spider then ctrial = 1;
else ctrial + 1;
run;

/* create a "cross-sectional" baseline dataset */
data finalunique;
set final;
where trial in (1);
run;

/* rename final to cfs and prepare for export */
data cfs;
set final;
run;

proc sort data=cfs; by id trial;
proc sort data=dcm; by id trial;
proc sort data=rand; by id;


data rand;
set rand;
by id;
if first.id then trial=1;
else trial+1;
run;

data cfs;
merge cfs dcm;
by id trial;
run;


/* drop unnecessary variables */
data cfs;
set cfs(drop=order included_in_final_data 
full_spider_breakthrough partial_spider_breakthrough
marker no_prior_reconsolidation_conditi stim block);
run;

/* rename some variables */
data cfs;
set cfs(rename=(stimulus=stim pref_order=order));
run;

/* delete all unwanted tables */
proc datasets library=work;
delete dcm dcm3 dcm4 eprime final finalprobs finalunique rand scr;
run;

/**
/* Plotting
/**
/* plot data to make sure everything is in place */

title 'Original group assignment';
proc freq data=cfs;
tables groupc;
where trial = 1 and not missing(dcm1) and spider not eq "";
run;

proc sort data=cfs; by id trial;
proc sgpanel data=cfs;
panelby groupc / novarname;
vline ctrial / response=dcm2 stat=mean group=spider limitstat=stderr;
where ctrial < 17;
run;

proc sgpanel data=cfs;
panelby groupc;
vbar stage / response=dcm2 stat=mean limitstat=stderr group=spider groupdisplay=cluster;
run;
title;




title 'Clear bleed-through subjects excluded';
proc freq data=cfs;
tables groupc;
where trial in (1) and not missing(dcm1) and bleedstatus not in (2) and spider not eq "" ;
run;

proc sort data=cfs; by id trial;
proc sgpanel data=cfs;
panelby groupc / novarname;
vbar stage / response=dcm2 stat=mean group=spider limitstat=stderr groupdisplay=cluster;
where ctrial < 17 and bleedstatus not in (2);
run;
title;




title 'Clear and possible bleed-through subjects excluded';
proc freq data=cfs;
tables groupc;
where trial in (1) and not missing(dcm1) and bleedstatus not in (1 2) and spider not eq "";
run;

proc sort data=cfs; by id trial;
proc sgpanel data=cfs;
panelby groupc / novarname;
vbar stage / response=dcm2 stat=mean group=spider limitstat=stderr groupdisplay=cluster;
where ctrial < 17 and bleedstatus not in (1 2);
run;
title;




title 'Subjects re-assigned as in previous submission';
proc freq data=cfs;
tables prevgroup;
where trial in (1) and not missing(dcm1) and not prevexcl=1 and spider not eq "";
run;

proc sort data=cfs; by id trial;
proc sgpanel data=cfs;
panelby prevgroup / novarname;
vbar stage / response=dcm2 stat=mean group=spider limitstat=stderr groupdisplay=cluster;
where ctrial < 17 and not prevexcl=1;
run;
title;




title 'Subjects pooled';
proc freq data=cfs;
tables confresp;
where not missing(dcm1) and spider not eq "";
run;

proc sort data=cfs; by id trial;
proc sgpanel data=cfs;
panelby stage / novarname;
vbar confresp / response=dcm2 stat=mean group=spider limitstat=stderr groupdisplay=cluster;
where not missing(confresp) and ctrial < 17;
colaxis values=("1" "2" "3") valuesdisplay=("Guess" "Unsure" "Confident");
run;
title;



/**/
/* export data to final database (excel sheet) */
/**/
proc export data=cfs
dbms=xlsx
  outfile="/folders/myshortcuts/mssm/tmp/cfs/cfs.xlsx" 
  replace;
run;

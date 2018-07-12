/**
/* CFS Analysis
/* 
/* PH, 09/2016
/*
/**

/* import masterfile */
proc import datafile="/folders/myshortcuts/projects/cfs/preproc/cfs.xlsx" 
dbms=xlsx
  out=cfs
  replace;
  sheet="cfs";
run;

/**
/* Plotting
/**
/* plot data to make sure everything is in place */

data cfs;
set cfs;
logdcm2 = log(dcm2+1);
run;

title 'Original group assignment';
proc freq data=cfs;
tables groupc;
where trial = 1 and not missing(dcm1) and spider not eq "";
run;

proc sort data=cfs; by id trial;
proc sgpanel data=cfs;
panelby groupc / novarname;
vline ctrial / response=logdcm2 stat=mean group=spider limitstat=stderr;
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
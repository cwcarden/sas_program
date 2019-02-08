*------------- bloom graphs.sas;
options noxwait obs=max ;

%let base = \\USPLA1FS001\Groups\PRPPLAIN\Bloom Graphs ;
filename macloc 'C:\SAS\Macros_Wall' ;
options sasautos=(sasautos,macloc) validvarname=V7 ;

%MACRO BLOOMIT_PREP(HYBRID=,MALE=,FEMALE=);

   %let master_file =_inbred_bloom.xls ;

   %let pdate = '14-Apr-2002'd ;   *--- THESE ARE CONSTANTS;
   %let fdate = '01-Jul-2002'd ;
   %let curryr = 2002 ;

   %let hybrid = %upcase(&hybrid) ;
   %let male = %upcase(&male) ;
   %let female = %upcase(&female) ;

   %importxls(file=&base\Data Files\&master_file,
              data=master,
              dbtype=excel2000,
              getnames=yes)

   proc sort data=master nodupkey ;
      by yr inbred pd bd ;
   run ;

   data female ;
      set master ;

      if (inbred="&female") ;

      Sex = 'F' ;
   run ;

   data male ;
      set master ;

      if (inbred="&male") ;

      Sex = 'M' ;
   run ;

   %sort(data=female,by=yr pd bd)
   %sort(data=male,by=yr pd bd)

   data unchecked ;
      set female male ;
      by yr pd bd;
   run ;

   %addstat(data=unchecked,var=pd,stat=N,newvar=parent_cnt,grpby=%str(yr,pd))

   %sort(data=unchecked,by=yr pd parent_cnt sex)

   proc transpose data=unchecked out=tran ;
      by yr pd parent_cnt ;
      var bd ;
      id sex ;
   run ;

   data balanced_pairs(keep=yr_pd) ;
      set tran ;

      if (M > . AND F > .) ;

      length yr_pd $5 ;
      yr_pd = compress(yr || '_' || pd) ;
   run ;

   %quotedta(data=balanced_pairs,var=yr_pd,listname=paired_yps,
             upcase=y,comma=y,quoted=y,countvar=listN)

   data matched ;
      set unchecked ;

      length yr_pd $5 ;
      yr_pd = compress(yr || '_' || pd) ;

      if (parent_cnt > 1 AND yr_pd IN (&paired_yps)) ;
   run ;

   %sort(data=matched,by=yr pd inbred)
   data trouble ;
      set matched ;
      by yr pd ;

      if (first.pd) then cnt = 0 ;

      cnt + 1 ;

      if (last.pd AND cnt > 2) then flag4delete = 'Y' ;
   run ;

   proc print data=trouble ;
      by yr ;
      id yr ;
      var inbred Sex pd bd flag4delete ;
      title "Trouble for &hybrid?" ;
   run ;

   data paired(keep=Hybrid Sex pd bd
               rename=(pd=plant bd=bloom)) ;
      set trouble ;

      length Hybrid $20 ;
      Hybrid = "&Hybrid" ;

      if (flag4delete='Y') then delete ;
   run ;

   %numobs(data=paired,count=obs4graph)

   %if (&obs4graph NE 0) %then %do ;

      data current ;
         length file $50 ;
         file = "&hybrid (&female).gif" ;
      run ;
  
      proc print data=paired ;
         var Hybrid Sex plant bloom ;
         title "Paired Data for &Hybrid" ;
      run ;

      proc sort data=paired ;
         by Hybrid Sex ;
      run ;

      ods output parameterestimates=parmests ;
      proc glm data=paired ;
         by Hybrid Sex ;
         model bloom=plant ;
         output out=sf1 predicted=yhat ;
         title "GLM Results for &Hybrid" ;
      run ;
      quit ;

      proc transpose data=parmests out=tran_parms ;
         by sex ;
         var estimate ;
         id parameter ;
      run ;

      filename fregeq temp ;
      filename mregeq temp ;
      data _null_ ;
         set tran_parms ;

         if (sex='F') then do ;
            file fregeq ;
            put 'yhat = ' intercept '+ (' plant +(-1) '*plant) ;' ;
         end ;
         else
         if (sex='M') then do ;
            file mregeq ;
            put 'yhat = ' intercept '+ (' plant +(-1) '*plant) ;' ;
         end ;
      run ;

      data sf2 ;
         set sf1 end=lastrec ;

         pdate = datejul(plant + juldate(&pdate)) ;
         fdate = datejul(yhat + juldate(&fdate)) ;

         label pdate='PLANTING DATE' fdate='FLOWERING DATE' ;

         output ;

         *------------- Adding beginning and ending points for the regression lines;
         if (lastrec) then do ;
            sex = 'F' ;
            hybrid = "&hybrid" ;
            plant = 0 ;
            %include fregeq ;
            pdate = datejul(plant + juldate(&pdate)) ;
            fdate = datejul(yhat + juldate(&fdate)) ;
            output ;

            sex = 'F' ;
            hybrid = "&hybrid" ;
            plant = 80 ;
            %include fregeq ;
            pdate = datejul(plant + juldate(&pdate)) ;
            fdate = datejul(yhat + juldate(&fdate)) ;
            output ;

            sex = 'M' ;
            hybrid = "&hybrid" ;
            plant = 0 ;
            %include mregeq ;
            pdate = datejul(plant + juldate(&pdate)) ;
            fdate = datejul(yhat + juldate(&fdate)) ;
            output ;

            sex = 'M' ;
            hybrid = "&hybrid" ;
            plant = 80 ;
            %include mregeq ;
            pdate = datejul(plant + juldate(&pdate)) ;
            fdate = datejul(yhat + juldate(&fdate)) ;
            output ;
         end ;
      run ;

      data Hybrid90_&hybrid ;
         set sf2 end=lastrec ;
         by Hybrid ;

         ptscnt + 1 ;

         if (last.Hybrid) then do ;
            npts = (ptscnt / 2) - 2 ;
            ptscnt = 0 ;
         end ;

         format fdate pdate date5. ;
      run ;

      proc sort data=Hybrid90_&hybrid ;
         by Hybrid ;
      run ;


      %sort(data=hybrid90_&hybrid,by=sex pdate)

      proc datasets lib=work nolist ;
         delete master female male unchecked tran
                balanced_pairs matched trouble
                paired current sf1 sf2 ;
      run;
      quit ;

   %end ;
   %else %do ;
      data _null_ ;
             file print ;
                 title "No Pairs Report" ;
             put "There were no paired observations to graph for:" //
                     "Hybrid: &hybrid   Female: &female   Male: &male" ;
          run ;
      %end ;


%MEND BLOOMIT_PREP;


%MACRO BLOOMIT_POST(HYBRID=,MALE=,FEMALE=);

   %let master_file =_inbred_bloom.xls ;

   %let pdate = '14-Apr-2002'd ;   *--- THESE ARE CONSTANTS;
   %let fdate = '01-Jul-2002'd ;
   %let curryr = 2002 ;

   %let hybrid = %upcase(&hybrid) ;
   %let male = %upcase(&male) ;
   %let female = %upcase(&female) ;

   %importxls(file=&base\Data Files\&master_file,
              data=master,
              dbtype=excel2000,
              getnames=yes)

   proc sort data=master nodupkey ;
      by yr inbred pd bd ;
   run ;

   data female ;
      set master ;

      if (inbred="&female") ;

      Sex = 'F' ;
   run ;

   data male ;
      set master ;

      if (inbred="&male") ;

      Sex = 'M' ;
   run ;

   %sort(data=female,by=yr pd bd)
   %sort(data=male,by=yr pd bd)

   data unchecked ;
      set female male ;
      by yr pd bd;
   run ;

   %addstat(data=unchecked,var=pd,stat=N,newvar=parent_cnt,grpby=%str(yr,pd))

   %sort(data=unchecked,by=yr pd parent_cnt sex)

   proc transpose data=unchecked out=tran ;
      by yr pd parent_cnt ;
      var bd ;
      id sex ;
   run ;

   data balanced_pairs(keep=yr_pd) ;
      set tran ;

      if (M > . AND F > .) ;

      length yr_pd $5 ;
      yr_pd = compress(yr || '_' || pd) ;
   run ;

   %quotedta(data=balanced_pairs,var=yr_pd,listname=paired_yps,
             upcase=y,comma=y,quoted=y,countvar=listN)

   data matched ;
      set unchecked ;

      length yr_pd $5 ;
      yr_pd = compress(yr || '_' || pd) ;

      if (parent_cnt > 1 AND yr_pd IN (&paired_yps)) ;
   run ;

   %sort(data=matched,by=yr pd inbred)
   data trouble ;
      set matched ;
      by yr pd ;

      if (first.pd) then cnt = 0 ;

      cnt + 1 ;

      if (last.pd AND cnt > 2) then flag4delete = 'Y' ;
   run ;

   data paired(keep=Hybrid Sex pd bd
               rename=(pd=plant bd=bloom)) ;
      set trouble ;

      length Hybrid $20 ;
      Hybrid = "&Hybrid" ;

      if (flag4delete='Y') then delete ;
   run ;

   %numobs(data=paired,count=obs4graph)

   %if (&obs4graph NE 0) %then %do ;

      %G2F

          data _null_ ;
         set Hybrid90_&hybrid ;

         if (npts > 0) then call symput('dp_per_line',left(put(npts,3.0))) ;
      run ;

      proc template ;
         define statgraph bloom ;
            begingraph ;

               entrytitle "&female/&male" / textattrs=(color=blue family="Calibri" size=16) ;
               entrytitle "Hybrid=&hybrid" / textattrs=(color=blue family="Calibri" size=12) ;

               discreteattrmap name="lines" / ignorecase=true ;
                  value "F" / lineattrs=(color=red pattern=1) ;
                  value "M" / lineattrs=(color=blue pattern=2) ;
               enddiscreteattrmap ;

               discreteattrvar attrvar=grouplines var=sex attrmap="lines" ;

               layout overlay / yaxisopts=(label="Flowering Date"
                                           labelattrs=(FAMILY="Calibri")
                                           griddisplay=on
                                           timeopts=(viewmin="30JUN&curryr"d viewmax="15SEP&curryr"d
                                                     interval=week tickvalueformat=date5.
                                                     minorticks=true minortickinterval=day))
                                xaxisopts=(label="Planting Date"
                                           labelattrs=(FAMILY="Calibri")
                                           griddisplay=on
                                           timeopts=(viewmin="17APR&curryr"d viewmax="03JUL&curryr"d
                                                     interval=week tickvalueformat=date5.
                                                     minorticks=true minortickinterval=day)) ;
                  regressionplot y=fdate x=pdate / group=grouplines name="bloom";

                  discretelegend "bloom" / title="Sex" ;
               endlayout ;

               entryfootnote "Number of data points per line : &dp_per_line" ;

            endgraph ;
         end ;
      run ;

      proc sgrender data=Hybrid90_&hybrid template=bloom ;
      run ;

   %end ;
   %else %do ;
      data _null_ ;
             file print ;
                 title "No Pairs Report" ;
             put "There were no paired observations to graph for:" //
                     "Hybrid: &hybrid   Female: &female   Male: &male" ;
      run ;
      title ;
   %end ;

%MEND BLOOMIT_POST;


%importxls(file=&base\Data Files\hybrids_to_graph.xls,
                data=bloomhybs,
                dbtype=excel2000,
                getnames=yes)

data bloomhybs ;
   set bloomhybs ;

   if (hybrid='') then delete ;
run;

options obs=max noxwait ;
filename tmpcall1 temp ;
data _null_ ;
   set bloomhybs ;

   *------------- FEMALES WITH AN 'F' TOWARD THE END OF THEIR NAME
                  MUST HAVE THOSE 'F's CHANGED TO 'G's, BUT FEMALES
                  THAT *START* WITH 'F' MUST BE HANDLED CAREFULLY;

   fpos = indexc(female,'F') ;
   put female= fpos= ;

   if (fpos > 1) then female = tranwrd(female,'F','G') ;
   else
   if (fpos=1) then female = compress('F' || tranwrd(substr(female,2),'F','G')) ;

   file tmpcall1 ;
   *file print ;


   put '%' 'bloomit_prep(' hybrid= +(-1) ', ' male= +(-1) ', ' female= +(-1) ')' ;
run ;
options obs=max ;

%include tmpcall1 ;

options obs=max noxwait ;
filename tmpcall2 temp ;
data _null_ ;
   set bloomhybs ;

   *------------- FEMALES WITH AN 'F' TOWARD THE END OF THEIR NAME
                  MUST HAVE THOSE 'F's CHANGED TO 'G's, BUT FEMALES
                  THAT *START* WITH 'F' MUST BE HANDLED CAREFULLY;

   fpos = indexc(female,'F') ;
   put female= fpos= ;

   if (fpos > 1) then female = tranwrd(female,'F','G') ;
   else
   if (fpos=1) then female = compress('F' || tranwrd(substr(female,2),'F','G')) ;

   file tmpcall2 ;
   *file print ;


   put '%' 'bloomit_post(' hybrid= +(-1) ', ' male= +(-1) ', ' female= +(-1) ')' ;
run ;
options obs=max ;

title;

ods graphics / width=10in reset=index ;
options orientation=landscape ;

ods powerpoint file="&base\Graphs\Summary GTL_2016.pptx" ;

   %include tmpcall2 ;

ods powerpoint close ;

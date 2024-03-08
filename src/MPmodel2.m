classdef MPmodel2 < handle
   %MPmodel2 class
   %PUBLIC FUNCTIONS:
   %  CONSTRUCTOR: obj = MPmodel2(MP_XYY, dataname)
   %  showplate(obj)
   %...
   %DEPENDENCIES
   %GRmodelOD.m class
   %----------------------------------------------------------------------
   properties (Access = private)
      GR_M        %GR model class     
      MP_XYY      %Microplate data matrix 
      MP_XYN      %Microplate in 8x12xN data matrix format 
      timevec     %Nx1 time points 0..hrs
      TPs         %number of time points (including time zero)
      firstN      %points used for initial OD calculation
      dataname    %string with name of data object
   end% properties

   methods %public methods------------------------------------------------
      function obj = MPmodel2(MP_XYYtable, dataname)
         %take Nx97 Table and make MP_XYY matrix, MP_XYN 3D matrix, etc. 
         obj.MP_XYY = MP_XYYtable{:,:};      %save original data (8 x 97)
         obj.timevec = MP_XYYtable{:,1};     %extract time vector
         obj.TPs = length(obj.timevec);      %number of data points
         obj.MP_XYN = zeros(8,12,obj.TPs);   %create XYN 3D matrix
         obj.setXYN(0);                      %Init XYN with zero blank 
         obj.firstN = 5;                     %pts used for minOD calc
            %NOTE: firstN = 0 means just get minimum value of all data pts
         obj.GR_M = GRmodelOD();             %new GRmodel class object
         obj.dataname = dataname;            %set data name
      end %function

      function dataname = getDataName(app)
         dataname = app.dataname;
      end %function

      function XYY = getXYY(obj)
         XYY = obj.MP_XYY;                   %return original XYY matrix
      end

      function tilegraph = showplate(obj)
         close all
         t = tiledlayout(8,12,'TileSpacing','none'); %8x12 grid, no spacing   
         xlabel(t,"Time (h)");               %x axis label
         ylabel(t,"Optical Density (600 nm)"); %y axis label
         Ydata = obj.MP_XYY(:,2:end);        %just get Y data values
         ymax_val = max(Ydata,[],"all");     %max Y val for Y axis scale
         for i=2:97                          %for each data col
            nexttile                         %add new tile to graph
            plot(obj.MP_XYY(:,1),obj.MP_XYY(:,i),'.k'); %plot in tile
            ylim([0,ymax_val]);              %set Y limit
            set(gca,'xtick',[],'ytick',[]);  %add tick marks
         end %for
         rowvals = ['A','B','C','D','E','F','G','H'];
         %rowvals = ["A","B","C","D","E","F","G","H"];
         colvals = ["1","2","3","4","5","6","7","8","9","10","11","12"];
         t.Children = flipud(t.Children);
         for j=1:8
            rnum = tilenum(t,j,1);
            ylabel(t.Children(rnum),rowvals(j),'Rotation',0, ...
               'FontWeight','bold');
         end %for
         for k=1:12
            cnum = tilenum(t,1,k);
            title(t.Children(cnum),colvals(k),'FontWeight','bold');
         end %for
         title(t,obj.dataname,'FontSize',12,'FontWeight','bold'); %title
         tilegraph = t;
      end %function

      function closegraph(~)
         close all;
      end

      function blankdata = getBlankData(obj,rowvec,col)
         %returns [mean stdev] 
         nrows = length(rowvec);             %row vector
         tempmat = zeros(obj.TPs,nrows);     %matrix for zero data
         for i=1:nrows                       %each row
            for j=1:obj.TPs                  %each TP
               tempmat(j,i) = obj.MP_XYN(rowvec(i),col,j); %build matrix
            end %for  
         end %for
         blankdata = [mean(tempmat,"all") std(tempmat,0,"all")]; %mean, std
      end %end

      function resetFirstN(obj,firstN)
         %NOTE: firstN = 0 means minOD is minimum from all Y values
         % otherwies from the first N values. 
         obj.firstN = firstN;
      end %function

      function curveData = getCurveData(obj, rowval, colvec, blankval)
         obj.setXYN(blankval);               %subtract blank OD from data
         
         crvmat = obj.getSubXYY(rowval, colvec);  %curve data (XYY);

         % matsize = size(crvmat);
         % if matsize(1) < obj.firstN
         %    crvmat = zeros(obj.firstN,);
         % end

         %crvmat has Y cols with blank subtracted and no zero or negative y
         % values included
         obj.setXYN(0);                      %reset blank to zero
         tempY = crvmat(:,2:end);            %Y col only matrix
         %with blank = 0 should be original XYY:
         curveData.oriXYY = obj.getSubXYY(rowval,colvec);
         
         oriX = curveData.oriXYY(:,1);
         trimX = crvmat(:,1);
         
         missvec = setdiff(oriX,trimX);
         curveData.missingTPs = missvec;
         
         % obj.firstN
         
         %NOTE: This line fails if background is too high? 
         curveData.OKfirstN = isequal(oriX(1:obj.firstN), ...
            trimX(1:obj.firstN));
         
         curveData.firstNvals = crvmat(1:obj.firstN,:);
         curveData.XYY = crvmat;             %XYY data
         curveData.Xcol = crvmat(:,1);       %X values only
         curveData.Ycols = tempY;            %Y colmatrix
         [Npts,NYcols] = size(tempY);        %dimensions of Y col matrix
         if Npts == obj.TPs 
            curveData.trimbool = false;
         else
            curveData.trimbool = true;
         end %if
         curveData.rowval = rowval;          %record XYY row used
         curveData.colvec = colvec;          %record XYY col(s) used
         curveData.blankval = blankval;      %record blank used
         curveData.Ndp = Npts;               %number of data points
         curveData.NYcols = NYcols;          %number of Y columns
         curveData.LnYcols = log(tempY);     %convert to Ln vals
         LnXYY = [curveData.Xcol, curveData.LnYcols]; %LnXYY data
         curveData.LnXYY = LnXYY;            %record LnXYY
         curveData.LnXYall = obj.getCombinedXYY(LnXYY); %XY with all reps
         curveData.meanYcol = mean(tempY,2); %mean of each row 
         curveData.meanLnYcol = mean(curveData.LnYcols,2); %col mean of row
         %assume firstN is less than total number of rows
         tempLnY = curveData.LnYcols(1:obj.firstN,:); %first five Ln rows
         curveData.LnMinVec = min(tempLnY,[],1); %minimum of each col
         curveData.LnMaxVec = max(curveData.LnYcols,[],1); %max each col
         curveData.meanLnMin = mean(curveData.LnMinVec); %mean for min/col
         curveData.meanLnMax = mean(curveData.LnMaxVec); %mean max
         %assume firstN is less than total number of rows
         curveData.meanfirstNLnMinOD = ...
            min(curveData.meanLnYcol(1:obj.firstN));
         curveData.meanfirstNMinOD = ...
            min(curveData.meanYcol(1:obj.firstN));
      end % function
       
      function plotCurves(obj,rowval,colvec,blankval,lngraphbool)

         %plot curves for debugging, Lngraphbool = 1 plot Ln values,
         % Lngraphbool = 0, plot untransformed OD
         clf("reset");
         crv = obj.getCurveData(rowval,colvec,blankval); %get curve(s)
         hold on                             %plot multiple XY as needed
         if lngraphbool
            for i=1:crv.NYcols                  %for each curve
               plot(crv.Xcol,crv.LnYcols(:,i),'-ok'); %plot data
            end %for
            ylabel("Ln Optical Density (600 nm)"); %Y axis label
         else
            for i=1:crv.NYcols                  %for each curve
               plot(crv.Xcol,crv.Ycols(:,i),'-ok'); %plot data
            end %for
            ylabel("Optical Density (600 nm)"); %Y axis label
         end
         titlestr = strcat(obj.dataname," [",string(rowval),", ", ...
            string(colvec), "]");            %make title
         title(titlestr);                    %set title
         xlabel("Time (h)");                 %X axis lable
         
         hold off                            %end hold
      end %function

      function repTable = getRepTable(obj,colvec,blankval,processbyN)
         postError = 0;
         datamat = zeros(8,12);               %matrix for data table
         Ncols = length(colvec);             %number of data sets
         tempColData = zeros(Ncols,4);       %data: lag,min,maxOD, LnBool
         for i=1:8                           %for each row
            for j=1:Ncols                    %for each col
               tempcrv = obj.getCurveData(i,colvec(j),blankval);
               crv = [tempcrv.Xcol, tempcrv.LnYcols];
               tempdata = obj.GR_M.getMaxGRdata(crv,processbyN, ...
                  obj.firstN); %get GR data structure 
               tempColData(j,1) = tempdata.LagTime; %save lag for each col
               tempColData(j,2) = tempdata.MinOD; %save minOD for each col
               tempColData(j,3) = tempdata.MaxOD; %save maxOD for each col
               tempColData(j,4) = tempcrv.OKfirstN;
            end %
            LagTime = mean(tempColData(:,1)); %mean of lag data from cols
            LagTimeStdErr = std(tempColData(:,1))/sqrt(Ncols); %stderr lag
            MinOD = mean(tempColData(:,2));  %mean of minOD from cols
            MinODStdErr = std(tempColData(:,2))/sqrt(Ncols); %stderr
            MaxOD = mean(tempColData(:,3));  %mean of maxOD from cols
            MaxODStdErr = std(tempColData(:,3))/sqrt(Ncols); %stderr
            NoBlankErrors = sum(tempColData(:,4));
            %get growth rate and doubling time from regression with stats
            processbyNr = processbyN * Ncols; %reps for processing by N
            rowcrv = obj.getCurveData(i,colvec,blankval); %all col data
            tempdata = obj.GR_M.getMaxGRdata(rowcrv.LnXYall,processbyNr, ...
               obj.firstN); %get GR data from all cols in one XY set
            datamat(i,1) = i;                %row number
            datamat(i,2) = tempdata.GR;      %get GRdata
            datamat(i,3) = tempdata.StdErrGR; %Gr std err
            datamat(i,4) = tempdata.DblTime; %get doubling time
            datamat(i,5) = tempdata.StdErrDblTime; %doubling time stderr
            datamat(i,6) = LagTime;          %assign col rep lag data
            datamat(i,7) = LagTimeStdErr;    %std err lag from col reps
            datamat(i,8) = MinOD;            %min OD col rep data
            datamat(i,9) = MinODStdErr;      %std err minOD from col reps
            datamat(i,10) = MaxOD;           %max OD col rep data
            datamat(i,11) = MaxODStdErr;     %std err maxOD from col reps
            if ~NoBlankErrors || MinOD < -6.2 %blank error, min OD too low
                datamat(i,6) = 0;            %lag undefined
                postError = 1;               %set errop
            end %if
            datamat(i,12) = postError;       %post error value
         end %for                            %Make table:
         repTable = array2table(datamat,"VariableNames",{'PlateRow', ...
            'GR','StdErrGR','DblTime','StdErrDblTime','LagTime', ...
            'StdErrLagTime','MinOD','StdErrMinOD','MaxOD' ...
            'StdErrMaxOD','BlankError'});    %assign table and heading
      end %function

      function ODcurveTable = getODcurveTable(obj,rowval,colvec, ...
            blankval, processbyN)
         crvdata = obj.getCurveData(rowval,colvec,blankval);
         %GRdata = obj.getLnGRdata(rowval,colvec,blankval,processbyN);
         Rtab = obj.getRepTable(colvec,blankval,processbyN);
         Xcol = crvdata.Xcol;
         meanYcol = crvdata.meanYcol;
         meanLnYcol = crvdata.meanLnYcol;
         GR = Rtab{rowval,"GR"};
         Lag = Rtab{rowval,"LagTime"};
         modelXY = obj.getGompertz(crvdata,GR,Lag);
         modelY = modelXY(:,2);
         Rep = crvdata.Ycols;
         Reps = array2table(Rep);
         MLGmat = [meanYcol,meanLnYcol,modelY];
         MLG = array2table(MLGmat,"VariableNames", ...
            {'MeanOD','LogOD','Gompertz'});
         Time = array2table(Xcol,"VariableNames",{'Time'});
         ODcurveTable = [Time,Reps,MLG];   
      end %function

      function modelXY = getGompertz(~,curveData,GR,Lag)
         %use curveData object, a GR value, and a Lag value to plot
         % Gompertz curve
         modelXY = zeros(curveData.Ndp,2);
         modelXY(:,1) = curveData.Xcol;
         start = min(curveData.LnMinVec);  
         A = curveData.meanLnMax + abs(start);
         e = exp(1);
         mu = GR;
         Lambda = Lag;
         for i=1:curveData.Ndp
            modelXY(i,2) = start + ...
               A*exp(-1*exp((mu*e/A)*(Lambda - modelXY(i,1))+1));
         end %for 
         if isnan(modelXY(:,2))
            modelXY(:,2) = 0;
         end %if

      end %function

   end %public methods----------------------------------------------------

   methods (Access = private)

      function setXYN(obj, blankval)
         %make XYN from XYY data, subtract blank value from Y data
         colindex = 2;                       %temp col index, first Y col
         for i=1:8                           %each row
            for j=1:12                       %each col 
               %subtract blank value form XYY OD data for each Y col
               % and load XYN matrix
               obj.MP_XYN(i,j,:) = obj.MP_XYY(:,colindex)-blankval; %load Y
               colindex = colindex + 1;      %advance to next col (2-97)
            end %for
         end %for 
      end %function

      function XYY = getSubXYY(obj,rowval,colvec)
         %extract a set of XYY data from XYN matrix
         ncols = length(colvec);             %number of cols to process
         XYY = zeros(obj.TPs,ncols+1);       %empty submatrix
         XYY(:,1) = obj.timevec;             %add times to first col
         colindex = 2;                       %first Y column index
         for i=1:ncols                       %for each Y col
            for j=1:obj.TPs                  %get tp data from XYN
               XYY(j,colindex) = obj.MP_XYN(rowval,colvec(i),j); 
            end %for
            colindex = colindex + 1;         %advance col index
         end %for
         %-----------------------------------------------------------------
         %test: remove any rows with negative Y values
         indx = any(XYY(:,2:end)<=0,2); %logical col vec, 1 means Ycol has <= 0
         XYY(indx,:) = [];  %remove all rows with Y col having val <= 0
         %check to see that there is data left
         matsize = size(XYY, 1);
         if matsize < obj.firstN %if data removed, something to process
            XYY = zeros(obj.firstN,ncols+1);
            XYY(:,1) = 1:obj.firstN;
            XYY(:,2:ncols+1) = 0.00001;
         end
         %-----------------------------------------------------------------
      end %function

      function XYall = getCombinedXYY(~,XYY)
         %return Nx2 with one X col and one Y col with all Y reps 
         [Nrows, Ncols] = size(XYY);         %get size of XYY
         NYcols = Ncols - 1;
         XYlen = Nrows*NYcols;               %XYlen = datapoints * Y cols
         resindex = 1;
         XYall = zeros(XYlen,2);             %result matrix with all zeros
         for i=1:NYcols                      %for each Y col
            for j=1:Nrows                    %each row in a Y col
               XYall(resindex,1) = XYY(j,1);
               XYall(resindex,2) = XYY(j,i+1); %get Y data for row,col
               resindex = resindex + 1;      %advance XYall row index
            end %for
         end %for
         XYall = sortrows(XYall,1);
      end %function

   end %private properties -----------------------------------------------

end %class def
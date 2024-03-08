classdef GRmodelOD < handle
   %Growth rate calculations class
   %Constructor requires XY data matrix (Nx2)
   %NOTE: Program converts OD data to Ln for processing
   %PUBLIC METHODS:
   %  GRmodelOD(), returns obj, constructor
   %  getMaxGRdata(processbyN, firstN), returns GRdata structure or 
   %   empty struct if 
   %GRdata structure:
   %  GRdata.XY = original XY data
   %  GRdata.Npts = number of data points
   %  GRdata.GR = max growth rate
   %  GRdata.StdErrGR = stderr GR
   %  GRdata.DblTime = doubling time
   %  GRdata.LagTime = lag time
   %  GRdata.MinOD = minimum OD from Y vals (min of firstN)
   %  GRdata.MaxOD = maximum OD from Y vals
   %  GRdata.DPused = corresponding OD data points for set with max slope 
   %  GRdata.startTime = time val for first data point of max slope set
   %  GRdata.startIndex = index for startTime
   %  GRdata.endTime = time val for end data point of max slope set
   %  GRdata.endIndex = index for endTime
   %  GRdata.firstNforMinOD = first N data points for calc of Min value
   %  GRdata.LnIntercept = intercept of regression line for GR (slope)
   %  GRdata.LnStdErrIntercept = stderr of intercept
   %  GRdata.LnRsq = Rsquared value for GR calc regression
   %  GRdata.LnResiduals = residuals of regression line fit for GR calc
   %----------------------------------------------------------------------

   properties %(Access = private)
      XY       %XY data matrix (Nx2) with time, OD values
      Nrow     %number of X vals 
   end %properties

   methods %public methods
       
      function obj = GRmodelOD()
         %constructor
         obj.XY = [];                        %empty matrix
         obj.Nrow = 0;                       %number of rows
      end %constructor

      function GRdataStruct = getMaxGRdata(obj, XY, processbyN, ...
            firstN)
         %Calculates sequential slopes with Ln values
         %NOTE: GRdata has results for Ln OD values not OD
         %ASSUMPTION: nsets is greater than 1
         XY(any(isnan(XY),2),:) = [];           %
         obj.setODdata(XY)                   %XY for processing
         startindex = 1;                     %save the first start index
         istart = 1;                         %initial start index
         iend = processbyN;                  %initial end index
         LRdata = obj.getLRegdata(istart,iend); %initial regression
         istart = istart + 1;                %advance start
         iend = iend + 1;                    %advance end
         while iend <= obj.Nrow          %check remaining sets  
            LRtemp = obj.getLRegdata(istart,iend); %calc GR data
            if LRtemp.Slope > LRdata.Slope   %check for largest slope
               startindex = istart;          %save start index
               LRdata = LRtemp;              %save regression data
            end
            istart = istart + 1;             %advance start
            iend = iend + 1;                 %advance end
         end %for
         %fill GRdata object with final results from best data set
         endval = startindex + processbyN - 1;
         GRdataStruct = obj.setGRdataStruct(LRdata, startindex, ...
            endval, firstN);
      end %function

   end %end public methods-----------------------------------------------

   methods (Access = private)

      function setODdata(obj, XY)
         obj.XY = XY;                        %assign raw data (Nx2)
         [obj.Nrow, ~] = size(XY);           %set number of rows
      end %function

      function GRdata = setGRdataStruct(obj,linRegData, ...
            indexstart,indexend,firstN)
         Ymin = obj.getMinOD(firstN);
         GRdata.XY = obj.XY;                 %original XY data
         GRdata.Npts = obj.Nrow;             %number of row values
         GRdata.GR = linRegData.Slope;       %max slope of regression
         GRdata.StdErrGR = linRegData.StdErrSlope; %stderr from LinReg
         GRdata.DblTime = 0;                 %default
         if GRdata.GR > 0                    %never divide by zero
             GRdata.DblTime = log(2)/GRdata.GR; %reset to Ln2/GR
             se_x2 = GRdata.StdErrGR^2;      %square of std err
             GRsq = GRdata.GR^2;             %square of GR
             Firstderiv = -1*log(2)/GRsq;
             %based on 'delta method' (first order approximation only)
             %see Gustaf Hendeby and Fredrik Gustafsson white paper
             GRdata.StdErrDblTime = sqrt(se_x2*(Firstderiv^2)); % + ...
         else
            GRdata.DblTime = 0;
            GRdata.StdErrDblTime = 0;
         end %if
         GRdata.LagTime = 0;                 %minimum lag allowed
         if linRegData.Slope > 0             %never divide by zero
            GRdata.LagTime = (Ymin - ...
               linRegData.Intercept)/linRegData.Slope; %reset
            if GRdata.LagTime < 0
               GRdata.LagTime = 0;
            end %if
         end %if
         GRdata.MinOD = Ymin;                %min OD
         GRdata.MaxOD = max(obj.XY(:,2));    %max OD
         GRdata.PtsUsed = [obj.XY(indexstart:indexend,1), ...
            obj.XY(indexstart:indexend,2)];  %OD pts used  
         GRdata.Times = [obj.XY(indexstart,1), ...
            obj.XY(indexend,1)];             %time pts used for GRmax
         GRdata.Indices = [indexstart, indexend]; %indices used for GRmax
         GRdata.LnIntercept = linRegData.Intercept; %assign intercept
         GRdata.LnStdErrIntercept = linRegData.StdErrIntercept;
         GRdata.LnRsq = linRegData.Rsquared; %set Rsq for regression
         GRdata.LnResiduals = linRegData.Residuals; %Residuals
      end %function

      function LRegdata = getLRegdata(obj,starttime,endtime)
         %do linear regression for sub set of data (start to end X vals)
         % using Ln OD values
         dataXvals = obj.XY(starttime:endtime,1);  %x vals
         dataYvals = obj.XY(starttime:endtime,2);  %use Ln vals
         datamat = [dataXvals dataYvals];    %make Nx2 matrix
         LRegdata = LinReg(datamat);         %get regression data
      end %function

      function minOD = getMinOD(obj,firstN)
         %get only min of firstN data values in matrix 
         if firstN == 0
            minOD = min(obj.XY(:,2),[],"all");
         else
            minOD = min(obj.XY(1:firstN,2:end),[],"all");
         end
      end %function

      function maxOD = getMaxOD(obj)
         maxOD = max(obj.XY(:,2),[],"a");
      end %function
   
   end %end private methods


end %class
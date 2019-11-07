function exportCTF(ebsd,fName,varargin)
%function exportCTF(ebsd,fName,varargin)
%Exporting EBSD data from mTex ebsd-object to Channel 5 Text File (ctf)
%The data in 'ctf' can for instance be opened with Channel 5 and Atex or
%further converted to 'ang' format for opening with Edax OIM
% *** Input arguments ***
%ebsd:      mTex ebsd-object
%fName:     Filename, optionally including relative or absolute path
%varargin:  Variable argument in
%           "'params',cprStruct" - Structure with properties from cpr-file import
%           "'params','manual'"  - String indicating prompt for manual
%                                  import of microscopy parameters
%           "'flipud',[0|1]"     - Flag 0 or 1 for flipping ebsd spatial
%                                  upside down (default: 0)
%           "'flipud',[0|1]"     - Flag 0 or 1 for flipping ebsd spatial
%                                  left right (default: 0)
% ************************************************************************
% Dr. Frank Niessen, University of Wollongong, Australia, 2019
% contactnospam@fniessen.com (remove the nospam to make this email address work)
% Acknowledgements go to Dr. Azdiar A. Gazder, University of Wollongong, Australia
% Version 1.0 - Published 18/04/2019
% ************************************************************************
scrPrnt('SegmentStart','Exporting ''ctf'' file');
%% Varargin
round0Thrsh = 1e-6;                                                        %Initialize threshold for rounding negative close to 0 x and y coordinates
params = [];                                                               %Initialize as empty
flip.ud = 0;                                                               %No flipping upside down
flip.lr = 0;                                                               %No flipping left right
for argidx = 2:2:(nargin + nargin(mfilename)+1) %Scan Varargin parameters
   switch varargin{argidx-1}
       case 'params'
            params = varargin{argidx};
       case 'flipud'
            flip.ud = varargin{argidx};
       case 'fliplr'
            flip.lr = varargin{argidx};
   end
end
%% Pre-processing
scrPrnt('Step','Collecting data');
if flip.ud %Flip spatial ebsd data upside down
    ebsd = flipud(ebsd);                                                   %Flipping ebsd data upside down
    scrPrnt('Step','Flipping EBSD spatial data upside down');
end
if flip.lr %Flip spatial ebsd data left right
    ebsd = fliplr(ebsd);                                                   %Flipping ebsd data left right
    scrPrnt('Step','Flipping EBSD spatial data left right');
end

ebsdGrid = ebsd.gridify;                                                   %Get gridified version of ebsd map
Laue = {'-1','2/m','mmm','4/m','4/mmm',...
        '-3','-3m','6/m','6/mmm','m3','m3m'};                              %11 Laue Groups
Laue = regexprep(Laue,'-','');                                             %Remove minuses in Laue Group names
% *** Get microscope acquisition parameters
AcquParam.Str = {'Mag','Coverage','Device','KV','TiltAngle','TiltAxis',...
                'DetectorOrientationE1','DetectorOrientationE2','DetectorOrientationE3',...
                'WorkingDistance','InsertionDistance'};                    %Cell-String array with Acquisition parameters
AcquParam.Fmt = {'%.4f','%.0f','%s','%.4f','%.4f','%.0f','%.4f','%.4f',...
                 '%.4f','%.4f','%.4f'};                                    %Cell-String array with Acquisition parameter formats                
%% Get Microscopy Acquisition Parameters
if isempty(params) %No parameters available
    scrPrnt('SubStep','Microscope acquisition parameters not available');
    AcquParam.Data(1:11) = {0};                                            %Filling in zeros
elseif isstruct(params) && isfield(params,'job') &&  isfield(params,'semfields')  %Cpr-file parameter structure from import
    scrPrnt('SubStep','Microscope acquisition parameters imported from Cpr-parameter structure');
    AcquParam.Data{1} = params.job.magnification;                          %Magnification
    AcquParam.Data{2} = params.job.coverage;                               %Coverage
    AcquParam.Data{3} = params.job.device;                                 %Device
    AcquParam.Data{4} = params.job.kv;                                     %Acceleration Voltage
    AcquParam.Data{5} = params.job.tiltangle;                              %Tilt angle
    AcquParam.Data{6} = params.job.tiltaxis;                               %Tilt axis
    AcquParam.Data{7} = params.semfields.doeuler1;                         %Detector orientation Euler 1   
    AcquParam.Data{8} = params.semfields.doeuler2;                         %Detector orientation Euler 2
    AcquParam.Data{9} = params.semfields.doeuler3;                         %Detector orientation Euler 3
    AcquParam.Data{10} = 0;                                                %Working Distance (information not available)
    AcquParam.Data{11} = 0;                                                %Insertion Distance (information not available)
elseif strcmpi(params,'Manual') %Manual prompt
    scrPrnt('SubStep','Insert microscope acquisition parameters manually');
    answer = inputdlg(strcat(AcquParam.Str,':'),'Input parameters - numeric only',...
                      [1 100],sprintfc('%d',zeros(1,11)));                 %Input dialog box
    if isempty(answer); error('Terminated by user'); end                   %Check if terminated
    AcquParam.Data = arrayfun(@str2double, answer, 'Uniform', false);   %Convert to numbers
end
%% Open ctf file
scrPrnt('Step','Opening file for writing');
filePh = fopen(fName,'w');                                                 %Open new ctf file for writing
%% Write header
scrPrnt('Step','Writing file header');
% *** Write File Info
fprintf(filePh,'Channel Text File\r\n');
fprintf(filePh,'Prj %s\r\n',fName);
fprintf(filePh,'Author\t%s\r\n',getenv('USERNAME'));
fprintf(filePh,'JobMode\tGrid\r\n');
% *** Write Grid Info
fprintf(filePh,'XCells\t%.0f\r\n',size(ebsdGrid,2));
fprintf(filePh,'YCells\t%.0f\r\n',size(ebsdGrid,1));
fprintf(filePh,'XStep\t%.4f\r\n',ebsdGrid.dx);
fprintf(filePh,'YStep\t%.4f\r\n',ebsdGrid.dy);
fprintf(filePh,'AcqE1\t%.4f\r\n',0);
fprintf(filePh,'AcqE2\t%.4f\r\n',0);
fprintf(filePh,'AcqE3\t%.4f\r\n',0);
% *** Write Acquisition parameters
fprintf(filePh,'Euler angles refer to Sample Coordinate system (CS0)!\t');	
for i = 1:length(AcquParam.Str) %Loop over aquisition parameters
   if ~strcmp(AcquParam.Fmt{i},'%s') %if numeric format is required
       AcquParam.Data{i} = num2str(AcquParam.Data{i},AcquParam.Fmt{i});    %Convert number to string
   elseif strcmp(AcquParam.Fmt{i},'%s') %if string format is required
      if ~ischar(AcquParam.Data{i}) %check if manual input is numeric
          AcquParam.Data{i} = num2str(AcquParam.Data{i}); %Convert to string
      end
   end
   fprintf(filePh,'%s\t%s\t',AcquParam.Str{i},AcquParam.Data{i});          %Write parameter
end
fprintf(filePh,'\r\n');                                                    %New line
% *** Write Phase Info
CSlst = ebsd.CSList(unique(ebsd('indexed').phaseId));                      %List of crystal systems
nrPhases = length(CSlst);                                                  %Number of phases
fprintf(filePh,'Phases\t%.0f\r\n',nrPhases);                               %Write nr of phases
for i = 1:nrPhases %Loop over phases
    mineral = CSlst{i}.mineral;                                            %Mineral name
    a = CSlst{i}.aAxis.abs;                                                %Parameter 'a'
    b = CSlst{i}.bAxis.abs;                                                %Parameter 'b'
    c = CSlst{i}.cAxis.abs;                                                %Parameter 'c'
    alpha = CSlst{i}.alpha*180/pi;                                         %Parameter 'alpha'
    beta = CSlst{i}.beta*180/pi;                                           %Parameter 'beta'
    gamma = CSlst{i}.gamma*180/pi;                                         %Parameter 'gamma'
    laueGr = find(strcmp(Laue,regexprep(CSlst{i}.pointGroup,'-','')));     %Get Laue Group
    spaceGr = 0;                                                           %Space Group (information not available)
    comment = 'Created from mtex';                                         %Phase information comment
    fprintf(filePh,'%.3f;%.3f;%.3f\t%.3f;%.3f;%.3f\t%s\t%.0f\t%.0f\t\t\t%s\r\n',...
                    a,b,c,alpha,beta,gamma,mineral,laueGr,spaceGr,comment);%Write phase info
end

%% Assemble data array
scrPrnt('Step','Assembling data array');
% *** Write Data header
fprintf(filePh,'Phase\tX\tY\tBands\tError\tEuler1\tEuler2\tEuler3\tMAD\tBC\tBS\r\n'); %Data header

%Get data order x
if ebsdGrid.x(1,1)< ebsdGrid.x(1,2)
   dim.x = 2;  
elseif ebsdGrid.x(1,1)> ebsdGrid.x(1,2)
   dim.x = -2;  
elseif ebsdGrid.x(1,1)< ebsdGrid.x(2,1)
   dim.x = 1; 
elseif ebsdGrid.x(1,1)> ebsdGrid.x(2,1)
   dim.x = -1;  
end
%Get data order y
if ebsdGrid.y(1,1)< ebsdGrid.y(1,2)
   dim.y = 2;  
elseif ebsdGrid.y(1,1)> ebsdGrid.y(1,2)
   dim.y = -2;  
elseif ebsdGrid.y(1,1)< ebsdGrid.y(2,1)
   dim.y = 1; 
elseif ebsdGrid.y(1,1)> ebsdGrid.y(2,1)
   dim.y = -1;  
end
%Gather data
flds{1} = ebsdGrid.phase;
flds{2} = ebsdGrid.x;
flds{3} = ebsdGrid.y;
flds{4} = ebsdGrid.prop.bands;
flds{5} = ebsdGrid.prop.error;
flds{6} = ebsdGrid.rotations.phi1/degree;
flds{7} = ebsdGrid.rotations.Phi/degree;
flds{8} = ebsdGrid.rotations.phi2/degree;
flds{9} = ebsdGrid.prop.mad;
flds{10} = ebsdGrid.prop.bc;
flds{11} = ebsdGrid.prop.bs;

%Write data
A = zeros(ebsdGrid.length,11); %initialize
for i = 1:length(flds)
    temp = flds{i};
    %Transpose matrices if required
    if abs(dim.x == 2) && abs(dim.y) == 1
        temp = temp';
    end
    %Flip matrices if required
    if dim.x < 0
       temp = temp(end:-1:1,:);
    end
    if dim.y < 0
       temp = temp(end:-1:1,:);
    end
    %Make vector
    A(:,i) = reshape(temp,ebsdGrid.length,1);    
end    
  
A(find(all([A(:,2)>-round0Thrsh,A(:,2)<round0Thrsh],2)),2) = 0;            %Rounding close to 0 X coordinates
A(find(all([A(:,3)>-round0Thrsh,A(:,3)<round0Thrsh],2)),3) = 0;            %Rounding close to 0 Y coordinates
A(isnan(A)) = 0;                                                           %Set NaN to 0
A(:,2) =  A(:,2) - A(1,2);                                                 %Set first x-value to 0
A(:,3) =  A(:,3) - A(1,3);                                                 %Set first y-value to 0
%% Write data array
scrPrnt('Step','Writing data array to ''ctf'' file');
outputCnt = 10;                                                            %Nr of screen output updates
k = 0;                                                                     %Counter
for i = 1:size(A,1) %Loop over all data points
   if ~mod(i,round(size(A,1)/outputCnt)) %Update screen very 'outputCnt'th iteration
        k = k+1;
        scrPrnt('SubStep',sprintf('%.0f/%.0f data lines written - %.0f percent',i,size(A,1),k*outputCnt));
   end
   fprintf(filePh,'%.0f\t%.4f\t%.4f\t%.0f\t%.0f\t%.4f\t%.4f\t%.4f\t%.4f\t%.0f\t%.0f\r\n',A(i,1),A(i,2),A(i,3),A(i,4),A(i,5),A(i,6),A(i,7),A(i,8),A(i,9),A(i,10),A(i,11)); %Write data line
end
%% Close ctf file
scrPrnt('Step','Closing file');
fclose(filePh);                                                            %Close file after writing 
scrPrnt('Step','All done',fName);    
end
%% *** Function scrPrnt - Screen Printing
function scrPrnt(mode,varargin)
%function scrPrnt(mode,varargin)
switch mode
    case 'SegmentStart'
        titleStr = varargin{1};
        fprintf('\n------------------------------------------------------');
        fprintf(['\n     ',titleStr,' \n']);
        fprintf('------------------------------------------------------\n'); 
   case 'Step'
        titleStr = varargin{1};
        fprintf([' -> ',titleStr,'\n']);
   case 'SubStep'
        titleStr = varargin{1};
        fprintf(['    - ',titleStr,'\n']);
end
end
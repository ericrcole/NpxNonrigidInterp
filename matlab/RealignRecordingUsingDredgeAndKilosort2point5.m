function RealignRecordingUsingDredgeAndKilosort2point5(basepath,config_version,dredgeParams)
% -------------------------------------------------------------------------------------------
% MODIFIED COPY. Part of the NPIX_MSIT DREDge + kriging pipeline; see matlab/README.md for what
% changed and ../CREDITS.md for provenance.
% Original: InterpolationAfterDREDge by Richard Hardstone and William Munoz (CC0 1.0),
%   https://github.com/williamunoz/InterpolationAfterDREDge
% -------------------------------------------------------------------------------------------
cd(basepath)

%% Config File
disp(['Running Kilosort with user specific settings ' config_version]);
config_string = str2func(['Kilosort25Configuration_' config_version]);
ops = config_string(dredgeParams);
ops.fbinary = 'continuous.dat';

%% GPU
gpuDeviceNum = 1;
disp(['Initializing GPU: ' num2str(gpuDeviceNum)])
gpuDevice(gpuDeviceNum); % initialize GPU (will erase any existing GPU arrays)

%% Preprocess Data
% ERC: preprocessDataSub exists here only to produce ops.fproc, which is read by
% standalone_detector inside applyDredge to draw the pre-realignment drift map.
% applyDredge's actual interpolation reads the raw file directly and never touches
% PreprocessedData.dat. So when the drift map is skipped, this whole pass over the
% recording is skipped too, and a minimal rez is built straight from the chan map.
if getOr(dredgeParams, 'skipDriftMap', false)
    disp('skipDriftMap set: bypassing preprocessDataSub, building rez from chanMap');
    [chanMap, xc, yc, kcoords, ~] = loadChanMap(ops.chanMap);
    rez.ops     = ops;
    rez.xc      = xc(:);
    rez.yc      = yc(:);
    rez.xcoords = xc(:);
    rez.ycoords = yc(:);
    rez.ops.chanMap  = chanMap;
    rez.ops.kcoords  = kcoords;
    rez.temp.Nbatch  = 0;   % overwritten by applyDredge
else
    rez = preprocessDataSub(ops);
end

%% Realign RAW file using dredge
applyDredge(rez); % last input is for shifting data

%% Delete temporary PreprocessedData.dat
if exist('PreprocessedData.dat','file')
    delete('PreprocessedData.dat')
end

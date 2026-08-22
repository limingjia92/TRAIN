function [] = aps_era5_files(orderflag_ECMWF_website)
%APS_ERA5_FILES Determine ERA5 epochs required by TRAIN and prepare downloads.
%   APS_ERA5_FILES() or APS_ERA5_FILES(0) creates the traditional TRAIN
%   ERA5_files.txt and ERA5_dates.txt inventories without submitting data.
%
%   APS_ERA5_FILES(1) additionally prepares the ERA5 download workflow.
%   The preferred path is the stabilized batch generator
%   generate_era5_train_script (Phase 3).  Until that generator is present,
%   this function safely falls back to generating one current-CDS Python
%   script per missing ERA5 epoch; the scripts are NOT launched automatically.
%
%        ****** Specific for ECMWF/CDS data website ******
%        orderflag_ECMWF_website = 0/1 - inventory only (0, default),
%                                          prepare download scripts (1)
%
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License along
%     with this program; if not, write to the Free Software Foundation, Inc.,
%     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%
% By David Bekaert - June 2017
%
% MJ 22/08/2026  Phase-2 stabilization for TRAIN fork:
%   - fix the historical primary-function mismatch (aps_era_files ->
%     aps_era5_files) while retaining the original one-argument interface;
%   - keep the dedicated era5_datapath parameter introduced in TRAIN-li;
%   - preserve the original TRAIN inventory formats: ERA5_files.txt contains
%     one filename per line with no header, and ERA5_dates.txt contains one
%     YYMMDD date per line;
%   - replace fragile duplicated char-array list construction with explicit
%     unique ERA5 datetimes, while keeping ggapYYYYMMDDHHMM.nc names;
%   - handle midnight crossover from the later interpolation epoch correctly;
%   - generate dates.txt from the actual required calendar days, including
%     next-day dates required by late-night satellite acquisitions;
%   - add ERA5_requests.txt containing exact missing YYYYMMDDHHMM epochs for
%     the stabilized batch downloader, avoiding ambiguity at midnight;
%   - retain the TRAIN-linux batch-download design, but only call the new
%     stable six-input generate_era5_train_script interface; this prevents an
%     older experimental generator from being used silently;
%   - if the stable batch generator is not yet available, generate individual
%     current-CDS scripts with aps_era5_ECMWF_Python instead of launching
%     hundreds of background requests as the historical implementation did;
%   - do not delete/overwrite existing NetCDF files automatically;
%   - validate UTC_sat, data ranges, output paths, and file creation errors;
%   - use floor/ceil region bounds plus the historical 2-degree margin so the
%     requested ERA5 area never becomes smaller than the intended buffer.

% -------------------------------------------------------------------------
% Input and TRAIN parameters
% -------------------------------------------------------------------------
if nargin<1 || isempty(orderflag_ECMWF_website)
    orderflag_ECMWF_website = 0;
end
if ~isnumeric(orderflag_ECMWF_website) && ~islogical(orderflag_ECMWF_website)
    error('aps_era5_files:InvalidOrderFlag', ...
        'orderflag_ECMWF_website must be 0 or 1.');
end
if ~isscalar(orderflag_ECMWF_website) || ...
        ~(orderflag_ECMWF_website==0 || orderflag_ECMWF_website==1)
    error('aps_era5_files:InvalidOrderFlag', ...
        'orderflag_ECMWF_website must be 0 or 1.');
end

workdir = pwd;
stamps_processed = getparm_aps('stamps_processed');
UTC_sat = getparm_aps('UTC_sat');
era5_datapath = getparm_aps('era5_datapath');

if isempty(era5_datapath)
    error('aps_era5_files:MissingDataPath', ...
        'Please specify era5_datapath in the TRAIN parameters.');
end
if isa(era5_datapath,'string')
    era5_datapath = char(era5_datapath);
end
if ~ischar(era5_datapath)
    error('aps_era5_files:InvalidDataPath', ...
        'era5_datapath must be a character/string path.');
end

% -------------------------------------------------------------------------
% Load acquisition dates using the original TRAIN/Stamps logic
% -------------------------------------------------------------------------
if strcmpi(stamps_processed,'y')
    ll_matfile = getparm_aps('ll_matfile');
    ps = load(ll_matfile);
    if ~isfield(ps,'day') || isempty(ps.day)
        error('aps_era5_files:MissingDates', ...
            'The StaMPS file %s does not contain a non-empty day variable.', ...
            ll_matfile);
    end
    dates = double(ps.day(:));
else
    ifgday_matfile = getparm_aps('ifgday_matfile');
    ifgs_dates = load(ifgday_matfile);
    if ~isfield(ifgs_dates,'ifgday') || isempty(ifgs_dates.ifgday)
        error('aps_era5_files:MissingDates', ...
            'The file %s does not contain a non-empty ifgday variable.', ...
            ifgday_matfile);
    end
    dates_yyyymmdd = unique(double(ifgs_dates.ifgday(:)));
    dates = datenum(num2str(dates_yyyymmdd),'yyyymmdd');
end

dates = unique(floor(dates(isfinite(dates))));
if isempty(dates)
    error('aps_era5_files:MissingDates','No valid acquisition dates were found.');
end

% -------------------------------------------------------------------------
% Determine the two hourly ERA5 epochs bracketing UTC_sat
% -------------------------------------------------------------------------
utc_hour = parse_utc_sat(UTC_sat);
t_before = floor(utc_hour);
t_after  = ceil(utc_hour);

before_hour = mod(t_before,24);
after_hour  = mod(t_after,24);
time1 = sprintf('%02d00',before_hour);
time2 = sprintf('%02d00',after_hour);

fprintf('Satellite pass is %.6g UTC. ERA5 bracketing epochs: %s and %s.\n', ...
    utc_hour,time1,time2);

% Build exact required datetimes. Using date + hour/24 handles 24:00 as the
% next calendar day automatically. Exact-integer UTC naturally collapses to
% one unique ERA5 epoch per acquisition.
required_dt = zeros(numel(dates)*2,1);
for k = 1:numel(dates)
    required_dt(2*k-1) = dates(k) + t_before/24;
    required_dt(2*k)   = dates(k) + t_after/24;
end
required_dt = unique(required_dt);

n_required = numel(required_dt);
filelist = cell(n_required,1);
request_ids = cell(n_required,1);
for k = 1:n_required
    request_ids{k} = datestr(required_dt(k),'yyyymmddHHMM');
    filelist{k} = ['ggap' request_ids{k} '.nc'];
end

required_days = unique(floor(required_dt));

% -------------------------------------------------------------------------
% Preserve traditional TRAIN inventory files
% -------------------------------------------------------------------------
write_text_lines(fullfile(workdir,'ERA5_files.txt'),filelist);

legacy_dates = cell(numel(required_days),1);
for k = 1:numel(required_days)
    legacy_dates{k} = datestr(required_days(k),'yymmdd');
end
write_text_lines(fullfile(workdir,'ERA5_dates.txt'),legacy_dates);

fprintf('Required ERA5 files for all interferograms: %d\n',n_required);
for k = 1:n_required
    fprintf('%s\n',filelist{k});
end
fprintf('Generated inventory: %s\n',fullfile(workdir,'ERA5_files.txt'));
fprintf('Generated legacy date list: %s\n',fullfile(workdir,'ERA5_dates.txt'));

if orderflag_ECMWF_website==0
    return
end

% -------------------------------------------------------------------------
% Prepare downloads: never overwrite existing valid-path files implicitly
% -------------------------------------------------------------------------
if exist(era5_datapath,'dir')~=7
    [status,msg] = mkdir(era5_datapath);
    if ~status
        error('aps_era5_files:CreateDirectoryFailed', ...
            'Cannot create era5_datapath %s: %s',era5_datapath,msg);
    end
end

missing_ids = {};
for k = 1:n_required
    day_dir = fullfile(era5_datapath,request_ids{k}(1:8));
    nc_path = fullfile(day_dir,filelist{k});
    if exist(nc_path,'file')~=2
        missing_ids{end+1,1} = request_ids{k}; %#ok<AGROW>
    end
end

if isempty(missing_ids)
    fprintf('All %d required ERA5 files already exist. No download scripts needed.\n', ...
        n_required);
    return
end

fprintf('%d of %d required ERA5 files are missing.\n', ...
    numel(missing_ids),n_required);

% dates.txt is retained for compatibility with the TRAIN-linux batch idea.
% It contains all calendar days that have at least one missing exact epoch,
% including a next day required by a 23:xx satellite pass.
missing_days_num = zeros(numel(missing_ids),1);
for k = 1:numel(missing_ids)
    missing_days_num(k) = datenum(missing_ids{k}(1:8),'yyyymmdd');
end
missing_days_num = unique(missing_days_num);
missing_days = cell(numel(missing_days_num),1);
for k = 1:numel(missing_days_num)
    missing_days{k} = datestr(missing_days_num(k),'yyyymmdd');
end

dates_filepath = fullfile(era5_datapath,'dates.txt');
write_text_lines(dates_filepath,missing_days);

% Exact request manifest for the stabilized Phase-3 batch downloader.
requests_filepath = fullfile(era5_datapath,'ERA5_requests.txt');
write_text_lines(requests_filepath,missing_ids);

fprintf('Generated batch date list: %s\n',dates_filepath);
fprintf('Generated exact request manifest: %s\n',requests_filepath);

% -------------------------------------------------------------------------
% Region for CDS request (historical TRAIN margin = 2 degrees)
% -------------------------------------------------------------------------
region_lat_range = getparm_aps('region_lat_range');
region_lon_range = getparm_aps('region_lon_range');
if isempty(region_lat_range) || isempty(region_lon_range)
    error('aps_era5_files:MissingRegion', ...
        'Specify both region_lat_range and region_lon_range for ERA5 downloads.');
end

region_lat_range = double(region_lat_range(:));
region_lon_range = double(region_lon_range(:));
if any(~isfinite(region_lat_range)) || any(~isfinite(region_lon_range))
    error('aps_era5_files:InvalidRegion', ...
        'ERA5 latitude/longitude ranges must contain finite numeric values.');
end

crop_range_in = 2;
S = max(-90,floor(min(region_lat_range)) - crop_range_in);
N = min( 90, ceil(max(region_lat_range)) + crop_range_in);
W = floor(min(region_lon_range)) - crop_range_in;
E = ceil(max(region_lon_range)) + crop_range_in;

weatherregstr = sprintf('%.10g/%.10g/%.10g/%.10g',N,W,S,E);
fprintf('ERA5 request region N/W/S/E: %s\n',weatherregstr);

% -------------------------------------------------------------------------
% Preferred Phase-3 stable batch generator; safe Phase-2 fallback otherwise
% -------------------------------------------------------------------------
use_stable_batch = false;
if exist('generate_era5_train_script','file')==2
    try
        generator_nargin = nargin('generate_era5_train_script');
        % Stable Phase-3 interface planned as:
        % generate_era5_train_script(output_dir,N,W,S,E,request_manifest)
        use_stable_batch = (generator_nargin==6);
    catch
        use_stable_batch = false;
    end
end

if use_stable_batch
    fprintf('Generating stabilized ERA5 batch downloader...\n');
    generate_era5_train_script(era5_datapath,N,W,S,E,requests_filepath);
    fprintf('Batch downloader prepared in: %s\n',era5_datapath);
    fprintf('Run it manually after reviewing the request: python download_ERA5_train.py\n');
else
    if exist('generate_era5_train_script','file')==2
        warning('aps_era5_files:LegacyBatchGeneratorIgnored', ...
            ['A generate_era5_train_script implementation was found, but it does ' ...
             'not expose the stabilized six-input interface. It is ignored to ' ...
             'avoid silently using the older experimental batch generator.']);
    else
        fprintf(['Stable batch generator is not present yet (Phase 3). ' ...
            'Generating individual CDS scripts as a safe fallback.\n']);
    end

    n_scripts = 0;
    for k = 1:numel(missing_ids)
        timing = missing_ids{k};
        subdirpath = fullfile(era5_datapath,timing(1:8));
        if exist(subdirpath,'dir')~=7
            [status,msg] = mkdir(subdirpath);
            if ~status
                error('aps_era5_files:CreateDirectoryFailed', ...
                    'Cannot create directory %s: %s',subdirpath,msg);
            end
        end
        generate_single_epoch_script(subdirpath,timing,weatherregstr);
        n_scripts = n_scripts + 1;
    end

    fprintf(['Generated %d individual ERA5 CDS Python scripts. They were NOT ' ...
        'started automatically; this avoids flooding the CDS queue.\n'],n_scripts);
end

end


function utc_hour = parse_utc_sat(UTC_sat)
% Convert historical TRAIN UTC_sat representations to decimal UTC hours.
if isnumeric(UTC_sat)
    if ~isscalar(UTC_sat) || ~isfinite(UTC_sat)
        error('aps_era5_files:InvalidUTC', ...
            'Numeric UTC_sat must be one finite decimal-hour value.');
    end
    utc_hour = double(UTC_sat);
else
    if isa(UTC_sat,'string')
        if ~isscalar(UTC_sat)
            error('aps_era5_files:InvalidUTC','UTC_sat must be scalar.');
        end
        UTC_sat = char(UTC_sat);
    end
    if ~ischar(UTC_sat)
        error('aps_era5_files:InvalidUTC', ...
            'UTC_sat must be numeric or a character/string time.');
    end

    value = strtrim(UTC_sat);
    if isempty(value)
        error('aps_era5_files:InvalidUTC','UTC_sat is empty.');
    end

    if ~isempty(strfind(value,':')) %#ok<STREMP>
        fields = strsplit(value,':');
        if numel(fields)<2 || numel(fields)>3
            error('aps_era5_files:InvalidUTC', ...
                'UTC_sat must have HH:MM or HH:MM:SS form.');
        end
        hh = str2double(fields{1});
        mm = str2double(fields{2});
        ss = 0;
        if numel(fields)==3
            ss = str2double(fields{3});
        end
    elseif length(value)==4 && all(value>='0' & value<='9')
        hh = str2double(value(1:2));
        mm = str2double(value(3:4));
        ss = 0;
    else
        error('aps_era5_files:InvalidUTC', ...
            'UTC_sat must be HH:MM, HH:MM:SS, HHMM, or decimal hours.');
    end

    if any(~isfinite([hh mm ss])) || hh<0 || hh>=24 || ...
            mm<0 || mm>=60 || ss<0 || ss>=60
        error('aps_era5_files:InvalidUTC','Invalid UTC_sat value: %s',value);
    end
    utc_hour = hh + mm/60 + ss/3600;
end

if utc_hour<0 || utc_hour>=24
    error('aps_era5_files:InvalidUTC', ...
        'UTC_sat decimal hours must satisfy 0 <= UTC_sat < 24.');
end
end


function write_text_lines(filepath,lines)
fid = fopen(filepath,'w');
if fid<0
    error('aps_era5_files:FileOpenFailed','Cannot create file: %s',filepath);
end
try
    for k = 1:numel(lines)
        fprintf(fid,'%s\n',lines{k});
    end
catch ME
    fclose(fid);
    rethrow(ME);
end
fclose(fid);
end


function generate_single_epoch_script(subdirpath,timing,weatherregstr)
% Keep the historical aps_era5_ECMWF_Python calling convention while
% containing the temporary directory change locally and restoring it even
% if script generation fails.
old_dir = pwd;
cleanup_obj = onCleanup(@() cd(old_dir)); %#ok<NASGU>
cd(subdirpath);
aps_era5_ECMWF_Python(timing,weatherregstr);
end

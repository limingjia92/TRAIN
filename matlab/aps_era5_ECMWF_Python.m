function [] = aps_era5_ECMWF_Python(timing,weatherregstr)
%APS_ERA5_ECMWF_PYTHON Create a CDS API Python script for one ERA5 epoch.
%   APS_ERA5_ECMWF_PYTHON(TIMING,WEATHERREGSTR) writes a Python download
%   script in the current directory. TIMING may be 'YYYYMMDDHH' or
%   'YYYYMMDDHHMM'; ERA5 pressure-level data are hourly, so MM must be 00.
%   WEATHERREGSTR is the TRAIN region string 'N/W/S/E'.
%
%   The generated file is named <YYYYMMDDHHMM>.py and downloads
%   ggapYYYYMMDDHHMM.nc, preserving the filename/layout expected by TRAIN.
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
% By Bekaert Davide -- June 2017
% Modified by Marin Govorcin -- May 2019 to use CDS API
% Modified by Jianbao Sun -- September 2025 for regional ERA5 downloads
%
% MJ 22/08/2026  Phase-2 stabilization for TRAIN fork:
%   - keep the historical two-input TRAIN calling interface;
%   - accept both YYYYMMDDHH and YYYYMMDDHHMM timing strings and validate
%     that the requested ERA5 epoch is an integer hour;
%   - replace the obsolete CDS request key "format" with the current
%     "data_format" request syntax and request an unarchived NetCDF file;
%   - generate the current hourly time form HH:00 and retain the full set of
%     37 ERA5 pressure levels required by TRAIN;
%   - validate the N/W/S/E region instead of relying on direction-letter
%     parsing of numeric strings;
%   - remove global urllib3 SSL-warning suppression: certificate/proxy
%     problems should remain visible rather than being hidden;
%   - update the generated script to use cdsapi.Client() with the user's
%     current CDS credentials (.cdsapirc);
%   - write to a temporary .part target and atomically rename it only after
%     a successful download, reducing the chance of treating a partial file
%     as a valid TRAIN input;
%   - add fopen/error checks and remove unused legacy variables/comments.

% -------------------------------------------------------------------------
% Validate and normalize timing
% -------------------------------------------------------------------------
if isnumeric(timing)
    if ~isscalar(timing) || ~isfinite(timing)
        error('aps_era5_ECMWF_Python:InvalidTiming', ...
            'timing must be a scalar number or a character/string value.');
    end
    timing = sprintf('%.0f',timing);
elseif isa(timing,'string')
    if ~isscalar(timing)
        error('aps_era5_ECMWF_Python:InvalidTiming', ...
            'timing must be a scalar string.');
    end
    timing = char(timing);
elseif ~ischar(timing)
    error('aps_era5_ECMWF_Python:InvalidTiming', ...
        'timing must be ''YYYYMMDDHH'' or ''YYYYMMDDHHMM''.');
end

timing = strtrim(timing);
if length(timing)==10
    timing = [timing '00'];
elseif length(timing)~=12
    error('aps_era5_ECMWF_Python:InvalidTiming', ...
        'timing must contain 10 or 12 digits: YYYYMMDDHH[MM].');
end

if any(timing<'0' | timing>'9')
    error('aps_era5_ECMWF_Python:InvalidTiming', ...
        'timing must contain digits only: YYYYMMDDHH[MM].');
end

year_str   = timing(1:4);
month_str  = timing(5:6);
day_str    = timing(7:8);
hour_str   = timing(9:10);
minute_str = timing(11:12);

month_num = str2double(month_str);
day_num   = str2double(day_str);
hour_num  = str2double(hour_str);

if month_num<1 || month_num>12 || day_num<1 || day_num>31 || ...
        hour_num<0 || hour_num>23
    error('aps_era5_ECMWF_Python:InvalidTiming', ...
        'Invalid date/hour fields in timing %s.',timing);
end
if ~strcmp(minute_str,'00')
    error('aps_era5_ECMWF_Python:NonHourlyTiming', ...
        'ERA5 pressure-level data are hourly; timing minute must be 00.');
end

% -------------------------------------------------------------------------
% Parse and validate TRAIN weather region: N/W/S/E
% -------------------------------------------------------------------------
if isa(weatherregstr,'string')
    if ~isscalar(weatherregstr)
        error('aps_era5_ECMWF_Python:InvalidRegion', ...
            'weatherregstr must be one N/W/S/E string.');
    end
    weatherregstr = char(weatherregstr);
end
if ~ischar(weatherregstr)
    error('aps_era5_ECMWF_Python:InvalidRegion', ...
        'weatherregstr must be a character/string value in N/W/S/E form.');
end

parts = strsplit(strtrim(weatherregstr),'/');
if numel(parts)~=4
    error('aps_era5_ECMWF_Python:InvalidRegion', ...
        'weatherregstr must contain exactly four values: N/W/S/E.');
end

coordinates = zeros(1,4);
for k = 1:4
    coordinates(k) = str2double(strtrim(parts{k}));
end
if any(~isfinite(coordinates))
    error('aps_era5_ECMWF_Python:InvalidRegion', ...
        'Unable to parse finite numeric N/W/S/E coordinates from "%s".', ...
        weatherregstr);
end

N = coordinates(1);
W = coordinates(2);
S = coordinates(3);
E = coordinates(4);

if N<S || N>90 || S<-90
    error('aps_era5_ECMWF_Python:InvalidRegion', ...
        'Invalid latitude bounds in N/W/S/E region: %s.',weatherregstr);
end

% -------------------------------------------------------------------------
% Generate current CDS API Python request
% -------------------------------------------------------------------------
target = ['ggap' timing '.nc'];
pythonsc_path = [timing '.py'];

fid = fopen(pythonsc_path,'w');
if fid<0
    error('aps_era5_ECMWF_Python:FileOpenFailed', ...
        'Cannot create Python script: %s',pythonsc_path);
end

try
    fprintf(fid,'#!/usr/bin/env python3\n');
    fprintf(fid,'\n');
    fprintf(fid,'import os\n');
    fprintf(fid,'import cdsapi\n');
    fprintf(fid,'\n');
    fprintf(fid,'# Configure current CDS credentials in ~/.cdsapirc before running.\n');
    fprintf(fid,'dataset = "reanalysis-era5-pressure-levels"\n');
    fprintf(fid,'request = {\n');
    fprintf(fid,'    "product_type": ["reanalysis"],\n');
    fprintf(fid,'    "variable": ["geopotential", "relative_humidity", "temperature"],\n');
    fprintf(fid,'    "pressure_level": [\n');
    fprintf(fid,'        "1", "2", "3", "5", "7", "10", "20", "30", "50", "70",\n');
    fprintf(fid,'        "100", "125", "150", "175", "200", "225", "250", "300",\n');
    fprintf(fid,'        "350", "400", "450", "500", "550", "600", "650", "700",\n');
    fprintf(fid,'        "750", "775", "800", "825", "850", "875", "900", "925",\n');
    fprintf(fid,'        "950", "975", "1000"\n');
    fprintf(fid,'    ],\n');
    fprintf(fid,'    "year": ["%s"],\n',year_str);
    fprintf(fid,'    "month": ["%s"],\n',month_str);
    fprintf(fid,'    "day": ["%s"],\n',day_str);
    fprintf(fid,'    "time": ["%s:00"],\n',hour_str);
    fprintf(fid,'    "area": [%.10g, %.10g, %.10g, %.10g],  # N, W, S, E\n',N,W,S,E);
    fprintf(fid,'    "data_format": "netcdf",\n');
    fprintf(fid,'    "download_format": "unarchived",\n');
    fprintf(fid,'}\n');
    fprintf(fid,'\n');
    fprintf(fid,'target = r"%s"\n',target);
    fprintf(fid,'partial_target = target + ".part"\n');
    fprintf(fid,'\n');
    fprintf(fid,'if os.path.exists(partial_target):\n');
    fprintf(fid,'    os.remove(partial_target)\n');
    fprintf(fid,'\n');
    fprintf(fid,'client = cdsapi.Client()\n');
    fprintf(fid,'client.retrieve(dataset, request, partial_target)\n');
    fprintf(fid,'os.replace(partial_target, target)\n');
    fprintf(fid,'print(f"ERA5 download complete: {target}")\n');
catch ME
    fclose(fid);
    if exist(pythonsc_path,'file')==2
        delete(pythonsc_path);
    end
    rethrow(ME);
end

fclose(fid);
fprintf('ERA5 CDS Python script generated: %s\n',pythonsc_path);

end

function [Temp,WVapour,Geopot,Pressure,longrid,latgrid,xx,yy,lon0360_flag] = aps_load_era(file,era_data_type)
% loading ERA-I or ERA5 data from ECMWF website or ERA-I from BADC website
% Bekaert David
% modifications
% DB    10/04/2016  extract code from aps_era_SAR.m to make code modular
% DB    07/06/2017  Update syntax to include ERA5 model
% MJ    22/08/2024  Adapt to new ERA5 NetCDF coordinate names
% MJ    30/08/2024  Add compatibility attempt for newer NetCDF files
% MJ    21/08/2026  Phase-1 stabilization for TRAIN fork:
%                   - remove global warning('off') side effect
%                   - keep low-level netcdf.getVar reading so scale_factor/
%                     add_offset are applied exactly once
%                   - read only required ECMWF/ERA5 variables and support
%                     legacy/current aliases (level/pressure_level, etc.)
%                   - ignore singleton time-like dimensions safely and
%                     merge expver slices by taking the first finite value
%                   - derive pressure-level count dynamically (no hard-coded 37)
%                   - explicitly order pressure levels from high to low
%                     (surface to upper atmosphere), preserving TRAIN convention
%                   - restore the original BADC latitude/longitude loading path
%                   - fail explicitly on unsupported non-singleton dimensions
%
% Output arrays Temp, WVapour, Geopot and Pressure are [lat, lon, level].
% Pressure levels are ordered from high pressure to low pressure.

% debug figure to test and validate data loading
debug_fig = 0;

% open the NetCDF file and guarantee closure if an error occurs
ncid = netcdf.open(file,'NC_NOWRITE');
cleanupObj = onCleanup(@() close_netcdf_safely(ncid));

%% Swap between BADC and ECMWF website data
if strcmpi(era_data_type,'ECMWF')
    % Read only the variables required by TRAIN.  This avoids failures on
    % auxiliary NC_STRING variables present in some newer CDS NetCDF files.
    longitude = read_coordinate(ncid, {'longitude','lon'});
    latitude  = read_coordinate(ncid, {'latitude','lat'});
    Plevs     = read_coordinate(ncid, {'pressure_level','level','isobaricinhpa'});

    Temp    = read_ecmwf_field(ncid, {'t','temperature'}, ...
                               numel(latitude),numel(longitude),numel(Plevs));
    Hum     = read_ecmwf_field(ncid, {'r','relative_humidity','relativehumidity'}, ...
                               numel(latitude),numel(longitude),numel(Plevs));
    Geopot  = read_ecmwf_field(ncid, {'z','geopotential'}, ...
                               numel(latitude),numel(longitude),numel(Plevs));

    % TRAIN historically expects level index 1 near the surface and the last
    % index in the upper atmosphere.  Do not assume the order provided by CDS.
    [Plevs,level_ix] = sort(double(Plevs(:)),'descend');
    if numel(unique(Plevs)) ~= numel(Plevs)
        error('aps_load_era:DuplicatePressureLevels', ...
              'Duplicate pressure levels found in %s.',file);
    end
    Temp   = Temp(:,:,level_ix);
    Hum    = Hum(:,:,level_ix);
    Geopot = Geopot(:,:,level_ix);

elseif strcmpi(era_data_type,'BADC')
    % Preserve the original TRAIN BADC layout.  The historical BADC reader
    % identifies these variables by ID, so keep that behavior unchanged.
    longitude = double(netcdf.getVar(ncid,0));
    latitude  = double(netcdf.getVar(ncid,1));
    Plevs     = double(netcdf.getVar(ncid,2));
    Geopot    = double(netcdf.getVar(ncid,4));
    Temp      = double(netcdf.getVar(ncid,5));
    Hum       = double(netcdf.getVar(ncid,11));

    Temp   = permute(Temp,[2,1,3]);
    Hum    = permute(Hum,[2,1,3]);
    Geopot = permute(Geopot,[2,1,3]);

    longitude = double(longitude(:));
    latitude  = double(latitude(:));
    Plevs     = double(Plevs(:));

else
    error('aps_load_era:UnsupportedDataType', ...
          'Unsupported era_data_type: %s',era_data_type);
end

% Close the NetCDF file now that all required variables have been loaded.
clear cleanupObj

%% Validate dimensions and build coordinate grids dynamically
longitude = double(longitude(:));
latitude  = double(latitude(:));
Plevs     = double(Plevs(:));

n_latitude_points  = numel(latitude);
n_longitude_points = numel(longitude);
n_pressure_levels  = numel(Plevs);

expected_size = [n_latitude_points,n_longitude_points,n_pressure_levels];
if ~isequal(size3(Temp),expected_size) || ...
   ~isequal(size3(Hum),expected_size) || ...
   ~isequal(size3(Geopot),expected_size)
    error('aps_load_era:DimensionMismatch', ...
          ['Weather-model field dimensions do not match latitude/longitude/' ...
           'pressure coordinates in %s.'],file);
end

Pressure = repmat(reshape(Plevs,1,1,n_pressure_levels), ...
                  n_latitude_points,n_longitude_points,1);
latgrid  = repmat(reshape(latitude,n_latitude_points,1,1), ...
                  1,n_longitude_points,n_pressure_levels);
longrid  = repmat(reshape(longitude,1,n_longitude_points,1), ...
                  n_latitude_points,1,n_pressure_levels);

% Get list of points to look at in analysis
[xx,yy] = meshgrid(1:n_longitude_points,1:n_latitude_points);

%% Calculate water-vapour pressure
% (see IFS documentation part 2: Data assimilation (CY25R1)).
% Calculate saturated water vapour pressure (svp) for water (svpw) using
% Buck (1981) and for ice (svpi) from Alduchow and Eskridge (1996).
svpw = 6.1121.*exp((17.502.*(Temp-273.16))./(240.97+Temp-273.16));
svpi = 6.1121.*exp((22.587.*(Temp-273.16))./(273.86+Temp-273.16));
tempbound1 = 273.16; % 0 C
tempbound2 = 250.16; % -23 C

wgt = (Temp-tempbound2)/(tempbound1-tempbound2);
svp = svpi + (svpw-svpi).*wgt.^2;
ix_bound1 = Temp > tempbound1;
svp(ix_bound1) = svpw(ix_bound1);
ix_bound2 = Temp < tempbound2;
svp(ix_bound2) = svpi(ix_bound2);
WVapour = Hum./100.*svp;
clear Hum

% Inform about the organisation of longitudes
if any(longitude > 180)
    lon0360_flag = 'y';
else
    lon0360_flag = 'n';
end

%% Validation plots
if debug_fig == 1
    figure('position',[3 628 1402 586]);
    subplot(2,5,1)
    imagesc(Temp(:,:,end)); colorbar; axis xy; axis equal; axis tight
    title('temp upper atmo')
    subplot(2,5,2)
    imagesc(Temp(:,:,1)); colorbar; axis xy; axis equal; axis tight
    title('temp lower atmo')

    subplot(2,5,3)
    imagesc(Pressure(:,:,end)); colorbar; axis xy; axis equal; axis tight
    title('pressure upper atmo')
    subplot(2,5,4)
    imagesc(Pressure(:,:,1)); colorbar; axis xy; axis equal; axis tight
    title('pressure lower atmo')

    subplot(2,5,6)
    imagesc(WVapour(:,:,end)); colorbar; axis xy; axis equal; axis tight
    title('Water vapour upper atmo')
    subplot(2,5,7)
    imagesc(WVapour(:,:,1)); colorbar; axis xy; axis equal; axis tight
    title('Water vapour lower atmo')

    subplot(2,5,8)
    imagesc(Geopot(:,:,end)); colorbar; axis xy; axis equal; axis tight
    title('geopotential upper atmo')
    subplot(2,5,9)
    imagesc(Geopot(:,:,1)); colorbar; axis xy; axis equal; axis tight
    title('geopotential lower atmo')

    subplot(2,5,5)
    imagesc(latgrid(:,:,1)); colorbar; axis xy; axis equal; axis tight
    title('lat')
    subplot(2,5,10)
    imagesc(longrid(:,:,1)); colorbar; axis xy; axis equal; axis tight
    title('lon')
end

end


function coord = read_coordinate(ncid,aliases)
% Read a numeric coordinate variable using a list of accepted aliases.
[data,~,~,~] = read_scaled_variable(ncid,aliases);
if ~isvector(data)
    error('aps_load_era:CoordinateNotVector', ...
          'Coordinate variable %s is not one-dimensional.',aliases{1});
end
coord = double(data(:));
end


function field = read_ecmwf_field(ncid,aliases,nlat,nlon,nlev)
% Read a packed ECMWF/ERA5 field and return it as [lat, lon, level].
[data,varname,dim_names,dim_lengths] = read_scaled_variable(ncid,aliases);

% Identify the three physical dimensions.  netcdf.inqVar returns dimension
% IDs in MATLAB/Fortran order, which matches the order of netcdf.getVar data.
is_lon = cellfun(@(x) is_alias(x,{'longitude','lon'}),dim_names);
is_lat = cellfun(@(x) is_alias(x,{'latitude','lat'}),dim_names);
is_lev = cellfun(@(x) is_alias(x,{'pressure_level','level','isobaricinhpa'}),dim_names);

if sum(is_lon)~=1 || sum(is_lat)~=1 || sum(is_lev)~=1
    error('aps_load_era:UnknownDimensions', ...
          'Cannot identify lon/lat/pressure dimensions for variable %s.',varname);
end

% Remove dimensions not used by TRAIN.  Singleton time/number dimensions are
% harmless.  expver can contain multiple complementary slices; merge those
% by taking the first finite value at each grid point.
target_mask = is_lon | is_lat | is_lev;
extra_ix = find(~target_mask);
for k = numel(extra_ix):-1:1
    d = extra_ix(k);
    dname = lower(dim_names{d});
    dlen = dim_lengths(d);

    if dlen == 1
        data = select_dimension(data,d,1,numel(dim_lengths));
    elseif strcmp(dname,'expver')
        data = merge_expver_dimension(data,d,dlen,numel(dim_lengths));
    else
        error('aps_load_era:UnsupportedDimension', ...
              ['Variable %s contains unsupported non-singleton dimension ' ...
               '%s (length %d). TRAIN expects one time step per file.'], ...
              varname,dim_names{d},dlen);
    end

    dim_names(d) = [];
    dim_lengths(d) = [];
    data = reshape_to_dimensions(data,dim_lengths);
end

% Re-identify positions after removing auxiliary dimensions.
lon_ix = find(cellfun(@(x) is_alias(x,{'longitude','lon'}),dim_names));
lat_ix = find(cellfun(@(x) is_alias(x,{'latitude','lat'}),dim_names));
lev_ix = find(cellfun(@(x) is_alias(x,{'pressure_level','level','isobaricinhpa'}),dim_names));

if numel(dim_names)~=3 || isempty(lon_ix) || isempty(lat_ix) || isempty(lev_ix)
    error('aps_load_era:UnexpectedFieldLayout', ...
          'Unexpected dimension layout for variable %s.',varname);
end

field = permute(data,[lat_ix,lon_ix,lev_ix]);
field = reshape(field,nlat,nlon,nlev);
end


function [data,varname,dim_names,dim_lengths] = read_scaled_variable(ncid,aliases)
% Read one required numeric variable with low-level NetCDF calls.
% Unlike ncread, netcdf.getVar does not apply CF scale_factor/add_offset;
% therefore those attributes are applied explicitly exactly once here.
[varid,varname] = find_variable(ncid,aliases);
[~,~,dimids,numatts] = netcdf.inqVar(ncid,varid);

raw = netcdf.getVar(ncid,varid);
if ~isnumeric(raw)
    error('aps_load_era:NonNumericVariable', ...
          'Required variable %s is not numeric.',varname);
end
raw = double(raw);

scale = 1;
offset = 0;
fill_values = [];
for j = 0:numatts-1
    attname = netcdf.inqAttName(ncid,varid,j);
    attval = netcdf.getAtt(ncid,varid,attname);
    if ~isnumeric(attval)
        continue
    end

    switch lower(attname)
        case 'scale_factor'
            scale = double(attval(1));
        case 'add_offset'
            offset = double(attval(1));
        case {'_fillvalue','missing_value'}
            fill_values = [fill_values; double(attval(:))]; %#ok<AGROW>
    end
end

for j = 1:numel(fill_values)
    raw(raw == fill_values(j)) = NaN;
end
data = raw.*scale + offset;

% netcdf.inqVar returns dimids in MATLAB/Fortran order, i.e. the same order
% as the dimensions of data returned by netcdf.getVar.
dim_names = cell(1,numel(dimids));
dim_lengths = zeros(1,numel(dimids));
for j = 1:numel(dimids)
    [dim_names{j},dim_lengths(j)] = netcdf.inqDim(ncid,dimids(j));
end

data = reshape_to_dimensions(data,dim_lengths);
end


function [varid,varname] = find_variable(ncid,aliases)
% Case-insensitive variable-name lookup, preserving alias priority.
[~,numvars,~,~] = netcdf.inq(ncid);
all_names = cell(1,numvars);
for i = 0:numvars-1
    [all_names{i+1},~,~,~] = netcdf.inqVar(ncid,i);
end

varid = [];
varname = '';
for a = 1:numel(aliases)
    ix = find(strcmpi(all_names,aliases{a}),1,'first');
    if ~isempty(ix)
        varid = ix-1;
        varname = all_names{ix};
        return
    end
end

error('aps_load_era:MissingVariable', ...
      'Required NetCDF variable not found. Accepted names: %s', ...
      join_aliases(aliases));
end


function tf = is_alias(name,aliases)
tf = any(strcmpi(name,aliases));
end


function data = select_dimension(data,dim,index,nd)
subs = repmat({':'},1,nd);
subs{dim} = index;
data = data(subs{:});
end


function data = merge_expver_dimension(data,dim,nexpver,nd)
% ERA5 files can contain expver=1 and expver=5 slices with complementary
% valid values.  Keep the first finite value at each location.
subs = repmat({':'},1,nd);
subs{dim} = 1;
data_out = data(subs{:});
for k = 2:nexpver
    subs{dim} = k;
    data_k = data(subs{:});
    ix = isnan(data_out) & ~isnan(data_k);
    data_out(ix) = data_k(ix);
end
data = data_out;
end


function data = reshape_to_dimensions(data,dim_lengths)
if isempty(dim_lengths)
    data = data(1);
elseif numel(dim_lengths) == 1
    data = reshape(data,[dim_lengths(1),1]);
else
    data = reshape(data,dim_lengths);
end
end


function sz = size3(data)
sz = [size(data,1),size(data,2),size(data,3)];
end


function txt = join_aliases(aliases)
txt = aliases{1};
for k = 2:numel(aliases)
    txt = [txt ', ' aliases{k}]; %#ok<AGROW>
end
end


function close_netcdf_safely(ncid)
try
    netcdf.close(ncid);
catch
    % Do not mask the original error during cleanup.
end
end

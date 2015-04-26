function [ data, properties ] = load_bbdata( bbfile, t0, duration )
%
% load_bbdata()  - allocate a data vector and load a VFBB data file into it
% merge_bbdata() - merge a VFBB data file into an existing data vector
%
% [ data, properties ] = load_bbdata(  bbfile, t0,   duration   )
%   merged_data        = merge_bbdata( bbfile, data, properties )
%
% input arguments:
%
%    bbfile      is the (string) name of a VFBB data file
%    t0          is an array encoding the start time of the returned data
%                vector: [ year day hour minute second msec usec ]
%    duration    is the (scalar) length of the returned data vector, in seconds
%    data        is a data vector created by load_bbdata() or a merged_data
%                vector returned by merge_bbdata()
%    properties  is the properties for the data vector input argument
%
% output arguments:
%
%    data        is the data vector created by load_bbdata(), filled with any
%                data from the bbfile that fit within the time range of the
%                data vector; empty data elements contain NaNs
%    properties  is created by load_bbdata() to describe the properties of the
%                data vector to merge_bbdata()
%    merged_data is the input data vector merged with any data from the bbfile
%                that fit within the time range of the input data vector; empty
%                data elements contain NaNs
%
% The elements of t0 are optional; missing values are set to 0.  If all the
% values are missing, i.e., the t0 array is empty ([]), t0 is set to the time
% of the first data sample from the VFBB data file integer header.
%
% duration must be positive.  duration is truncated (not rounded) to an integral
% number of sample intervals.
%
% load_bbdata() creates an empty data vector from the specified start time and
% duration, then loads VFBB data from the specified bbfile that fit within the
% time range of the data vector.  Additional VFBB data can then be merged into
% the same data vector with one or more subsequent calls to merge_bbdata().  The
% time range of the merged data vector is exactly the same as the time range of
% the initial data vector.  That is, the start time and duration of a data vec-
% tor are fixed by the call to load_bbdata() and are not altered by any subse-
% quent calls to merge_bbdata().
%
% The output data vector is initialized with NaNs.  Data from a bbfile that fit
% within the specified start time and duration are loaded/merged into the data
% vector.
%
% The properties structure includes .t0, the start time vector extended to
% microseconds, .dt, the sample interval, in seconds, .srate, the sample rate,
% in samples per second, and orientation (.vertical, .horizontal).
%
% After an initial call to load_bbdata() to load VFBB data, merge_bbdata()
% merges additional VFBB data.  A warning is given if existing (non-Nan) data
% values differ.  The sample rate and orientation from the properties structure
% must match the values in the bbfile headers.  (VFBB data file headers can be
% viewed using the dhead command-line program.)
%
% The data vector may cross the end of a calendar year, including a leap year.
% However, when comparing the data vector start time (t0) and the VFBB file
% start time, leap seconds are NOT taken into account.  This should not be a
% problem for data files using GPS time, since GPS time is not adjusted for
% leap seconds.  (Ref.: http://tycho.usno.navy.mil/leapsec.html.)
%
% load_bbdata() is useful to extract the initial triggering event or an embedded
% event from within a recording containing multiple events (e.g., the Parkfield
% earthquake recorded at UPSAR).  load_bbdata() combined with merge_bbdata() is
% useful to stitch together multiple recordings into a single time series (e.g.,
% the three records of the start of the San Simeon earthquake recorded at UPSAR,
% or data from the experiment to capture San Andreas tremor events at UPSAR).
%

% 1969 is used for the pivot year of a 2-digit year (XPG5-compatible)

%
% Author:  Lawrence M. Baker
%          U.S. Geological Survey
%          345 Middlefield Road  MS977
%          Menlo Park, CA  94025
%          baker@usgs.gov
%
%                                 Disclaimer
%
% Although  this  program  has  been  used by the U.S. Geological Survey, no
% warranty, expressed or implied, is made by the USGS as to the accuracy and
% functioning  of  the  program  and related program material, nor shall the
% fact of distribution constitute any such warranty, and  no  responsibility
% is assumed by the USGS in connection therewith.
%
%
% Modification History:
%
%  9-Feb-2005  L. M. Baker      Original version.
% 10-Feb-2005  L. M. Baker      Add warnings if orientations don't match.
%  1-Nov-2006  L. M. Baker      Allow the start time vector, t0, to be specified
%                                  to microsecond precision.
%  1-Nov-2006  L. M. Baker      Allow an empty start time vector, t0, ([]) to
%                                  imply the time of the first data sample from
%                                  the VFBB data file integer header.
% 23-Aug-2007  L. M. Baker      Add more comments.
%                               Fix MatLab complaint when duration is not an
%                                  integer (i.e., for a fractional duration).
%                               Help load_bbdata and help merge_bbdata print the
%                                  same text now.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exactly 3 input arguments
if ( nargin ~= 3 )
   error( 'load_bbdata requires 3 input arguments; type help load_bbdata.' );
end

% Exactly 2 output arguments
if ( nargout ~= 2 )
   error( 'load_bbdata returns 2 output results; type help load_bbdata.' );
end

% Initialize default t0 values
year   = 0;
day    = 0;
hour   = 0;
minute = 0;
second = 0;
msec   = 0;
usec   = 0;

% Verify t0 values (if the t0 array is not empty)
ntimes = numel( t0 );
if ( ntimes > 0 )
   if ( ntimes >= 1 )
      year   = t0(1);
   end
   if ( ntimes >= 2 )
      day    = t0(2);
   end
   if ( ntimes >= 3 )
      hour   = t0(3);
   end
   if ( ntimes >= 4 )
      minute = t0(4);
   end
   if ( ntimes >= 5 )
      second = t0(5);
   end
   if ( ntimes >= 6 )
      msec   = t0(6);
   end
   if ( ntimes >= 7 )
      usec   = t0(7);
   end
   if ( year >= 0 & year <= 99 )
      error( 'Invalid t0; use 4-digit year.' );
   end
   if ( ntimes > 7 | ...
        year   < 0 | ...
        day    < 1 | day    > 365 + leap_days( year) | ...
        hour   < 0 | hour   >  23 | ...
        minute < 0 | minute >  59 | ...
        second < 0 | second >  59 | ... % We don't handle leap seconds
        msec   < 0 | msec   > 999 | ...
        usec   < 0 | usec   > 999 )
      error( 'Invalid t0; must be [ year day hour minute second msec usec ].' );
   end
end

% Verify seconds duration
if ( duration <= 0 )
   error(' Invalid duration; must be > 0.' );
end

% Read bbfile integer and real headers
ihdr = bbihdr( bbfile );
rhdr = bbrhdr( bbfile );
inull = ihdr(3);
rnull = rhdr(2);

% Set t0 to the time of the first sample for an empty t0 array ([])
if ( ntimes == 0 )
   year   = ihdr(10);
   if ( year <= 99 )
      if ( year >= 69 )
         year   = 1900 + year;
      else
         year   = 2000 + year;
      end
   end
   day    = ihdr(11);
   hour   = ihdr(12);
   minute = ihdr(13);
   second = ihdr(14);
   msec   = ihdr(15);
   usec   = ihdr(16);
end

vertical   = ihdr(41);
horizontal = ihdr(42);
if ( ( vertical == inull ) | ( horizontal == inull ) )
   orientation = 'unknown';
else
   switch ( mod( vertical, 360 ) )
      case { 0, 180 }
         orientation = 'vertical (positive ';
         if ( vertical == 0 )
            orientation = [ orientation 'up)' ];
         else
            orientation = [ orientation 'down)' ];
         end
      case  90
         orientation = 'horizontal (positive ';
         switch ( mod( horizontal, 360 ) )
            case   0
               orientation = [ orientation 'N)' ];
            case  90
               orientation = [ orientation 'E)' ];
            case 180
               orientation = [ orientation 'S)' ];
            case 270
               orientation = [ orientation 'W)' ];
            otherwise
               orientation = [ orientation int2str( horizontal ) ...
                               ' degrees clockwise from N)' ];
         end
      otherwise
         orientation = [ 'positive ' int2str( vertical ) ...
                         ' degrees down from vertical, ' ...
                         int2str( horizontal ) ...
                         ' degrees clockwise from N' ];
   end
end
motion = ihdr(254);
switch ( motion )
   case 1
      units = 'acceleration (cm/s/s)';
   case 2
      units = 'velocity (cm/s)';
   case 3
      units = 'displacement (cm)';
   case 50
      units = 'volumetric strain';
   otherwise
      units = 'unknown units';
end
samples_per_second = rhdr(5);

disp( [ 'VFBB data is ' orientation ' ' units ', ' ...
        num2str( samples_per_second ) ' samples per second.' ] );

% Allocate the output data vector and fill it with NaNs
no_of_samples = floor( duration * samples_per_second );
data = NaN * ones( 1, no_of_samples );

% Save the start time to microseconds, the sample interval, the sample rate,
% and the orientation
properties.t0 = [ year day hour minute second msec usec ];
properties.dt = 1 / samples_per_second;
properties.srate = samples_per_second;
properties.vertical   = vertical;
properties.horizontal = horizontal;

% Call merge_bbdata() to load the bbfile data
data = merge_bbdata( bbfile, data, properties );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function days = leap_days( year )

days = ( rem( year,   4 ) == 0 ) & ...
       ( rem( year, 100 ) ~= 0 | rem( year, 400 ) == 0 );

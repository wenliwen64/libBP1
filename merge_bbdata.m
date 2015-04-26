function merged_data = merge_bbdata( bbfile, data, properties )
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
% 13-Jun-2006  L. M. Baker      Correct error when offset < 0 (VFBB file
%                                  start time is before data window start time).
%  1-Nov-2006  L. M. Baker      Allow an empty start time vector, t0, ([]) to
%                                  imply the time of the first data sample from
%                                  the bbfile integer header.
% 23-Aug-2007  L. M. Baker      Fix completely wrong determination of sign in
%                                  microseconds_difference()
%                               Help load_bbdata and help merge_bbdata print the
%                                  same text now.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Exactly 3 input arguments
if ( nargin ~= 3 )
   error( 'merge_bbdata requires 3 input arguments; type help merge_bbdata.' );
end

% Exactly 1 output argument
if ( nargout ~= 1 )
   error( 'merge_bbdata returns 1 output result; type help merge_bbdata.' );
end

duration = length( data );
%bbfile
%properties.t0
%properties.dt
%properties.srate
%duration

% Read the bbfile integer and real headers
ihdr   = bbihdr( bbfile );
rhdr   = bbrhdr( bbfile );

% VFBB "null" value markers
inull  = ihdr(3);
% This is a bit of a guess in MatLab.  The VFBB i4null value is most likely
% 0x80008000 (the normal 16-bit inull value, -32768, concatenated to itself).
i4null = ( inull * 65536 ) + mod( inull, 65536 );
rnull  = rhdr(2);

year = ihdr(10);
if ( year <= 99 )
   if ( year >= 69 )
      year = 1900 + year;
   else
      year = 2000 + year;
   end
end
bbfile_time = [ year ihdr(11:16) ];
clock_correction = rhdr(60);
if ( clock_correction ~= rnull )
   microseconds_clock_correction = round( 1000 * 1000 * clock_correction );
   subtract_microseconds( bbfile_time, microseconds_clock_correction );
end

% Verify the sample rate
if ( rhdr(5) ~= properties.srate )
   error( [ bbfile ' has the wrong sample rate: ' ...
            num2str( bbfile_rhdr(5) ) '.' ] );
end

% Verify the orientation
if ( ( ihdr(41) ~= properties.vertical ) | ( ihdr(42) ~= properties.horizontal ) )
   disp( [ 'Warning: ' bbfile ' has the wrong orientation: ' ...
           'vertical=' int2str( ihdr(41) ) ', horizontal=' ...
           int2str( ihdr(42) ) '.' ] );
end

% Determine the VFBB data type (+=Real, -=Integer, ABS()=bytes/sample)
data_type = ihdr(4);
if ( data_type == inull )
   data_type = 0;
end
switch ( data_type )
   case {  0, -2 }
      missing_data = inull;
      block_size   = 256;
    case   -4
      missing_data = i4null;
      block_size   = 128;
    case {  1,  4 }
      missing_data = rnull;
      block_size   = 128;
   otherwise
      error( [ bbfile ' has an unsupported data format.' ] );
end

% Number of samples (data points)
data_blocks = ihdr(31);
last_index  = ihdr(32);
if ( ( data_blocks == inull ) | ( last_index == inull ) )
   np = ihdr(256);
else
   np = block_size * ( data_blocks - 1 ) + last_index;
end

disp( [ bbfile ': ' int2str( np ) ' samples (' num2str( np / rhdr(5) ) ...
        ' s) starting at ' ...
        sprintf( '%4u %3.3u-%2.2u:%2.2u:%2.2u.%3.3u%3.3u.', bbfile_time ) ] );

% Merge the VFBB data into the the output data vector
merged_data = data;

% The offset, in samples, is the offset, in microseconds, times the sample rate,
% in samples per microsecond (samples per second / 1000000 )
microseconds_offset = microseconds_difference( bbfile_time, properties.t0 );
offset = round( microseconds_offset * rhdr(5) / 1000000 );
%offset

if ( ( np + offset <= 0 ) | ( offset >= duration ) )
   disp( '   No VFBB data matches the selected time window.' );
else
   first_pt = max(  1,        1 - offset );
   last_pt  = min( np, duration - offset );
   offset   = max( 0, offset );
   disp( [ '   Merging VFBB data(' int2str( first_pt ) ':' ...
           int2str( last_pt ) ') into data(' ...
           int2str( offset + 1 ) ':' ...
           int2str( offset + last_pt - first_pt + 1 ) ').' ] );
    [ unscaled_data, scale ] = bbdata( bbfile, first_pt, last_pt );
% Serial version:
%    for i = 1 : last_pt - first_pt + 1
%       offset = offset + 1;
%       if ( unscaled_data( i ) ~= missing_data )
%          scaled_data = scale * unscaled_data( i );
%          if ( isnan( merged_data(offset) ) )
%             merged_data(offset) = scaled_data;
%          else
%             if ( merged_data(offset) ~= scaled_data )
%                disp( [ '   Warning: VFBB file data(' int2str( i ) ')=' ...
%                        num2str( scaled_data ) ...
%                        ' does not match data(' int2str( offset ) ')=' ...
%                        num2str( merged_data(offset) ) '.' ] );
%             end
%          end
%       end
%    end
% Vector version:
   scaled_data = scale * unscaled_data;
   missing = unscaled_data == missing_data;
   i = find( missing );
   scaled_data(i) = NaN;
   i = find( isnan( merged_data(offset+1:offset+last_pt-first_pt+1) ) );
   merged_data(offset+i) = scaled_data(i);
   i = find( ~missing );
   if ( merged_data(offset+i) ~= scaled_data(i) )
      disp( '   Warning: some VFBB data does not match the existing data.' );
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function add_microseconds( time, microseconds )

if ( microseconds < 0 )
   subtract_microseconds( time, -microseconds );
else
   time(7) = time(7) + microseconds;
%   Adjust the time for any carries required
   if ( time(7) > 999 )
      while ( time(7) > 999 )
         time(7) = time(7) - 1000;
         time(6) = time(6) +    1;
      end
      while ( time(6) > 999 )
         time(6) = time(6) - 1000;
         time(5) = time(5) +    1;
      end
      while ( time(5) > 59 )
         time(5) = time(5) - 60;
         time(4) = time(4) +  1;
      end
      while ( time(4) > 59 )
         time(4) = time(4) - 60;
         time(3) = time(3) +  1;
      end
      while ( time(3) > 23 )
         time(3) = time(3) - 24;
         time(2) = time(2) +  1;
      end
      while ( time(2) > 365 + leap_days( time(1) ) )
         time(2) = time(2) - ( 365 + leap_days( time(1) ) );
         time(1) = time(1) + 1;
      end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subtract_microseconds( time, microseconds )

if ( microseconds < 0 )
   add_microseconds( time, -microseconds );
else
   time(7) = time(7) - microseconds;
%   Adjust the time for any borrows required
   if ( time(7) < 0 )
      while ( time(7) < 0 )
         time(6) = time(6) -    1;
         time(7) = time(7) + 1000;
      end
      while ( time(6) < 0 )
         time(5) = time(5) -    1;
         time(6) = time(6) + 1000;
      end
      while ( time(5) < 0 )
         time(4) = time(4) -  1;
         time(5) = time(5) + 60;
      end
      while ( time(4) < 0 )
         time(3) = time(3) -  1;
         time(4) = time(4) + 60;
      end
      while ( time(3) < 0 )
         time(2) = time(2) -  1;
         time(3) = time(3) + 24;
      end
      while ( time(2) < 1 )
         time(1) = time(1) -   1;
         time(2) = time(2) + 365 + leap_days( time(1) );
      end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function microseconds = microseconds_difference( time1, time2 )
% microseconds_difference = time1 - time2, in microseconds

% Save the sign and compute the absolute value of the time difference
diff = 0;
for i = 1:7
   diff = time1(i) - time2(i);
   if ( diff ~= 0 )
      break;
   end
end
if ( diff >= 0 )
   sign  = +1;
   delta = time1 - time2;
   year  = time1(1) - 1;
else
   sign  = -1;
   delta = time2 - time1;
   year  = time2(1) - 1;
end

% Adjust the time difference for any borrows required
while ( delta(7) < 0 )
   delta(6) = delta(6) -    1;
   delta(7) = delta(7) + 1000;
end
while ( delta(6) < 0 )
   delta(5) = delta(5) -    1;
   delta(6) = delta(6) + 1000;
end
while ( delta(5) < 0 )
   delta(4) = delta(4) -  1;
   delta(5) = delta(5) + 60;
end
while ( delta(4) < 0 )
   delta(3) = delta(3) -  1;
   delta(4) = delta(4) + 60;
end
while ( delta(3) < 0 )
   delta(2) = delta(2) -  1;
   delta(3) = delta(3) + 24;
end
while ( delta(2) < 0 )
   delta(1) = delta(1) -   1;
   delta(2) = delta(2) + 365 + leap_days( year );
   year = year - 1;
end

% Compute the time difference in microseconds
days = delta(2);
while ( delta(1) > 0 )
   delta(1) = delta(1) - 1;
   days     = days + 365 + leap_days( year );
   year     = year - 1;
end
hours        = ( days         *   24 ) + delta(3);
minutes      = ( hours        *   60 ) + delta(4);
seconds      = ( minutes      *   60 ) + delta(5);
milliseconds = ( seconds      * 1000 ) + delta(6);
microseconds = ( milliseconds * 1000 ) + delta(7);

% Apply the correct sign to the difference
microseconds = sign * microseconds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function days = leap_days( year )

days = ( rem( year,   4 ) == 0 ) & ...
       ( rem( year, 100 ) ~= 0 | rem( year, 400 ) == 0 );
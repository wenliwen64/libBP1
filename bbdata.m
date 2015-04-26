% bbdata - MATLAB .mex function to read a VFBB file's time-series.
%
% [ unscaled_data, scale ] = bbdata( 'file' , first_pt , last_pt );
%     scaled_data          = bbdata( 'file' , first_pt , last_pt );
% [ unscaled_data, scale ] = bbdata( 'file' );
%     scaled_data          = bbdata( 'file' );
%
% scaled_data   = scaled time-series data (engineering units, e.g., cm/s)
% unscaled_data = unscaled time-series data (counts)
% scale         = scale factor to convert unscaled_data to scaled_data
% 'file'        = VFBB-format file name
% first_pt      = index of first data point (numbered from 1)
% last_pt       = index of last data point (0=end-of-file)
%
% bbdata returns the time-series from a VFBB-format (blocked-binary) file.
%
% To request a portion of the time-series, specify first_pt and last_pt.
% last_pt is always minimized with the number of data points in the file
% (without error).
%
% See also bbihdr, bbrhdr, bbthdr.

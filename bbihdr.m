% bbihdr - MATLAB .mex function to read a VFBB file's integer header(s).
%
% ihdr = bbihdr( 'file' );
% ihdr = bbihdr( 'file', ihdr_no );
%
% ihdr    = contents of integer header number ihdr_no
% ihdr_no = integer header number (default=1)
% 'file'  = VFBB-format file name
%
% bbihdr returns the contents of the requested integer header from a VFBB-format
% (blocked-binary) file.
%
% See also bbrhdr, bbthdr, bbdata.

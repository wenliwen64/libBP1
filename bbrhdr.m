% bbrhdr - MATLAB .mex function to read a VFBB file's real header(s).
%
% rhdr = bbrhdr( 'file' );
% rhdr = bbrhdr( 'file', rhdr_no );
%
% rhdr    = contents of real header number rhdr_no
% rhdr_no = real header number (default=1)
% 'file'  = VFBB-format file name
%
% bbrhdr returns the contents of the requested real header from a VFBB-format
% (blocked-binary) file.
%
% See also bbihdr, bbthdr, bbdata.

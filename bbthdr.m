% bbthdr - MATLAB .mex function to read a VFBB file's text header(s).
%
% thdr = bbthdr( 'file' );
% thdr = bbthdr( 'file', thdr_no );
%
% thdr    = contents of text header number thdr_no
% thdr_no = text header number (default=1)
% 'file'  = VFBB-format file name
%
% bbthdr returns the contents of the requested text header from a VFBB-format
% (blocked-binary) file as a MatLab string.
%
% Note:  Text headers are optional in VFBB-format files, and very few use them.
%
% See also bbihdr, bbrhdr, bbdata.

function of=osdep
% Returns the value of the read option to be used
% on the local operating system for files created
% on the SOLARIS operating system.

%if strcmp(getenv('OSTYPE'),'linux')
if strmatch('linux',getenv('OSTYPE'))
  of= 'b'; 
end
if strcmp(getenv('OSTYPE'),'solaris')
    of= 'l'; 
end

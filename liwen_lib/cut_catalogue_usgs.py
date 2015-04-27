ifile = open('/home/lwen/Documents/seis_research/nepal_2015/catalogue_nepal_2015_8hrs.usgs.raw.dat', 'r');
ofile = open('/home/lwen/Documents/seis_research/nepal_2015/catalogue_nepal_2015_8hrs.usgs.dat', 'w');
for line in ifile:
    substrings = line.split(',');

    time = substrings[0]; 
    lat = float(substrings[1]);
    lon = float(substrings[2]);
    dep = float(substrings[3]);
    mag = float(substrings[4]);

    subtime = time.split('T');
    date_str = subtime[0];
    time_str = subtime[1];
    year, mon, day = date_str.split('-');
    hour, minu, sec = time_str.split(':');
    sec = sec[:-1];

    newline = ''.join([year, ' ', mon, ' ', day, ' ', hour, ' ', minu, ' ', sec, ' ', str(lat), ' ', str(lon), ' ', str(mag), '\n']);
    ofile.write(newline);
    print mag 

// This script is used to parse the catalog file provided by Lingsen Meng which is from Caltech network.
void cut_catalog(){
    ifstream file;
    ofstream ofile("./eqinfo_processed.dat");
    file.open("/Users/lwen/Downloads/eqinfo_2014.nchile.origin.dat");//Input file
    char line[500];
    if(file.is_open()) cout<<"open"<<endl;
    cout<<"happy"<<endl;
    char part1[100], part2[100], part3[100], part4[100], lat[100], lon[100], dummy1[100], mag[100], mag_str[100];
    while(file>>part1>>part2>>part3>>part4>>lat>>lon>>dummy1>>mag>>mag_str){
        int cst_day, cst_mon, cst_year, cst_hour, cst_min, cst_sec,  utc_day, utc_mon, utc_year, utc_hour, utc_min, utc_sec;
        
        sscanf(part3, "%d/%d/%d", &utc_day, &utc_mon, &utc_year);
        sscanf(part4, "%d:%d:%d", &utc_hour, &utc_min, &utc_sec);
        ofile << utc_year << " " << utc_mon << " " << utc_day << " " << utc_hour << " " << utc_min << " " << utc_sec << " " << lat << " " << lon << " " << mag << std::endl;
    }
}

figure(2);
% name=[71042612 elt31 78123001 76101000 72081900 epll01wt dpsn01wt     dsdp08gc kk831116 87001611 v2113      erdc11wt NBP0207 v2403 aria01wt]
name=['aria01wt ';'erdc11wt ';'v2403    ';'dsdp08gc ';'kk831116a';'kk831116b';'87001611 ';'kk831116c';'kk831116d';'v2113    ';'dpsn01wta';'epll01wta';'epll01wtb';'dpsn01wtb';'71042612 ';'elt31    ';'78123001 ';]
lat =[1.4          1.68        2.8            6.1          8.0        9.0        9.95         10.07        12.0       12.8        18.5       18.2      19.5        19.3          19.5          20.5        21]
skw= [224        229         218           180           178        179        180         175          173        170          63         69        79         62            95          120           94]
L=   [220        219         222           175           173        174        173          165          156        160         61          66       77         55            89           116          75]
H=   [228        234         224           187           190        185        183          183          197        173         66          72       86         76            98           124          101]
namep=['dsdp08gc';'87001611';'elt31   ';'v2113   ';'dpsn01wt';'epll01wt';'71042612']
latp  =[6.1         9.95      12.4        12.8       17         17.5       19]
skwp  =[202         200       134         172        95         160        132 ]
%figure(1);
hold on;
errorbar(lat, skw,skw-L,H-skw,'b.','MarkerSize',5);
%plot(lat(12:17), skw(12:17),'r*');
for kk=1:17
temp=name(kk,:)
%text(lat(kk),L(kk)+5,temp,'FontSize',5)
end
 %ylim([0 270]);

% plot(latp, skwp,'b-*');

for kk=1:7
temp=namep(kk,:)
%text(latp(kk),skwp(kk)+5,temp)
end
title('latitude vs skewness (for chron 30)');
ylabel('skewness (degree)');
xlabel('latitude (degree)');
% for i=1:16
% ddl(i)=(skw(i+1)-skw(i))/(lat(i+1)-lat(i));
% end
%plot(lat(1:16),ddl,'-*');
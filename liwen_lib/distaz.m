function [dk,dd,daze,dazs] = distaz(sta,sto,epa,epo)
% This function is provided by specfem3d_glob, 
% copied from https://github.com/geodynamics/specfem3d_globe/blob/master/utils/Visualization/VTK_ParaView/matlab/distaz.m
% input: 
%     sta, sto = staion latitude, longitude
%     epa, epo = epicenter latitude, longitude
% output:
%     dk = distance in kms
%     dd = distance in degree
%     daze = back azimuthal angle (from north to epicenter@station)
%     dazs = azimuthal angle (from north to station@epicenter)
% Example:
%     [dk,dd,daze,dazs] = distaz(0, 0, 30, 0)
        rad=pi/180.0d0;

        sa  = atan(.993270*tan(sta*rad));
        ea  = atan(.993270*tan(epa*rad));
        ssa = sin(sa);
        csa = cos(sa);
        so  = sto*rad;
        eo  = epo*rad;
        sea = sin(ea);
        cea = cos(ea);
        ces = cos(eo-so);
        ses = sin(eo-so);
        
        if  (sa==ea)
           if (sto==epo)
              	dk =0.00;
					dd =0.00;
	  				daze=0.0;
               dazs=0.0;
                 return
  end
end

if sta==90.
   if epa==90.0
     dk =0.00;
	  dd =0.00;
	  daze=0.00;
     dazs=0.00;
     return
  end
end

      
      
if sta==-90.0
   if epa==-90.0
	  dk =0.00;
	  dd =0.00;
	  daze=0.00;
	  dazs=0.00;
	  return
	end
end
        dd = ssa*sea+csa*cea*ces;
        if dd~=0. , dd=atan(sqrt(1.0-dd*dd)/dd); end
        if dd==0. , dd=pi/2.0; end
        if dd<0.0 , dd=dd+pi; end
        dd = dd/rad;
        dk = dd*111.19;

        dazs = atan2(-ses,(ssa/csa*cea-sea*ces));
        daze = atan2(ses,(sea*csa/cea-ssa*ces));
        dazs = dazs/rad;
        daze = daze/rad;
        if (dazs<0.00), dazs=dazs+360.0; end
        if (daze<0.00), daze=daze+360.0; end


return
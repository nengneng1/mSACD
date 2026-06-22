function [raw_m, raw_phi] = polar_convert(phasor_struct)

z_struct=phasor_struct;

for pp=1:numel(z_struct.G)
if (z_struct.G(pp)>0 && z_struct.S(pp)>0)
raw_phi(pp)=atan(z_struct.S(pp)./z_struct.G(pp));
raw_m(pp)=sqrt(((z_struct.G(pp)).^2)+((z_struct.S(pp)).^2));

elseif (z_struct.G(pp)<0 && z_struct.S(pp)>0)
raw_phi(pp)=pi+atan(z_struct.S(pp)./z_struct.G(pp));
raw_m(pp)=sqrt(((z_struct.G(pp)).^2)+((z_struct.S(pp)).^2));

elseif (z_struct.G(pp)>0 && z_struct.S(pp)<0)
raw_phi(pp)=atan(z_struct.S(pp)./z_struct.G(pp));
raw_m(pp)=sqrt(((z_struct.G(pp)).^2)+((z_struct.S(pp)).^2));

elseif (z_struct.G(pp)<0 && z_struct.S(pp)<0)
raw_phi(pp)=pi+atan(z_struct.S(pp)./z_struct.G(pp));
raw_m(pp)=sqrt(((z_struct.G(pp)).^2)+((z_struct.S(pp)).^2));

elseif (z_struct.G(pp)==0 && z_struct.S(pp)>0)
raw_phi(pp)=pi/2;
raw_m(pp)=sqrt(((z_struct.G(pp)).^2)+((z_struct.S(pp)).^2));

elseif (z_struct.G(pp)==0 && z_struct.S(pp)<0)
raw_phi(pp)=3*pi/2;
raw_m(pp)=sqrt(((z_struct.G(pp)).^2)+((z_struct.S(pp)).^2));

elseif (z_struct.G(pp)>0 && z_struct.S(pp)==0)
raw_phi(pp)=0;
raw_m(pp)=sqrt(((z_struct.G(pp)).^2)+((z_struct.S(pp)).^2));

else
raw_phi(pp)=pi; 
raw_m(pp)=sqrt(((z_struct.G(pp)).^2)+((z_struct.S(pp)).^2));

end    
end




end
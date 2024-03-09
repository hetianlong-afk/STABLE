function [V_load,V_hom_load_0]=VoltageLoadCalc_matlab(V_hom_load_0,V_load_cpu,TbAng_coef_hom,pattern)
j=1;
n=length(pattern);
V_load=zeros(1,n);
for i =1:n
    if pattern(i)==1
        V_load(i) = V_hom_load_0*TbAng_coef_hom;
        V_hom_load_0 = V_load(i)+V_load_cpu(j);
        j=j+1;
    else
        V_hom_load_0 = V_hom_load_0*TbAng_coef_hom;
        V_load(i) = V_hom_load_0;
    end
end
end
clear all;

SF_STD = 6.0;
gate_count = 68215;
sram.ps_a = -614.9807;       
sram.ps_b = 3528.4329;

vdd_range = 0.3:0.05:1.1;
temp_range = 0:2:110;
vdd_index = 0;
P_s_x = [];
P_s_y = [];
P_s_z = [];
for vdd = vdd_range
    vdd_index = vdd_index + 1;
    tmpr_index = 0;
    for tmpr = temp_range
        tmpr_index = tmpr_index + 1;
        P_s(vdd_index, tmpr_index) = SF_STD * gate_count * vdd * ...
            (2.48325558 + 2.26321378 + 2.6529006) / 3.0 * 1e-8 * ...
            (tmpr + 273) * (tmpr + 273) * exp(-(sram.ps_a * vdd + ...
            sram.ps_b) / (tmpr + 273));
        P_s_1 = SF_STD * gate_count * vdd * ...
            (2.48325558 + 2.26321378 + 2.6529006) / 3.0 * 1e-8 * ...
            (tmpr + 273) * (tmpr + 273) * exp(-(sram.ps_a * vdd + ...
            sram.ps_b) / (tmpr + 273));
        
        P_s_x = [P_s_x vdd];
        P_s_y = [P_s_y tmpr];
        P_s_z = [P_s_z P_s_1];
        
        x_s(vdd_index, tmpr_index) = vdd;
        y_s(vdd_index, tmpr_index) = tmpr;
    end
end

tri = delaunay(P_s_x, P_s_y,{'Qt','Qbb','Qc','Qz'});
trisurf(tri, P_s_x, P_s_y, P_s_z*60)

fid = fopen('leakage_surf_values.dat','w');
for i = 1:size(tri,1)
    fprintf(fid,'%f %f %f\n', P_s_x(tri(i,1)), P_s_y(tri(i,1)), ...
        P_s_z(tri(i,1))*60);
    fprintf(fid,'%f %f %f\n', P_s_x(tri(i,2)), P_s_y(tri(i,2)), ...
        P_s_z(tri(i,2))*60);
    fprintf(fid,'%f %f %f\n\n', P_s_x(tri(i,3)), P_s_y(tri(i,3)), ...
        P_s_z(tri(i,3))*60);
end
fclose(fid);

P_s;

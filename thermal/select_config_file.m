function [hs_config_file, s_spreader, s_sink] = select_config_file(mpsoc_flp, nx, ny)
%Selects the configuration file for simulation
%The dimensions of 
%  Die = TIM < Speader< Heat Sink 
%         # ev6 die side in meters
%         -s_die                0.016
%         # ev6 tim side in meters
%         -s_tim                0.016
%		 # package substrate side in meters
%		  -s_sub                0.021
%		 #	solder ball side in meters
%         -s_solder             0.021
%         # spreader side in meters
% 		  -s_spreader			0.03
% 		  # heatsink side in meters
% 		  -s_sink				0.06
%		  # PCB side in meters
%		  -s_pcb                0.1

% keep the area ratio of the die:tim:spreader:heat_sink same for all
% architectures/ floorplans. The floor plan is assumed to be rectangular

% eg usage:
% [a b c]=select_config_file('mpsoc_ev6_hetero_2x2_std.flp', 2,2)

%mpsoc_flp='mpsoc_ev6_hetero_4x4.flp'

switch (nargin)
    case 1
        % if only one argument
        if isa(mpsoc_flp, 'char')
            [~, mpsoc_flp]= flp2struct(mpsoc_flp)
        elseif ~isa(mpsoc_flp, 'struct')
            error('Incorrect floorplan file')
        end
        
        ev6_die_area=0.016^2;
        
        die_area= mpsoc_flp.Width * mpsoc_flp.Height;
        
        tim_area= die_area;
        extra_hspdr_area=0.03^2-ev6_die_area;
        extra_hs_area =0.06^2 - ev6_die_area;
        extra_pck_sub_area=0.021^2-ev6_die_area;
        extra_pcb_side_area=0.1^2- ev6_die_area;
        
        spreader_area= die_area+extra_hspdr_area;
        heatsink_area= die_area+extra_hs_area
        pck_sub_area = die_area+ extra_pck_sub_area;
        pcb_side_area= die_area+ extra_pcb_side_area;
        
        spreader_area_ideal=0.03^2*die_area/ev6_die_area;
        heatsink_area_ideal=0.06^2*die_area/ev6_die_area;
        
        s_spreader=sqrt(spreader_area)
        s_sink=sqrt(heatsink_area)
        s_sub =sqrt(pck_sub_area)
        s_solder= s_sub
        s_pcb = sqrt(pcb_side_area)
        
     
        
        isfield(mpsoc_flp, 'TotalCores')
        mpsoc_flp.NoOfBlocks
        TotalCores=mpsoc_flp.NoOfBlocks/18
        
        
        switch (TotalCores)
            case 1
                hs_config_file ='hotspot.config';
            case 4
                hs_config_file ='hotspot_2x2.config';
            case 9
                hs_config_file ='hotspot_3x3.config';
            case 16
                hs_config_file ='hotspot_4x4.config';
            case 25
                hs_config_file ='hotspot_5x5.config';
            otherwise
                %dynamically generate the config file
                
                
        end
        
        
    case 2
        error('Incorrect arguments')
        
    case 3
        % if the first argument is char or struct 
        if isa(mpsoc_flp, 'char')
            [~, mpsoc_flp]=flp2struct(mpsoc_flp);
        elseif ~isa(mpsoc_flp, 'struct')
            error('Incorrect floorplan file')
        end
        
        ev6_die_area=0.016^2;
        die_area= mpsoc_flp.Width * mpsoc_flp.Height;
        
        tim_area= die_area;
        tim_area= die_area;
        extra_hspdr_area=0.03^2-ev6_die_area;
        extra_hs_area =0.06^2 - ev6_die_area;
        extra_pck_sub_area=0.021^2-ev6_die_area;
        extra_pcb_side_area=0.1^2- ev6_die_area;
        
        spreader_area= die_area+extra_hspdr_area;
        heatsink_area= die_area+extra_hs_area
        pck_sub_area = die_area+ extra_pck_sub_area;
        pcb_side_area= die_area+ extra_pcb_side_area;
        
        spreader_area_ideal=0.03^2*die_area/ev6_die_area;
        heatsink_area_ideal=0.06^2*die_area/ev6_die_area;
        
        s_spreader=sqrt(spreader_area)
        s_sink=sqrt(heatsink_area)
        s_sub =sqrt(pck_sub_area)
        s_solder= s_sub
        s_pcb = sqrt(pcb_side_area)
        
        
        if (nx==1)&& (ny==1)
            hs_config_file ='hotspot.config';
            
        elseif (nx==2) && (ny==2)
            hs_config_file ='hotspot_2x2.config';
            
        elseif (nx==3) &&( ny==3)
            hs_config_file ='hotspot_3x3.config';
            
        elseif (nx==4) &&( ny==4)
            hs_config_file ='hotspot_4x4.config';
        elseif (nx==5) &&( ny==5)
            hs_config_file ='hotspot_4x4.config';
              
            
        else
            %dynamically generate the config file
        end
end


end
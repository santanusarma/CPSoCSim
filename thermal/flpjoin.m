function  joinedcore_flp= flpjoin( core1_flp, core2_flp, position)

% core1_flp = floorplan structure of the core 1, core1 is the reference
% core2_flp = floorplan structure of the core 2
% position = is a string showing 'left', 'right', 'up', 'down',
% 'lowerleft', 'upperleft', 'lowerright', 'upperright' 
% 

% TO DO:
% Naming conventions and labels name to b eupdated.



switch position
    
    case 'left'
        %get the base cordinates of core 1 and core 2 and shift the base
        %corrdinates of core 2 to left by its width to abut core 1 
        BaseX=core1_flp.BaseX-core2_flp.BaseX-core2_flp.Width;
        BaseY=core1_flp.BaseY;
        
        [core2_cell, core2_flp1]=flp_oper(core2_flp, 'shift', [BaseX, BaseY])
        
        %rename the lables
        for i=1:length(core2_flp1.Labels)
        core2_flp1.Labels{i}=[core2_flp1.Labels{i}, '_L'];
        end
        for j=1:length(core1_flp.Labels)
        core1_flp.Labels{j} = [core1_flp.Labels{j}, '_R'];
        end
        
        joinedcore_flp.Labels={core1_flp.Labels{:}, core2_flp1.Labels{:}}
        
        joinedcore_flp.BaseX=core2_flp1.BaseX;
        joinedcore_flp.BaseY=core2_flp1.BaseY;
        
        joinedcore_flp.w =[core1_flp.w; core2_flp1.w]
        joinedcore_flp.h =[core1_flp.h; core2_flp1.h]
        joinedcore_flp.x =[core1_flp.x; core2_flp1.x]
        joinedcore_flp.y =[core1_flp.y; core2_flp1.y]
        
        joinedcore_flp.EmptyLabels={core1_flp.EmptyLabels{:}, ...
                                    core2_flp1.EmptyLabels{:}}
        
        joinedcore_flp.NoOfBlocks=core1_flp.NoOfBlocks+core2_flp1.NoOfBlocks
        joinedcore_flp.Width = core1_flp.Width+ core2_flp1.Width
        joinedcore_flp.Height= core1_flp.Height
        
    case 'right'
        %get the base cordinates of core 1 and core 2 and shift the base
        %corrdinates of core 2 to right by its width to abut core 1 
        BaseX=core1_flp.BaseX+core2_flp.BaseX+core2_flp.Width;
        BaseY=core1_flp.BaseY;
        
        [core2_cell, core2_flp1]=flp_oper(core2_flp, 'shift', [BaseX, BaseY])
        
        %rename the lables
        for i=1:length(core2_flp1.Labels)
        core2_flp1.Labels{i}=[core2_flp1.Labels{i}, '_R'];
        end
        for j=1:length(core1_flp.Labels)
        core1_flp.Labels{j} = [core1_flp.Labels{j}, '_L'];
        end

        
        joinedcore_flp.Labels={core1_flp.Labels{:}, core2_flp1.Labels{:}}
        
        joinedcore_flp.BaseX=core2_flp1.BaseX;
        joinedcore_flp.BaseY=core2_flp1.BaseY;
        
        joinedcore_flp.w =[core1_flp.w; core2_flp1.w]
        joinedcore_flp.h =[core1_flp.h; core2_flp1.h]
        joinedcore_flp.x =[core1_flp.x; core2_flp1.x]
        joinedcore_flp.y =[core1_flp.y; core2_flp1.y]
        
        joinedcore_flp.EmptyLabels={core1_flp.EmptyLabels{:}, ...
                                    core2_flp1.EmptyLabels{:}}
        
        joinedcore_flp.NoOfBlocks=core1_flp.NoOfBlocks+core2_flp1.NoOfBlocks
        joinedcore_flp.Width = core1_flp.Width+ core2_flp1.Width
        joinedcore_flp.Height= core1_flp.Height
        
        
        
    case 'up'
        %get the base cordinates of core 1 and core 2 and shift the base
        %corrdinates of core 2 to upside by its width to abut core 1 
        BaseX=core1_flp.BaseX;
        BaseY=core1_flp.BaseY+core2_flp.BaseY+core2_flp.Height;
        
        [core2_cell, core2_flp1]=flp_oper(core2_flp, 'shift', [BaseX, BaseY])
        
        %rename the lables
        for i=1:length(core2_flp1.Labels)
        core2_flp1.Labels{i}=[core2_flp1.Labels{i}, '_Up'];
        end
        for j=1:length(core1_flp.Labels)
        core1_flp.Labels{j} = [core1_flp.Labels{j}, '_Dn'];
        end

        
        
        joinedcore_flp.BaseX=core2_flp1.BaseX;
        joinedcore_flp.BaseY=core2_flp1.BaseY;

        
        joinedcore_flp.Labels={core1_flp.Labels{:}, core2_flp1.Labels{:}}
        
        joinedcore_flp.w =[core1_flp.w; core2_flp1.w]
        joinedcore_flp.h =[core1_flp.h; core2_flp1.h]
        joinedcore_flp.x =[core1_flp.x; core2_flp1.x]
        joinedcore_flp.y =[core1_flp.y; core2_flp1.y]
        
        joinedcore_flp.EmptyLabels={core1_flp.EmptyLabels{:}, ...
                                    core2_flp1.EmptyLabels{:}}
        
        joinedcore_flp.NoOfBlocks=core1_flp.NoOfBlocks+core2_flp1.NoOfBlocks
        joinedcore_flp.Width = core1_flp.Width
        joinedcore_flp.Height= core1_flp.Height+core2_flp.Height
       
        
    case 'down'
       %get the base cordinates of core 1 and core 2 and shift the base
        %corrdinates of core 2 to down by its width to abut core 1 
        BaseX=core1_flp.BaseX;
        BaseY=core1_flp.BaseY-core2_flp.BaseY-core2_flp.Height;
        
        [core2_cell, core2_flp1]=flp_oper(core2_flp, 'shift', [BaseX, BaseY])
        
        joinedcore_flp.Labels={core1_flp.Labels{:}, core2_flp1.Labels{:}}
        
        %rename the lables
        for i=1:length(core2_flp1.Labels)
        core2_flp1.Labels{i}=[core2_flp1.Labels{i}, '_Dn'];
        end
        for j=1:length(core1_flp.Labels)
        core1_flp.Labels{j} = [core1_flp.Labels{j}, '_Up'];
        end

        
        joinedcore_flp.BaseX=core2_flp1.BaseX;
        joinedcore_flp.BaseY=core2_flp1.BaseY;
        
        joinedcore_flp.w =[core1_flp.w; core2_flp1.w]
        joinedcore_flp.h =[core1_flp.h; core2_flp1.h]
        joinedcore_flp.x =[core1_flp.x; core2_flp1.x]
        joinedcore_flp.y =[core1_flp.y; core2_flp1.y]
        
        joinedcore_flp.EmptyLabels={core1_flp.EmptyLabels{:}, ...
                                    core2_flp1.EmptyLabels{:}}
        
        joinedcore_flp.NoOfBlocks=core1_flp.NoOfBlocks+core2_flp1.NoOfBlocks
        joinedcore_flp.Width = core1_flp.Width
        joinedcore_flp.Height= core1_flp.Height+core2_flp.Height
        
        
    case 'lowerleft'
        %tbd
        
    case 'upperleft'
        
    case 'lowerright'
        
    case 'upperright'
        
end
        


end
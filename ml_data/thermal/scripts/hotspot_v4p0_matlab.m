
%8x8 core simulation (check teh mode)
!./hotspot -c hotspot_2x2.config  -f mpsoc_8x8.flp -p mpsoc_8x8_749.ptrace -o mpsoc_8x8_749.ttrace -steady_file mpsoc_8x8_749.steady >>mpsoc_8x8_749_log.txt
!./hotspot -c hotspot_2x2.config  -f mpsoc_8x8.flp -p mpsoc_8x8_981_0_0625.ptrace -o mpsoc_8x8_981_0_0625.ttrace -steady_file mpsoc_8x8_981_0_0625.steady >>mpsoc_8x8_981_0_0625.ptrace_log.txt

%2x2 core simulation (check the model)
nx=2
ny=2
NoOfBenchmarks=25
NoOfTraces =3000
[mpsoc_flp mpsoc_flp_file]=create_mpsoc_flp('ev6.flp',nx,ny)
mpsoc_config_file=select_config_file(mpsoc_flp_file,nx,ny)

!./hotspot -c hotspot_2x2.config  -f mpsoc_ev6_2x2.flpp -p gcc_2x2.ptrace -o gcc_2x2.ttrace -steady_file gcc_2x2.steady >>gcc_2x2_log.txt



%hotspot execution from matlab in windows

!./hotspot -c hotspot.config -init_file gcc.init -f ev6.flp -p gcc.ptrace -o gcc.ttrace -steady_file gcc.steady

%!./hotspot -c hotspot.config  -f ev6.flp -p gcc.ptrace -o gcc.ttrace -steady_file gcc.steady

% Save the temperature tarce to  gcc1.ttrace
!./hotspot -c hotspot.config -init_file gcc.init -f ev6.flp -p gcc.ptrace -o gcc1.ttrace


%Run hotspot for ev6_1x1.flp for the 20 block floor plan



%Run hotspot with grid model ( output in BLOCK level):
!./hotspot -c hotspot.config -f ev6.flp -p gcc.ptrace -steady_file gcc.steady -model_type grid

%Run hotspot with grid model for specified grid dimensions (output in BLOCK level):
!./hotspot -c hotspot.config -init_file gcc.init -f ev6.flp -p gcc.ptrace -o gcc_grid.ttrace -model_type grid -grid_rows 64 -grid_cols 64

%Run hotspot with grid model for DEFAULT grid dimensions (output in GRID level):
%The grid model produces 50x50 gird by default

!./hotspot -c hotspot.config -f ev6.flp -p gcc.ptrace -steady_file gcc.steady -model_type grid -grid_steady_file gcc.grid.steady
!./grid_thermal_map.pl ev6.flp gcc.grid.steady > gcc2_v4p0.svg

%Run hotspot with grid model for SPECIFIED grid dimensions (output in GRID level):
!./hotspot -c hotspot.config -f ev6.flp -p gcc.ptrace -steady_file gcc.steady -model_type grid -grid_steady_file gcc.grid64x64.steady -grid_rows 64 -grid_cols 64

!./grid_thermal_map.pl ev6.flp gcc.grid64x64.steady > gcc2_64x64.svg

!./hotspot -c hotspot.config -f ev6.flp -p gcc.ptrace -steady_file gcc.steady -model_type grid -grid_steady_file gcc.grid55x55.steady -grid_rows 55 -grid_cols 55
!./grid_thermal_map.pl ev6.flp gcc.grid55x55.steady > gcc2_55x55.svg
!./grid_thermal_map_v5p0.pl ev6.flp gcc.grid55x55.steady > gcc2_55x55.svg


%Simulate for trcae file in sperate director by providing the path 

!./hotspot -c hotspot.config -f ev6.flp -p ../power_data/ptrace/gcc.ptrace -steady_file gcc.steady -model_type grid -grid_steady_file gcc.grid55x55.steady -grid_rows 55 -grid_cols 55

%The conver function is not working
% download imagemagic  http://www.imagemagick.org/script/binary-releases.php#macosx
% download imagemagic for mac http://cactuslab.com/imagemagick/
% Use the InkScape as the default SVG viewer

% The convert function is not working in the imagemagic toolbox 
!./ImageMagick-6.8.7/bin/convert -font Helvetica svg:gcc.svg gcc2.pdf 

%Use the manual conversion in InkScape toolbox or the commandline options


% The user can select the mapping between the grid
% and block temperatures of the grid model through the command line
% option '-grid_map_mode <mode>'. The four mapping modes supported are:
% 'min', 'max', 'avg' and 'center'. The first three options respectively
% mean that the block's temperature is computed as the minimum, maximum,
% or average temperature of the grid cells in it. The 'center' option
% means that the block's temperature is given by the temperature of the
% grid cell at its center.

!./hotspot -c hotspot.config -f ev6.flp -p gcc.ptrace -steady_file gcc.steady -model_type grid -grid_steady_file gcc.grid.steady -grid_map_mode max 

!./hotspot -c hotspot.config -f ev6.flp -p gcc.ptrace -steady_file gcc.steady -model_type grid -grid_steady_file gcc.grid.steady -grid_map_mode min

!./hotspot -c hotspot.config -f ev6.flp -p gcc.ptrace -steady_file gcc.steady -model_type grid -grid_steady_file gcc.grid.steady -grid_map_mode avg

!./hotspot -c hotspot.config -f ev6.flp -p gcc.ptrace -steady_file gcc.steady -model_type grid -grid_steady_file gcc.grid.steady -grid_map_mode center




%Load the grid data and display the data
gcc_grid=load('gcc.grid.steady')


%Grid based simulation and display of data
x=vect2mat('gcc.grid.steady', 50, 50)
view_grid_ttrace('gcc.grid55x55.steady', 55, 55, 'waterfall')
view_grid_ttrace('gcc.grid55x55.steady', 55, 55, 'mesh')
view_grid_ttrace('gcc.grid55x55.steady', 55, 55, 'mesh', 'Centigrade')
view_grid_ttrace('gcc.grid55x55.steady', 55, 55, 'surf', 'Centigrade')
view_grid_ttrace('gcc.grid55x55.steady', 55, 55, 'waterfall', 'Centigrade')

%==========================================================================
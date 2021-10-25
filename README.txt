This is a plain text showing the usage of the package.

The core looping program consists of efficiency.h, efficiency.C, makehists.h, makehists.C and Control.C

To compile, you can simply use 
. gcc.sh

After that, you can run 
./efficiency.exe filelist category filenumber

After looping each file, scaling needs tp be applies to normalize the genWeight. Please go to the scaling+eff repository to do the scaling and generating efficiency plots.

Plotting programs are seperately written in plot_program directory, there is also a gcc_plot.sh to compile. Stack plots anc comparison plots are available.

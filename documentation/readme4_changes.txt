FIles that have been modified by F.Lavorenti during 1st year of PhD (01/2021 - 09/2021):

--> CMakeLists.txt: (i) added twice a library on irene which was difficult to find during th compilation - now commented

--> fields/EMfields3D: (i) added divergence cleaning for B at boundaries, (ii) added SAL boundary conditions for E,B, (iii) added Dipole Offset

--> include/Collective: (i) added variable Yes_sal, Nlayers_sal taken in inputfile and which define the absorbing layer at the boundaries, (ii) read DipoleOffset from input file

--> include/Particles3Dcomm: (i) added routine getRho() to compute total net charge over the box 

--> include/PSKOutput: (i) in parts for Ball, Bx, By, Bz use the variable getBx,getBy,getBz and NOT getBtot otherwis we sum two tims the xtrnal dipole!!! (this was giving problems in restart)

--> main/Collective: (i) read Yes_sal, Nlayers_sal, DipoleOffset (ii) added warning to print in output file "width of the SAL must be greater than few alfven wave times"

--> main/iPic3Dib: some changes, no major ones

--> particles/Particles3Dcomm: (i) added getRho at the end, (ii) added some comments and prints to check which routines the code was passing trough

--> particles/Particles3D: (i) in mover_PC_AoS_Relativistic changed a variable called "c" with "cc" since c is already speed of light, was giving warning, (ii) slightly changed populate_cell_with_particles which is using the n,v,T of solar wind at t=0 to create particles - the part using rotor is deprecated, (iii) populate_cell_with_particles_interp is a new routine which I added and is supposd to use a given value of n,v,T to repopulate particles in cell - which are not a priori the sw values - however this routine is no more used in th code, deprecated after testing, (iv) repopulate_particles has not been strongly modified after all

--> utility/Basic: (i) added routine extrapolate() to compute n,v,T given their values in neighbours cells close to the boundary - deprecated after testing 
 

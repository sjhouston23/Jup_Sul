# 2Stream

The two-stream code is the electron code that assumes two streams of electrons, one flowing up, another flowing down. It has a large explanation at the top of the code.

The command to run the code is:
gfortran -o elect.x elect_jup_sec.f08 && ./elect.x

Another important code for the use of the 2-stream is readall_jup.f08 code:
gfortran -o read.x readall_jup.f08 && ./read.x [energy]

The electron energy bins are visualized in  ./common/input/xgridnew_ext.jupiter.

The main file is elect_jup_sec.f08, which I have edited so there is also an older version titled elect_jup_sec_OLD.f08. The newer version just allows for the use of two different atmospheres, the original atmosphere and the well-mixed atmosphere. This can be edited at line 433, 436, and 439.

There are a few important things to note from this code. Firstly, it has been mangled together by a handful of people. So it is not very coherent. I have tried editing it enough to at least have somewhat consistent syntax throughout, but it's still not great. Next, I'll go over the important input and output files, then how to run a monoenergetic electron beam.

To run the code, need to use the parameters file located in ./input/elect/jelect.par (line 123). Each line of this file can be deduced by looking at the code, as it reads it in beneath the opening of the file. The 4th line is where you define the ion precipitation energy. I have defined it here so it can automatically read in the correct files for lines 167-180. Earlier editions had all the file names read in from this file, but I have changed it so they're read in within the code so the energy can easily be changed automatically. These files include:

flux3d_JUP...dat - An output file that contains all of the phi up and phi down electron fluxes vs. altitude (r/s in this case) and for each electron energy bin.

airglow_JUP...dat - An output file that contains all of the "airglow" information. That includes UV emission and atmospheric ionization. So, it begins with Lyman-Werner band emission from each atmospheric species then shows electron impact ion production rates. All vs. altitude.

../Output/[energy]/2Str_Elect_[Fwd/Bwd]-Comb.dat - Input files from the ion precipitation code. This obviously goes backward a directory into the ion precipitation model, then goes into the output there to grab the correct electron production files.

ForwardElecProd...dat - An output file that contains the forward electron production at a number of different altitudes (see Fig. 5 in Houston et al. 2018).

You'll notice that the output files may have an extension of "-eq". I have just used this to establish the well-mixed (equally mixed) atmosphere files. And I just remove it when doing runs for atmosphere 1.

If you want to do a JEDI run, you need to put the combined 2Str JEDI files into a directory in ../Output/ with some energy (I always used 111 since I don't actually run that energy). Then it can read those files in just as it would any other file and output files with "111" attached to it.

This is not all the 2-stream code can do, but the rest is accomplished with the readall_jup.f08 code, which I'll get to in a minute.

First I want to discuss how to run a monoenergetic electron beam. This is done with the very similar code named elect_jup_sec_ME.f08, which uses the parameters file ./input/elect/jelect_ME.par. This code can also use both atmospheres by changing lines 420, 422, and 426 and is only found in the Jup_Sul/2Stream directory. All of the monoenergetic parameters are defined in lines 458-476. Where I have commented it to try and make it clear. Line 461, the if(I.eq.###)then, statement is what determines the energy. Change ### to the electron energy bin you want your monoenergetic electron beam to be. I have commented in common energies on lines 465-470. PHIINF(I) stands for Phi (the flux) at infinity (i.e., the top of the atmosphere). It naturally assumes number of electrons/cm^2/s/eV. A lot of times we want the energy flux in terms of a power input of mW/m^2. I have made the conversion easy by outputting the required PHIINF value to input exactly 1 mW/m^2. I also documented the common values on lines 465-470.

--------------------------------------------------------------------

readall_JUP.f08 instructions:

This is the code that makes sense of all of the 2-Stream output files. I once again have updated this file, so there is an _OLD version. This program takes in all of the 2-Stream output and calculates different things. It also uses ./input/elect/jelect.par to find the energy and the SZA (Solar Zenith Angle, which I always have set to SE because I don't include solar photoionization - only considering the nightside case). It begins by reading in the 2-Stream output files and also some of the ionization files from the ion precipitation code. Then it outputs a variety of things:

secprd...dat - This is an output file containing all of the atmospheric ionization production rates.

escapeflux...dat - This is an output file containing the upward electron flux at 3000 km. (Fig. 11 Houston et al. 2018)

escapecurrent...dat - This is an output file that contains just 3 values that are commonly referenced. That is the escape flux in eV/cm^2/s, the escape electron flux in electrons/cm^2/s, and the escape current in A/m^2. (Tab. 1 in Houston et al. 2018).

ALTvsFLUX...dat - This is an output file containing the upward electron flux vs. altitude for a handful of different electron energies. (Fig. 12 Houston et al. 2018).

LymanWerner...dat - This is an output file containing the Lyman-Werner band emission vs. altitude.

LWSpectrum...dat - Very similar to the previous file... I don't actually know why I output this twice. But both this file is used in the program titled LWSpec.f08.
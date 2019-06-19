# Jup_Sul

The top portion is for running the sulfur precipitation code locally. The bottom discusses the sulfur precipitation code that has been edited for the computing cluster (CRC).

I have already run, and gotten results for, a bunch of different energies. The only reason to re-run any of the energies is if new cross-sections become available or you want a new energy that I haven't run yet. When using the computing cluster, I would run ~20,000 ions at one time.

Most codes I have written will have a description of what they do at the very top of the code.

------------------------------------------------------------------

To create the formatting module:
gfortran -c formatting.f08
This creates formatting.o to be used in the next command.

The following command is used to run the JupOxyPrecip.f08 program:
gfortran -o precip.x JupSulPrecip.f08 ranlux.f08 CollisionSim.f08 EjectedElectron.f08 EnergyLoss.f08 formatting.o && ./precip.x

JupSulPrecip.f08 is the main sulfur precipitation program. It loops through all of the desired energies (variable IonEnergyNorm or IonEnergyJuno, defined in the "Data Declaration" section) and calculates the results for them. The energies chosen to run is determined by the EnergySwitch on line 286. The number of ions it considered is defined in the variable "nIons", currently on line 301. The following line (302), defining "trial", is what separates the output files for different runs. Each run is labeled with the number that is set here. The next line (303) begins the do-loop for which energy to run. Each energy is numbered a few lines above this one, so you can set it to do a single energy. Or you can loop it to do multiple energies.

This code is well commented, so most of the things within it can be figured out by reading all of the comments. All of the output is put into the ./Output/ directory in the associated energy directory. The subroutines are the following:


ranlux.f08 - The random number generator. It is seeded by calling rluxgo (e.g. on line 313) and it generates an array of random numbers by calling ranlux (e.g. on line 315). In my case, I have seeded the random number generator with whatever the trial number is.

CollisionSim.f08 - Determines which collision will occur. It creates a probability distribution for a given energy and charge state for each type of collision then uses the RNG to choose a collision. Since the collision is chosen in this subroutine, it also determines the number of electrons ejected, whether dissociation will occur, etc... The integral cross-sections are read in on lines 218-247 of JupSulPrecip.f08.

EjectedElectron.f08 - If an electron is ejected, this subroutine is called so it can figure out the energy and angle of the ejected electron. It uses singly differential cross-sections for both energy and angle (these cross-sections are read in at the beginning of EjectedElectron.f08).

EnergyLoss.f08 - Determines the amount of energy that is lost after the collision occurs. It separates all of the collision energy models and calculates an energy loss for a specific energy.

formatting.o - Has all of the formatting code in it so it doesn't take up a bunch of space at the bottom of the code and is reused for various programs.

------------------------------------------------------------------

To run on the computing cluster, you first need to get an account there (https://crc.ku.edu/hpc/account).

In general, you can learn more about how to use the cluster here: https://crc.ku.edu/hpc/how-to

The code I have made for the cluster is JupSulPrecipCRC.f08. To compile this you use the following command:
gfortran -o precip.x JupSulPrecipCRC.f08 ranlux.f08 CollisionSim.f08 EjectedElectron.f08 EnergyLoss.f08 formatting.o

Notice that this command doesn't actually execute the compiled executable. This is because rather than executing the file once, you want to execute it 100 or 200 times. The number of ions currently defined in the program is 100. So if you execute it 200 times, you'll get the results for 20,000 ions. But there are a couple things to note before this can be done.

After you get access to the cluster you do the following:

1)Make sure you can ssh into the cluster, e.g. ssh s500h396@hpc.crc.ku.edu

2)Before moving all of the subroutines over, you want to change anything that writes out to the screen (i.e. write(*,*)) to write out to file 206 (i.e. write(206,*)) in all of the subroutines that are being moved over. JupSulPrecipCRC.f08 already does this, so look there for examples.

3)You then want to secure copy (scp) all of the relevant files onto the cluster. These need to be put in your home directory (~). So first create a directory titled "Jup_Sul" then scp all of the subroutines over into that directory (e.g. scp ~/Jup_Sul/JupSulPrecipCRC.f08 s500h396@hpc.crc.ku.edu:~/Jup_Sul/). (NOTE: Don't worry about moving executable files over, they will all have to be compiled on the cluster since it's a Linux system. They compile differently than on a Mac or Windows, or even different version of Linux.) Also, you need to move formatting.f08 over onto the cluster.

4)Assuming all of the files have been moved to the cluster, you then want to create the directories for the output files. I direct JupSulPrecipCRC.f08 to write the output files into ../scratch/Jup_Sul/Output/[energy]/. So you need to go into the scratch folder and create a directory titled Jup_Sul/, then in there and create a directory titled Output/, then in there create directories for whatever energies you're going to be running.

5)After this is all done, you should be set up to do a run. I have written this program such that when you call the executable (i.e. precip.x) it looks for a number after the execution command so it can initialize the RNG (e.g. ./precip.x 23). This allows me to use different numbers for each execution of the program so they're not all the same. From here you first need to create a new formatting file. This is done by running the command "gfortran -c formatting.f08" which creates a formatting.o file. Then you need to recompile all the programs for the executable.

6)You need "submission script" to actually use the cluster. I have created one in the same directory as this file titled "run.sh". You can scp this file over. If you look at it, you can see the line that says #SBATCH --array=1-200. This says that I want to do 200 different iterations. It then chooses a random number with R=$RANDOM and plugs that variable into the line ./precip.x $R. I then output all of those random numbers into ../scratch/Jup_Sul/Output/Trials.dat so they can be kept track of (I write these numbers out multiple times to different places).

7)To submit the run.sh script, you put the script in the same directory as the model on the cluster, ~/Jup_Sul/ and then run the command: sbatch run.sh
More of these commands can be found on the crc website and I've posted a screenshot of some common commands within this directory, titled "HRC Commands".

8)After all of the processes are finished running, the files need to be combined into a single file. This is done using the Combine_Rec.f08 program. So scp that file over to the cluster, ensure the energies match up with what you just ran with the JupSulPrecipCRC.f08 code, then this can be run by using the command: gfortran -o comb.x Combine_Rec.f08 formatting.o && ./comb.x
This program works by reading in the random numbers used and then reading in all of the output files which are attached to those random numbers and creating single files out of them. It places all of the combined files into a directory within ~/Jup_Sul/Output/[energy]. So you need to make sure you have those directories created, too.

9)These final combined files can then be moved back to your local computer.

This same process is used for oxygen precipitation.

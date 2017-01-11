#########################################
### SPPLot input file - Version: 140314
#########################################
### Jack Morford
#########################################

### SPPlot follows a 'simple' routine. If you choose to apply no flag table to your data 
# and specify over 5 baselines to be plotted, then the program will simply take your UV 
# dataset from your AIPS catalogue and proceed to append the visibilities into the memory.
# Once each source is read into the memory (assuming you choose multiple sources) then you 
# have the option to save the visibilities into a picklefile, for example if you wish to 
# plot the PHASE next time, this will save time compared to appending the visibilities into 
# the memory arrays again. The function plotspec() then pulls out the recquired infomation 
# from the visibility arrays such as your desired source and stoke(s) and proceeds to plot 
# them baseline by baseline. The program will then merge each single plot of each source 
# and stoke into a single pdf file hence the number of merged pdf files will depend upon 
# the number of baselines you have chosen. You can now choose to further combine these 
# merged baseline plots into a final single plot file.
### If you have selected to plot 5 or fewer baselines the program will automatically UVCOP 
# your UV data into a new dataset that includes only your specified baselines. This is to 
# save time when appending visibilities into the memory, then the program will continue as 
# described above.
### When applying a particular flagversion to the data, the SPPlot will UVCOP your dataset 
# regardless of the number of baselines in order to apply your chosen flag table. However, 
# again if 5 baselines or fewer are specified then the new UVCOP'd dataset will contain
# only those baselines specified.

### This program can be run on both single and multi source files and you may select the 
# sources you wish to plot.  

AIPS.userno = 3216           	# The AIPS user number the data is on.
Name = 'J0958+65'             	# The uvdata name of the catalogue.
Klass = 'SPLAT'                 # The uvdata klass of the catalogue.
Disk = 1                        # The uvdata disk number of the catalogue (integer).
Seq = 1


#path2folder = '/import/cobrasarchive1/data/COBRAS/LEGACY/L-band-May13/calibration2/'		
#path2folder = '/import/cobrasscratch1/dmf/python/lev2/'
#path2folder = '/import/cobrasarchive1/data/COBRAS/LEGACY/L-band-Jan14/calibration/flagging/'
path2folder = '/import/cobrasscratch1/dmf/python/VLBI/'
# Choose output directory for your plots to be stored. 

picklepath = str(path2folder)+'picklefiles'
# Choose output directory for your pickled dictionaries to be stored.

picklelookfor = False		      # If True the program will look for any previously created 
                                      # picklefile of the same outfilename+Name+flagver+outclass

picklesave = False		      # If False, the visibilities won't be saved - this will 
                                      # ensure no confusion with future runs *RECOMMENDED*
                                      # If True, the program will pickle your visibility data into 
                                      # a file - this is intended to save time if all baselines are 
                                      # choosen and further plots are wanted to be produced of the 
                                      # same data e.g. Amp, Phase, different stokes. NB if True is 
                                      # selected, keep track of the pickle file. For a future run 
                                      # for a different dataset (perhaps with a different flag table) 
                                      # then delete picklefile or set picklelookfor = False.

choosesources = 'all'		      # Please select either 'all' or 'choose'

specifysources = ['J0958+65']	      # Must be strings in python list format.


choosebaselines = 'all'	              # Please select either 'all' or 'choose'
	
specifybaselines = [[2,7]]	      # Specify your baselines in a tuple format, i.e. [[1,2],[1,5],etc]
                                      # NB: if 5 baselines or fewer are chosen the program will UVCOP 
                                      # the dataset for your chosen baselines - this is intended to save 
                                      # time when appending the visibilities into the memory. If more than 
                                      # 5 baselines are selected, the entire UV datafile including all 
                                      # baselines will be read into the memory and pickled into a pickle file.


outfilename = 'J0958+65'      	      # This is the base name of your output plots.
                                      # NB: output file will have structure of outfilename+sourcename+visstart
                                      # -visend+Amp/Phas+stoke+baseline+flagtablenumber.

deleteprev = True		      # Select True if you wish to delete any previously created single 
                                      # baseline UV datafiles within your AIPS disk 
                                      # THIS IS A GOOD IDEA FOR MULTIPLE RUNS OF SERPent!
                                      # Otherwise select False if you have previously created a UV data set 
                                      # that includes your desired baselines along with the correct flagtables 
                                      # applied!
                                      # NB: if False is selected, the program will tell you if the UV dataset 
                                      # in question doesn't include your desired baselines!

flagver = 0 			      # Please specify which FG version you would like to apply to the data - 
                                      # if 0 is selected this implies that no flag table is selected.

stokes = ['RR','LL']		      # Must be strings seperated by commas!

# Plotting Options:

scale = 'log'	         	      # Choose whether to scale the amplitudes linearly ('linear'), 
                                      # logarithmically ('log') or min/max within 3 std of mean ('std') - 
                                      # a log scale may show greater detail if the differences in min and max 
                                      # amplitudes are very large. If 'std' is chosen, the max and min amplitudes 
                                      # will be taken to be the mean +/- three stand deviations.

scale_over_all_IFs = 'no'	      # Choosing to scale over all IFs gives a true representation of your data 
                                      # but since RFI can vary severly from IF to IF it is likely that your plots 
                                      # will show less detail for those IFs less affected by RFI. If 'no' the 
                                      # program will scale each IF individually according to their own maximum 
                                      # and minimum amplitudes - this will give a good indication of RFI within 
                                      # individual IFs but note that the colorbar scale will correspond only to 
                                      # the last IF. Underneath each IF there is a min and max amplitude value 
                                      # for each IF - this will give you an idea of the differences between each 
                                      # IF if 'no' is selected. 
colourscheme = 'jet'                  # This allows you to choose the colour scheme adopted for the colour 
                                      # scaling of your plots. Currently supported options are jet, spectral, 
                                      # gist_rainbow, cool, summer or winter. To see what these colour schemes
                                      # look like see http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
                                      # This will default to the 'jet' colourscheme if none is specified.


amporphas = 'A'			      # Specify your plots to display either amplitude: 'A' or phase: 'P'. 
                                      # NB: Might be pretty buggy with phase plots.

timeperpage = 3200 		      # Set the amount of visibilities to use per page - i.e. the amount of time 
                                      # (y-axis) to squeeze onto a sheet of A4. Something between 1000-2000 is fine.

IF = 2				      # How many IFs are present in your dataset i.e. for e-Merlin L-Band = 8, 
                                      # C-Band = 16.

IF_start = 1 			      # Choose starting IF to plot.
IF_end = 2			      # Choose last IF to plot.


requestfullmerge = True               # Choose to merge all baseline plots into a single file.

removebslfile = True                 # If requestfullmerge = True, you can choose to remove the single baseline
                                      # merged files leaving only the a single (all baseline) pdf. 

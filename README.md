# Major_Project_Docking_Profiler
Program description: Docking_Profiler.py generates overview plots from a VLS-docking run with ICM, using a SDF-file containing ICM-VLS scores as input

Input requirements: An SDF file in which the first column contains molecule names and is called "NAME" and in which all other columns must only contain numerical values. A Samle_input.sdf is provided.

The user can specify the input and output file names and a threshold value at the top of the Docking_Profiler.py script.
The parameters that are to be plotted are specified at the end of the Docking_Profiler.py script.

Output: The program generates a PDF file (Plot_Docking_Profile.pdf) containing all plots. Furthermore the corresponding tables are written out (scorelist.tsv, scorelist_sorted.tsv, scoreist_sorted_topX.tsv and averages.tsv).

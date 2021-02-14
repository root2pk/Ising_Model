# Ising_Model
Monte-Carlo simulation of Ising model using Glauber and Kawasaki Dynamics

  • The file ’CP1.py’ contains the code to simulate the Ising model for one temperature value.
Usage : python CP1.py N T
Where N is the length of one side of the square lattice and T is the temperature.
The system them prompts the user to enter 1 for Glauber or 2 for Kawasaki Dynamics.
The program writes temperature, <E>, <|M|>, C, σC and χ to the file ”ObservableT.txt”

  • The file ’CP1 loop.py’ contains the code to simulate the Ising model for a range of temperature values.
Usage : python CP1 loop.py N T start T end number of temps
Where N is the length of one side of the square lattice, T start is the start temperature, T end is the end
temperature and number of temps is the number of equidistant temperature values between the starting
and ending temperature(including both). The system them prompts the user to enter 1 for Glauber or 2 for Kawasaki Dynamics.
The program writes temperature, <E>, <|M|>, C, σC and χ to the file ”Observables.txt” for each
temperature iteration.

  • The file ’methods.py’ contains all the functions used in the main program

  • The file ’plot program’ uses the data columns in ”Observables.txt” to plot the graphs for <E>, <|M|>,
C(with errors),and χ.

  • This program uses the bootstrap resampling method to calulate the error in C

#----------------------------------------------------------------------------------

#SimEfficienciesCalculation
#Written by Andrew Boldy
#University of South Carolina
#Spring 2025

#----------------------------------------------------------------------------------

#Plots the electron momentum for a given non-mixed generator filelist
#Must be run using C++ so compile mode in order to rectify issues that are generated when using common_cuts.hh

#----------------------------------------------------------------------------------


# --- Database Interaction ---
import DbService
# Establish a connection and retrieve simulation efficiencies from the database.

# Database Tools
db_tool = DbService.DbTool() #creates the database tool 
db_tool.init() #initializes the database tool 

# Define arguments as a list of strings for the database query
query_arguments = [
    "print-run", 
    "--purpose", "Sim_best",
    "--version", "v1_1",
    "--run", "1430",
    "--table", "SimEfficiencies2",
    "--content"
]

# Execute the database query
db_tool.setArgs(query_arguments) #sets the aguments for the database tool 
db_tool.run() #runs the database tool given these arguments given

# Store the raw result for further processing
rr = db_tool.getResult() #when the program has been run, we store the raw result as a variable called rr

# Fill varaibles associated with muon stops in target
lines = rr.split("\n") #split lines up

#initialize variables for rate and for target_stopped_muons_per_pot
rate = 1.0 
target_stopped_muons_per_pot = 0.0

#create a dictionary to store efficiencies
efficiencies = {}

for line in lines:
    words = line.split(",")
    if not line.strip(): #skipping empty lines
        continue
    if len(words) < 4: #preventing indexing errors
        continue

    name = words[0].strip()
    value = float(words[3])

    efficiencies[name] = value

    if name in ["MuminusStopsCat", "MuBeamCat"]:
        rate*=value
        
target_stopped_muons_per_pot = rate * 1000

#Output
print("\n--- Simulation Efficiencies ---")
for key, val in efficiencies.items():
    print(f"{key:30s} : {val:.6e}")

print("\n--- Derived Quantities ---")
print(f"Target stopped muons per POT: {target_stopped_muons_per_pot:.6e}")

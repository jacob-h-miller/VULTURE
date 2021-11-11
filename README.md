# VULTURE
MATLAB code for our biorefinery optimizer: VFA Upgrading to Liquid Transportation fUels Refinery Estimation

FuelCriteriaMain.m is the main function. A list of C2-8 VFA mole fractions (muast add up to 1) must be input into the vector "VFA_molfracs" (line 13) to initialize the code. Possible input acids (in the order used for the input vector) are: acetic, propionic, isobutyric, butyric, isovaleric, valeric, hexanoic, heptanoic, and octanoic. All othe MATLAB files in the directory are functions called by FuelCriteriaMain.m or sub-functions called by other functions.

Of the excel workbooks, Biorefinery_hydrocarbon_database.xlsx and Prediction_FP_Database.xlsx calculate and store the properties of neat n-alkanes and alcohols, respectively. Information about the sources and/or formulas used for these calculations is also contained in these spreadsheets. Neat_Petrofuel_Properties.xlsx lists the properties of RBOB gasoline, clay-treated diesel, and Jet A. Molecule_vapor_pressures.xlsx lists vapor pressures of all molecules used in VULTURE as a function of temperature and sources of these calculations. VPs_MWs.xlsx lists the molecular weights of all molecules used in VULTURE.

# main.cpp

Entry point for the program. Checks for runtime arguments and passes along to SABREsim class. Calls SABREsim::Run().

# SABREsim.cpp

'SABREsim' is the main driver class for the SABRE simulation. It is responsible for:

- Reading and parsing the simulation configuration file
- Initializing detectors, beam properties, and physics models
- Running the simulation for different reaction types (2-body, 3-body, 4-body)
- Managing input/output files
- Collecting and printing summary statstics

This class acts as the central coordinator that connects all major components of SABREsim.

In the constructor, the configuration file is read in to memory and class values for input and output files are set. RootWriter object is accordingly set up.

SABREsim::Run() is called in main.cpp and is responsible for actually setting up output to files defined in config file, initializes SABRE array, beam spot characteristics, and energy loss and angular straggling models. Then, according to the kinX specificed in the configuration file, SABREsim::Simulate2body(), SABREsim::Simulate3body(), or SABREsim::Simulate4body() is called. Once the SABREsim::Simulate2/3/4body() finishes the prepared ROOT TTree is written to disk and the simulation statistics are reported and printed out via SABREsim::PrintSummary().

SABREsim::CleanUp() is called in the destructor and just cleans up pointers.

SABREsim::Simulate2body() creates an object from class det2mc.cpp. This object handles the actual process for 2body cases and will be described in the de2mc.cpp class section.

SABREsim::Simulate3body() creates an object from class det3mc.cpp. This object handles the actual process for 3body cases and will be described in the de3mc.cpp class section.

SABREsim::Simulate4body() creates an object from class det4mc.cpp. This object handles the actual process for 4body cases and will be described in the det4mc.cpp class section.

# RootWriter.cpp

'RootWriter' is the class responsible for writing TTree output for all cases.

# HistoManager.cpp

# plot2mc.cpp

# plot3mc.cpp

# plot4mc.cpp

# det2mc.cpp

# det3mc.cpp

# det4mc.cpp

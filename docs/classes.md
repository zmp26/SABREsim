# main.cpp

Entry point for the program. Checks for runtime arguments and passes along to SABREsim class.

# SABREsim.cpp

'SABREsim' is the main driver class for the SABRE simulation. It is responsible for:

- Reading and parsing the simulation configuration file
- Initializing detectors, beam properties, and physics models
- Running the simulation for different reaction types (2-body, 3-body, 4-body)
- Managing input/output files
- Collecting and printing summary statstics

This class acts as the central coordinator that connects all major components of SABREsim.



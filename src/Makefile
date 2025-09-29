CC = g++
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
CFLAGS = -Wall -std=c++17 $(ROOTCFLAGS)
EXECUTABLE = SABREsim

OBJECTS = Rotation.o Vec3.o SABRE_Detector.o Beamspot.o SABRE_EnergyResolutionModel.o SABRE_DeadLayerModel.o ConsoleColorizer.o SABRE_ChannelMap.o TargetEnergyLoss.o SABREsim.o det2mc.o main.o

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(ROOTLIBS)

Rotation.o: Rotation.cpp Rotation.h
	$(CC) $(CFLAGS) -c $<

Vec3.o: Vec3.cpp Vec3.h
	$(CC) $(CFLAGS) -c $<

SABRE_Detector.o: SABRE_Detector.cpp SABRE_Detector.h
	$(CC) $(CFLAGS) -c $<

Beamspot.o: Beamspot.cpp Beamspot.h
	$(CC) $(CFLAGS) -c $<

SABRE_EnergyResolutionModel.o: SABRE_EnergyResolutionModel.cpp SABRE_EnergyResolutionModel.h
	$(CC) $(CFLAGS) -c $<

SABRE_DeadLayerModel.o: SABRE_DeadLayerModel.cpp SABRE_DeadLayerModel.h
	$(CC) $(CFLAGS) -c $<

ConsoleColorizer.o: ConsoleColorizer.cpp ConsoleColorizer.h
	$(CC) $(CFLAGS) -c $<

SABRE_ChannelMap.o: SABRE_ChannelMap.cpp SABRE_ChannelMap.h
	$(CC) $(CFLAGS) -c $<

TargetEnergyLoss.o: TargetEnergyLoss.cpp TargetEnergyLoss.h
	$(CC) $(CFLAGS) -c $<

SABREsim.o: SABREsim.cpp SABREsim.h
	$(CC) $(CFLAGS) -c $<

det2mc.o: det2mc.cpp det2mc.h
	$(CC) $(CFLAGS) -c $<

main.o: main.cpp SABREsim.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
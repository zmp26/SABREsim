CC = g++
CFLAGS = -Wall -std=c++17
EXECUTABLE = SABREsim

OBJECTS = Rotation.o Vec3.o SABRE_Detector.o Beamspot.o SABRE_EnergyResolutionModel.o SABRE_DeadLayerModel.o ConsoleColorizer.o main.o

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^

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

main.o: main.cpp SABRE_Detector.h Vec3.h Rotation.h Beamspot.h SABRE_EnergyResolutionModel.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

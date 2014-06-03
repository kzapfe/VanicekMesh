FUENTESTest=QuantumConstants.hpp BoltzmanforNelson.hpp
FUENTES=QuantumConstants.hpp BoltzmanforNelson.hpp PropagatorFreeParticle.hpp
LIBSGEN=-lgsl -larmadillo

TestBoltz01.x: TestingDistro01.cpp $(FUENTESTest)
	g++ -lgsl  $< -o $@

TestPropagator.x: ITNUMM_VetH_01.cpp $(FUENTES)
	 g++ $(LIBSGEN)  $< -o $@

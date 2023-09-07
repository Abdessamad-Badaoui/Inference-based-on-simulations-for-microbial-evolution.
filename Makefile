all: atreyu_forward_simulator atreyu_forward_simulator_ct gillespie

gillespie: gillespie.cpp forward_simulator.h
	g++ -Wall -O3 -ffast-math $< -o $@

atreyu_forward_simulator: atreyu_forward_simulator.cpp forward_simulator.h
	g++ -Wall -O3 -ffast-math $< -o $@

atreyu_forward_simulator_ct: atreyu_forward_simulator.cpp forward_simulator.h
	g++ -Wall -O3 -ffast-math -DCONSTANT_MUTATION_RATE $< -o $@

.PHONY: clean tests

tests:
	./run_comparaison_tests.py
	./run_multisegments_tests.py
	./run_replication_tests.py

clean:
	rm -f compare.test.*

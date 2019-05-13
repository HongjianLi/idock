CC=g++ -O3

bin/idock: obj/io_service_pool.o obj/safe_counter.o obj/array.o obj/atom.o obj/scoring_function.o obj/receptor.o obj/pka.o obj/residue.o obj/ligand.o obj/result.o obj/random_forest.o obj/random_forest_x.o obj/random_forest_y.o obj/stopwatch.o obj/main.o
	${CC} -o $@ $^ -pthread -lstdc++fs -L${BOOST_ROOT}/lib -lboost_program_options

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++17 -DNDEBUG -Wall -I${BOOST_ROOT}/include

clean:
	rm -f bin/idock obj/*.o

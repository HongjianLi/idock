CC=clang++ -std=c++11 -O2
NVCC=nvcc -std=c++11 -use_fast_math

all: bin/idock_cp bin/idock_cu bin/idock_cl src/kernel.fatbin

bin/idock_cp: obj/io_service_pool.o obj/safe_class.o obj/array.o obj/scoring_function.o obj/atom.o obj/receptor.o obj/ligand.o obj/random_forest.o obj/random_forest_x.o obj/random_forest_y.o obj/log.o obj/main_cp.o obj/kernel.o
	${CC} -o $@ $^ -pthread -L${BOOST_ROOT}/lib -lboost_system -lboost_program_options -lboost_filesystem

bin/idock_cu: obj/io_service_pool.o obj/safe_class.o obj/array.o obj/scoring_function.o obj/atom.o obj/receptor.o obj/ligand.o obj/random_forest.o obj/random_forest_x.o obj/random_forest_y.o obj/log.o obj/main_cu.o obj/source_cu.o
	${CC} -o $@ $^ -pthread -L${BOOST_ROOT}/lib -lboost_system -lboost_program_options -lboost_filesystem -L${CUDA_ROOT}/lib64 -lcuda -lcurand

bin/idock_cl: obj/io_service_pool.o obj/safe_class.o obj/array.o obj/scoring_function.o obj/atom.o obj/receptor.o obj/ligand.o obj/random_forest.o obj/random_forest_x.o obj/random_forest_y.o obj/log.o obj/main_cl.o obj/source_cl.o
	${CC} -o $@ $^ -pthread -L${BOOST_ROOT}/lib -lboost_system -lboost_program_options -lboost_filesystem -L${ICD_ROOT}/bin -L${AMDAPPSDKROOT}/lib/x86_64 -L${INTELOCLSDKROOT}/lib64 -lOpenCL

obj/main_cu.o: src/main_cu.cpp
	${CC} -o $@ $< -c -I${BOOST_ROOT} -I${CUDA_ROOT}/include

obj/main_cl.o: src/main_cl.cpp
	${CC} -o $@ $< -c -I${BOOST_ROOT} -I${ICD_ROOT}/inc -I${AMDAPPSDKROOT}/include -I${INTELOCLSDKROOT}/include

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -I${BOOST_ROOT}

src/%.fatbin: src/%.cu
	${NVCC} -o $@ $< -fatbin -gencode arch=compute_35,code=compute_35

clean:
	rm -f bin/idock_cp bin/idock_cu bin/idock_cl src/kernel.fatbin obj/*.o

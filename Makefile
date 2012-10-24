CXX=mpicxx -g -Wall -pedantic
INCS= -I. -I/home/tongfei/soft/PETSc/include -I/home/tongfei/soft/PETSc/arch-linux2-c-debug/include
LIBS = -L/home/tongfei/soft/PETSc/arch-linux2-c-debug/lib 

OBJS =   solver_petsc.o \
	 solver_lapack_cf.o \
	 octree.o \
	 DATA.o \
	 EllipticSolverBase.o \
	 MHDSolver.o \
	 HyperbolicSolverBase.o \
	 CompressibleSolver.o \
	 TimeIntegratorBase.o \
	 VerletScheme.o \
	 Controller.o \

all: libLSGFD.a LSGFD

LSGFD:

solver_petsc.o: solver_petsc.cpp solver_petsc.h
	${CXX} -c solver_petsc.cpp ${INCS}
solver_lapack_cf.o: solver_lapack_cf.cpp solver_lapack_cf.h
	${CXX} -c solver_lapack_cf.cpp ${INCS}
octree.o: octree.cpp octree.h
	${CXX} -c octree.cpp ${INCS}
DATA.o: DATA.cpp DATA.h
	${CXX} -c DATA.cpp ${INCS}
EllipticSolverBase.o: EllipticSolverBase.cpp EllipticSolverBase.h
	${CXX} -c EllipticSolverBase.cpp $(INCS)
MHDSolver.o: MHDSolver.cpp MHDSolver.h
	${CXX} -c MHDSolver.cpp ${INCS}
HyperbolicSolverBase.o: HyperbolicSolverBase.cpp HyperbolicSolverBase.h
	${CXX} -c HyperbolicSolverBase.cpp ${INCS}
CompressibleSolver.o: CompressibleSolver.cpp CompressibleSolver.h
	${CXX} -c CompressibleSolver.cpp ${INCS}
TimeIntegratorBase.o: TimeIntegratorBase.cpp TimeIntegratorBase.h
	${CXX} -c TimeIntegratorBase.cpp ${INCS}
VerletScheme.o: VerletScheme.cpp VerletScheme.h
	${CXX} -c VerletScheme.cpp ${INCS}
Controller.o: Controller.cpp Controller.h
	${CXX} -c Controller.cpp ${INCS}

LSGFD:	LSGFD.cpp libLSGFD.a
	${CXX} -c LSGFD.cpp ${INCS}
	${CXX} -o LSGFD LSGFD.o ${INCS} -L. -lLSGFD ${LIBS} -lpetsc -llapack -lblas -lX11

libLSGFD.a: ${OBJS}
	ar cru libLSGFD.a ${OBJS}
	ranlib libLSGFD.a;

clean:
	rm -f ${OBJS} LSGFD

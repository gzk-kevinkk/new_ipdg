.SUFFIXES: .cxx

CPP		= /usr/bin/mpicc -E
CPPFLAGS	= -I/home/kevinkk/Desktop/app/phg/cur_phg/phg-0.9.6-20220323/phg-0.9.6/include  -I/usr/local/metis/include -I/usr/local/parmetis/include -I/usr/local/hypre/include -I/usr/include/suitesparse -D__PHG__
BUILD_CPPFLAGS	= -DMPI_NO_CPPBIND   -DVTK_DIR="\"\""  -DGZIP_PROG="\"/usr/bin/gzip\"" -DBZIP2_PROG="\"/home/kevinkk/anaconda3/bin/bzip2\""

CC		= /usr/bin/mpicc 
CFLAGS		= -g -O2 -Wall -Wmissing-field-initializers -fopenmp
BUILD_CFLAGS	= 
USER_CFLAGS	=

CXX		= /usr/bin/mpicxx 
CXXFLAGS	= -g -O2 -Wno-deprecated -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX -fopenmp
BUILD_CXXFLAGS	= 
USER_CXXFLAGS	=

FC		= /usr/bin/mpif90 
FCFLAGS		= -g -O2
USER_FCFLAGS	=

F77		= /usr/bin/mpif77 
FFLAGS		= -g -O2
USER_FFLAGS	=

LINKER		= /usr/bin/mpicxx 
LDFLAGS		= -L/home/kevinkk/Desktop/app/phg/cur_phg/phg-0.9.6-20220323/phg-0.9.6/src  -fopenmp -L/usr/local/parmetis/lib -L/usr/local/hypre/lib
USER_LDFLAGS	=

LIBS		= -lphg  -lmpfr  -lgmp  -lumfpack -lamd  -Wl,-rpath,/usr/local/hypre/lib  -lHYPRE -lparmetis  /usr/local/metis/lib/libmetis.a /usr/local/parmetis/lib/libmetis.a  -llapack  -lblas  -lm  -L/usr/lib/x86_64-linux-gnu/openmpi/lib -L/usr/lib/gcc/x86_64-linux-gnu/9 -L/usr/lib/gcc/x86_64-linux-gnu/9/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/9/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/9/../../.. -lgfortran -lm -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lquadmath -lpthread  -lgfortran -lm -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lquadmath -lpthread -Wl,-rpath,/usr/local/lib
USER_LIBS	=

AR		= /usr/bin/ar 
RANLIB		= /usr/bin/ranlib 
BUILD_SHARED	= 
BUILD_SHARED_LIBS = 
LIB_SUFFIX	= .a

.c.o:
	${CC} ${CFLAGS} ${CPPFLAGS} ${USER_CFLAGS} -c $*.c

.cxx.o:
	${CXX} ${CXXFLAGS} ${CPPFLAGS} ${USER_CXXFLAGS} -c $*.cxx

.cpp.o:
	${CXX} ${CXXFLAGS} ${CPPFLAGS} ${USER_CXXFLAGS} -c $*.cpp

.C.o:
	${CXX} ${CXXFLAGS} ${CPPFLAGS} ${USER_CXXFLAGS} -c $*.C

.f.o:
	${FC} ${FCFLAGS} ${USER_FCFLAGS} -c $*.f

% : %.o
	${LINKER} ${USER_LDFLAGS} ${LDFLAGS} -o $@ $*.o ${USER_LIBS} ${LIBS}

% : %.c
	${CC} ${CFLAGS} ${CPPFLAGS} ${USER_CFLAGS} -c $*.c
	${LINKER} ${USER_LDFLAGS} ${LDFLAGS} -o $@ $*.o ${USER_LIBS} ${LIBS}

% : %.cxx
	${CXX} ${CXXFLAGS} ${CPPFLAGS} ${USER_CXXFLAGS} -c $*.cxx
	${CXX} ${USER_LDFLAGS} ${LDFLAGS} -o $@ $*.o ${USER_LIBS} ${LIBS}

% : %.cpp
	${CXX} ${CXXFLAGS} ${CPPFLAGS} ${USER_CXXFLAGS} -c $*.cpp
	${CXX} ${USER_LDFLAGS} ${LDFLAGS} -o $@ $*.o ${USER_LIBS} ${LIBS}

% : %.C
	${CXX} ${CXXFLAGS} ${CPPFLAGS} ${USER_CXXFLAGS} -c $*.C
	${CXX} ${USER_LDFLAGS} ${LDFLAGS} -o $@ $*.o ${USER_LIBS} ${LIBS}

% : %.f
	${FC} ${FCFLAGS} ${USER_FCFLAGS} -c $*.f
	${LINKER} ${USER_LDFLAGS} ${LDFLAGS} -o $@ $*.o ${USER_LIBS} ${LIBS}

CC:
	@echo $(CC)

CPP:
	@echo $(CPP)

CXX:
	@echo $(CXX)

FC:
	@echo $(FC)

CPPFLAGS:
	@echo $(CPPFLAGS)

CFLAGS:
	@echo $(CFLAGS)

CXXFLAGS:
	@echo $(CXXFLAGS)

FFLAGS:
	@echo $(FFLAGS)

FCFLAGS:
	@echo $(FCFLAGS)

LDFLAGS:
	@echo $(LDFLAGS)

LIBS:
	@echo $(LIBS)

LINKER:
	@echo $(LINKER)

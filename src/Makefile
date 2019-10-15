
#  Written by C. Froese Fischer and Y. Qiu
#  Modified and extended by O. Zatsarinny
#  
#  Modified by C. Froese Fischer
#  January, 2006

.SUFFIXES: .f90 .mod

LIBS   =   -llapack -lblas
MODDIR = .

# Environment variables needed for compilation
#FC = gfortran
FC_FLAGS = -g -Mbounds
FC_LD =  
#LAPACK_DIR= 

install:  sphf_scf sphf_all

My_LIBS= spline_param.o		 spline_grid.o	    spline_galerkin.o      \
         spline_hl.o		 spline_integrals.o spline_moments.o       \
         hf_atomic_state.o	 hf_orbitals.o	    hf_energy_expression.o \
         factorial.o av_energy.o hf_inout.o	    hf_param.o             \
         gauss.o		 vbsplvd.o	    mrm.o                  \
         bxv.o			 hlm.o		    azl.o                  \
         mrk_cell.o		 rk.o		    density.o              \
         yval.o			 minty.o	    bhwf.o                 \
         rk_moments.o		 vinty.o	    moments.o              \
         routines.o		 bhl.o		    bvmv.o                 \
         quadr.o		 el_nl.o	    bvalu2.o               \
         lapack_interface.o	 plot_bsw.o

My_APPS= sphf.o                  get_atom.o         get_case.o             \
	 get_energy.o            mkgrid.o           get_estimates.o        \
         get_spline_param.o      define_spline.o    apply_orthogonality.o  \
         coulom.o                rotate.o           read_bsw.o             \
         refine_grid.o           hfall_nr.o         hf_matrix.o            \
         ab.o                    energy.o           summry.o               \
         read_hf_param.o

sphf_scf : $(My_LIBS) $(My_APPS) solve_scf.o
	$(FC) -o sphf_scf $(FC_LD) $(My_LIBS) $(My_APPS) solve_scf.o \
	-L ${LAPACK_DIR} -llapack -lblas
	
sphf_all : $(My_LIBS) $(My_APPS) solve_all.o
	$(FC) -o sphf_all $(FC_LD) $(My_LIBS) $(My_APPS) solve_all.o \
	-L ${LAPACK_DIR} -llapack -lblas
	
.f90.o:
	$(FC) -c $(FC_FLAGS)  $< -I $(MODDIR) -o $@

clean:
	-rm -f *.o core  *.mod



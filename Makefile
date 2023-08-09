# Makefile 
# Nom du compilateur
FC = f95

# Options de compilation: optimisation, debug etc...
OPT = -fdefault-real-8 -fbounds-check   
# nom de l'executable
EXE = step
# Options de l'edition de lien..
LINKOPT =  

# Defining the objects (OBJS) variables
SRC_DIR = src
OBJ_DIR = obj

OBJS =  \
       $(OBJ_DIR)/parameters.o \
       $(OBJ_DIR)/main.o \
       $(OBJ_DIR)/algo_tdma.o \
       $(OBJ_DIR)/initialisation.o \
       $(OBJ_DIR)/boundary_condition.o \

# Linking object files
exe :   $(OBJS)
	$(FC) $(LINKOPT) $(MODS) $(OBJS)  -o $(EXE) 
	
$(OBJ_DIR)/parameters.o : $(SRC_DIR)/parameters.f90
	$(FC) -c $(OPT)  $< -o $@

$(OBJ_DIR)/initialisation.o : $(SRC_DIR)/initialisation.f90
	$(FC) -c $(OPT)  $< -o $@

$(OBJ_DIR)/boundary_condition.o : $(SRC_DIR)/boundary_condition.f90
	$(FC) -c $(OPT)  $< -o $@

$(OBJ_DIR)/main.o : $(SRC_DIR)/main.f90
	$(FC) -c $(OPT)  $< -o $@

$(OBJ_DIR)/algo_tdma.o : $(SRC_DIR)/algo_tdma.f90
	$(FC) -c $(OPT)  $< -o $@


# Creating object directory if it doesn't exist
$(shell mkdir -p $(OBJ_DIR))


# Removing object files
clean :
	/bin/rm -f $(OBJS) $(EXE)  *.mod  ini  velo* vitesse* *dat  *gif mesh



EXE_NAME = dio
EXE_DIR = ./exe
SRC_DIR = ./src
MOD_DIR = ./src/mod
OBJ_DIR = ./src/obj
SRC_FILE_PREFIX = DIRHBZ

SOURCES = $(wildcard $(SRC_DIR)/*.f90)
# SOURCES = 	./src/dio.f90\
# 		./src/module_constants.f90\
# 		./src/module_globals.f90 \
# 		./src/module_inout.f90\
# 		./src/module_force.f90\
# 		./src/module_basis.f90\
# 		./src/module_mathmethods.f90\
# 		./src/module_nucleus.f90\
# 		./src/module_field.f90\
# 		./src/module_matrix.f90\
# 		./src/module_dirac_equation.f90\
# 		./src/module_occupation.f90\
# 		./src/module_density.f90\
# 		./src/module_expectation.f90\
# 		./src/module_constraint.f90\

OBJECTS = $(patsubst $(SRC_DIR)/%.f90, ${OBJ_DIR}/%.o, $(SOURCES))


#default: gfortran
default:  gfortran #ifort

# compiled by gfortran
gfortran: FC = gfortran
gfortran: FFLAGS = -O2 -J ${MOD_DIR}  #-Wall -g
gfortran: printConfiguration ${EXE_NAME} printEndInformation

# compiled by ifort
ifort: FC = ifort
ifort: FFLAGS = -O2 -module ${MOD_DIR}
ifort: printConfiguration ${EXE_NAME} printEndInformation

printConfiguration:
	@echo "=============Compiling with ${FC}====================="
	@echo "src path: ${SRC_DIR}/"
	@echo "mod path: ${MOD_DIR}/"
	@echo "obj path: ${OBJ_DIR}/"
	@echo "exe path: ${EXE_DIR}/${EXE_NAME}"
	@echo "------------------------------------------------------"
	
printEndInformation:
	@echo "------------------------------------------------------"
	@echo -e "\033[32mCompilation Finished ! \033[0m "
	@echo -e "src path: \033[32m ${SRC_DIR}/ \033[0m"
	@echo -e "exe path: \033[32m ${EXE_DIR}/${EXE_NAME} \033[0m"
	@ENDPrintEndInformation=]]]]]] # make brackets can be matched

${EXE_NAME}:${OBJECTS} | ${EXE_DIR}
	@echo "compiling $@ ......"
	${FC} ${FFLAGS} -o ${EXE_DIR}/${EXE_NAME} $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | ${OBJ_DIR}  ${MOD_DIR}
	@echo "compiling  $@ ......"
	$(FC) $(FFLAGS) -c $< -o $@ 

${OBJ_DIR}:
	mkdir -p ${OBJ_DIR}

${MOD_DIR}:
	mkdir -p ${MOD_DIR}

${EXE_DIR}:
	mkdir -p ${EXE_DIR}

# Dependencies
${OBJ_DIR}/dio.o : $(filter-out ${OBJ_DIR}/dio.o, ${OBJECTS})

${OBJ_DIR}/${SRC_FILE_PREFIX}_globals.o : $(OBJ_DIR)/${SRC_FILE_PREFIX}_constants.o

${OBJ_DIR}/${SRC_FILE_PREFIX}_mathmethods.o : ${OBJ_DIR}/${SRC_FILE_PREFIX}_constants.o ${OBJ_DIR}/${SRC_FILE_PREFIX}_globals.o

${OBJ_DIR}/${SRC_FILE_PREFIX}_inout.o : ${OBJ_DIR}/${SRC_FILE_PREFIX}_constants.o ${OBJ_DIR}/${SRC_FILE_PREFIX}_globals.o \
										${OBJ_DIR}/${SRC_FILE_PREFIX}_expectation.o

${OBJ_DIR}/${SRC_FILE_PREFIX}_field.o : ${OBJ_DIR}/${SRC_FILE_PREFIX}_constants.o ${OBJ_DIR}/${SRC_FILE_PREFIX}_globals.o \
						   			    ${OBJ_DIR}/${SRC_FILE_PREFIX}_inout.o ${OBJ_DIR}/${SRC_FILE_PREFIX}_dirac_equation.o

${OBJ_DIR}/${SRC_FILE_PREFIX}_force.o : ${OBJ_DIR}/${SRC_FILE_PREFIX}_constants.o ${OBJ_DIR}/${SRC_FILE_PREFIX}_globals.o \
						   			    ${OBJ_DIR}/${SRC_FILE_PREFIX}_field.o

${OBJ_DIR}/${SRC_FILE_PREFIX}_nucleus.o : ${OBJ_DIR}/${SRC_FILE_PREFIX}_constants.o ${OBJ_DIR}/${SRC_FILE_PREFIX}_globals.o

${OBJ_DIR}/${SRC_FILE_PREFIX}_basis.o : ${OBJ_DIR}/${SRC_FILE_PREFIX}_constants.o ${OBJ_DIR}/${SRC_FILE_PREFIX}_globals.o \
						   			    ${OBJ_DIR}/${SRC_FILE_PREFIX}_mathmethods.o

${OBJ_DIR}/${SRC_FILE_PREFIX}_constraint.o : ${OBJ_DIR}/${SRC_FILE_PREFIX}_constants.o ${OBJ_DIR}/${SRC_FILE_PREFIX}_globals.o

${OBJ_DIR}/${SRC_FILE_PREFIX}_dirac_BCS.o : ${OBJ_DIR}/${SRC_FILE_PREFIX}_constants.o ${OBJ_DIR}/${SRC_FILE_PREFIX}_globals.o \
											${OBJ_DIR}/${SRC_FILE_PREFIX}_dirac_equation.o ${OBJ_DIR}/${SRC_FILE_PREFIX}_mathmethods.o

${OBJ_DIR}/${SRC_FILE_PREFIX}_dirac_equation.o : ${OBJ_DIR}/${SRC_FILE_PREFIX}_constants.o ${OBJ_DIR}/${SRC_FILE_PREFIX}_globals.o \
												 ${OBJ_DIR}/${SRC_FILE_PREFIX}_mathmethods.o

${OBJ_DIR}/${SRC_FILE_PREFIX}_density.o : ${OBJ_DIR}/${SRC_FILE_PREFIX}_constants.o ${OBJ_DIR}/${SRC_FILE_PREFIX}_globals.o

${OBJ_DIR}/${SRC_FILE_PREFIX}_expectation.o : ${OBJ_DIR}/${SRC_FILE_PREFIX}_constants.o ${OBJ_DIR}/${SRC_FILE_PREFIX}_globals.o

${OBJ_DIR}/${SRC_FILE_PREFIX}_broyden.o : ${OBJ_DIR}/${SRC_FILE_PREFIX}_constants.o ${OBJ_DIR}/${SRC_FILE_PREFIX}_globals.o \
										  ${OBJ_DIR}/${SRC_FILE_PREFIX}_dirac_equation.o ${OBJ_DIR}/${SRC_FILE_PREFIX}_RHB_delta_field.o \
										  ${OBJ_DIR}/${SRC_FILE_PREFIX}_mathmethods.o ${OBJ_DIR}/${SRC_FILE_PREFIX}_field.o

${OBJ_DIR}/${SRC_FILE_PREFIX}_RHB_delta_field.o : ${OBJ_DIR}/${SRC_FILE_PREFIX}_constants.o ${OBJ_DIR}/${SRC_FILE_PREFIX}_globals.o

${OBJ_DIR}/${SRC_FILE_PREFIX}_RHB_equation.o : ${OBJ_DIR}/${SRC_FILE_PREFIX}_constants.o ${OBJ_DIR}/${SRC_FILE_PREFIX}_globals.o \
											   ${OBJ_DIR}/${SRC_FILE_PREFIX}_mathmethods.o

${OBJ_DIR}/${SRC_FILE_PREFIX}_expectation_rotation.o : ${OBJ_DIR}/${SRC_FILE_PREFIX}_constants.o ${OBJ_DIR}/${SRC_FILE_PREFIX}_globals.o \

debug:
	@echo "src path: ${SRC_DIR}/"
	@echo "mod path: ${MOD_DIR}/"
	@echo "obj path: ${OBJ_DIR}/"
	@echo "exe path: ${EXE_DIR}/${EXE_NAME}"
	@echo "SOURCES : ${SOURCES}"
	@echo "OBJECTS : ${OBJECTS}"

clean:
	rm -f ${EXE_DIR}/${EXE_NAME} $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod

deepclean:
	rm -rf ${EXE_DIR} $(OBJ_DIR) $(MOD_DIR) ${SRC_DIR}/*.mod


CC		= g++
SRC_DIR	= src
OBJ_DIR	= obj
SRC_EXT	= .cpp
SRC		= $(SRC_DIR)/MainFilterKruskal.cpp
OBJ		= $(OBJ_DIR)/MainFilterKruskal.o
EXEC	= FilterKruskal
DEBUG	= -DNDEBUG
OPT		= -O3
WARN	= -Wall
MISC	= -c -ffast-math
CFLAG	= $(OPT) $(DEBUG) $(WARN) $(MISC)

LIB		= 
INC		= 

RM		= rm -rf

MKDIR	= mkdir -p

all: dir $(EXEC)
	
dir: 
	$(MKDIR) $(OBJ_DIR)
	
$(EXEC): $(OBJ)
	@echo "Building Target: $@"
	$(CC) -o "$@" $(OBJ) $(LIB)
	@echo "Finished Building Target: $@"
	
$(OBJ_DIR)/%.o: $(SRC_DIR)/%$(SRC_EXT)
	@echo "Building File: $<"
	$(CC) $(INC) $(CFLAG) -o "$@" "$<"
	@echo "Finished Building File: $<"
	
clean:
	$(RM) $(OBJ) $(EXEC)

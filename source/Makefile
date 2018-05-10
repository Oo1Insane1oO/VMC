INC = -I/usr/local/include -I/usr/local/include/eigen3
LIBS = -L/usr/local/lib

TESTS ?= OFF
ifeq ($(TESTS),ON)
	INC += -I/usr/local/include/UnitTest++ 
	LIBS += -L/usr/local/lib/libUnitTest++.a -lUnitTest++
endif

CC = mpic++
LDFLAGS	= -fopenmp -lpthread -lyaml-cpp $(LIBS)
CXXFLAGS = -O3 -fopenmp -std=c++17 $(INC)

MPI ?= ON
ifeq ($(MPI),ON)
	CXXFLAGS += -DMPI_ON
endif

DEBUG ?= OFF
ifeq ($(DEBUG),ON)
	CXXFLAGS := -g -std=c++17 -Wall -Wextra $(INC)
endif

ifeq ($(TESTS),ON) 
	CXXFLAGS += -DTESTS
endif

JASTROW = PADEJASTROW
WAVEFUNCTION = HARMONICOSCILLATOR 
CXXFLAGS += -D$(WAVEFUNCTION) -D$(JASTROW)

EXEC = main
SRC = src

SUFFIXES += .d

INPUTS_DIR = inputs
DATA_DIR = data

JAST_DIR = $(SRC)/jastrow
JAST_LOWER = $(shell echo $(JASTROW) | tr A-Z a-z)
JAST = $(JAST_DIR)/$(JAST_LOWER)

WAVE_DIR = $(SRC)/wavefunctions
WAVE_LOWER = $(shell echo $(WAVEFUNCTION) | tr A-Z a-z)
WAVE = $(WAVE_DIR)/$(WAVE_LOWER)
OBJ_DIR = bin
TEST_DIR = tests

BASIS = $(SRC)/basis/$(WAVE_LOWER)basis

OTHER_BASIS_SOURCE = $(SRC)/basisfunctions/cartesian.cpp

MAINSRC = $(SRC)/main.cpp

ifeq ($(TESTS), ON)
	TEST_SRC = $(wildcard $(TEST_DIR)/*.cpp)
	TEST_OBJ += $(patsubst $(TEST_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(TEST_SRC))
	TEST_DEPS = $(patsubst %.o,%.d,$(TEST_OBJ))
endif

TMPSRC = $(wildcard $(SRC)/*.cpp) $(WAVE).cpp $(BASIS).cpp $(OTHER_BASIS_SOURCE) $(wildcard $(SRC)/examples/*.cpp) $(JAST).cpp
SOURCES = $(filter-out $(MAINSRC), $(TMPSRC))
HEADERS = $(wildcard $(SRC)/*.h) $(WAVE).h $(BASIS).h $($(OTHER_BASIS_SOURCE):.cpp=.h) $(JAST).h
ifeq ($(JASTROW),PADENQS)
	SOURCES += $(JAST_DIR)/padejastrow.cpp $(JAST_DIR)/rbmjastrow.cpp
	HEADERS += $(JAST_DIR)/padejastrow.h $(JAST_DIR)/rbmjastrow.h
endif
ifeq ($(JASTROW),EXPNQS)
	SOURCES += $(JAST_DIR)/padejastrow.cpp $(JAST_DIR)/rbmjastrow.cpp $(JAST_DIR)/expnqs.cpp
	HEADERS += $(JAST_DIR)/padejastrow.h $(JAST_DIR)/rbmjastrow.h $(JAST_DIR)/expnqs.h
endif

SOURCES += $(MAINSRC)
OBJECTS = $(patsubst $(SRC)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))

DEPS = $(patsubst %.o,%.d,$(OBJECTS))

build: $(EXEC)

clean:
	rm -rf $(OBJ_DIR)
	rm $(EXEC)

createDataDir:
	mkdir $(DATA_DIR)

cleanData:
	rm -rf $(DATA_DIR)

cleanInputs:
	rm $(wildcard $(INPUTS_DIR)/*)

.PHONY: clean build

NODEPS := clean tags svn
ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
    -include $(DEPS) $(TEST_DEPS)
endif

$(EXEC): $(OBJECTS) | $(TEST_OBJ)
	$(CC) $^ $(LDFLAGS) -o $@

$(OBJ_DIR)/%.d: $(SRC)/%.cpp
	mkdir -p $(@D)
	$(CC) $(CXXFLAGS) -MM -MT '$(patsubst $(SRC)/%.cpp,$(OBJ_DIR)/%.o,$<)' $< -MF $@

$(OBJ_DIR)/%.d: $(TEST_DIR)/%.cpp
	mkdir -p $(@D)
	$(CC) $(CXXFLAGS) -MM -MT '$(patsubst $(TEST_DIR)/%.cpp,$(OBJ_DIR)/%.o,$<)' $< -MF $@

$(OBJ_DIR)/%.o: $(SRC)/%.cpp $(OBJ_DIR)/%.d $(HEADERS)
	mkdir -p $(@D)
	$(CC) $(CXXFLAGS) -o $@ -c $<

$(OBJ_DIR)/%.o: $(TEST_DIR)/%.cpp $(OBJ_DIR)/%.d
	mkdir -p $(@D)
	$(CC) $(CXXFLAGS) -o $@ -c $<

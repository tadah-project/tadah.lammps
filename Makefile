###############################################################################
# To compile Tadah.MD just type 'make' inside the lib root directory.
###############################################################################
CC:= g++
CXXFLAGS:= -std=c++11 -g -O3 -Wall -fPIC -Wpedantic -pedantic-errors -Wextra
###############################################################################

ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
BUILD_DIR:=$(ROOT_DIR)/build

CORE:=$(ROOT_DIR)/TADAH.MD/CORE
MODELS:=$(ROOT_DIR)/TADAH.MD/MODELS
TADAHMD:=$(ROOT_DIR)/TADAH.MD

CORE_SRC := $(shell find $(CORE) -name "*.cpp")
MODELS_SRC := $(shell find $(MODELS) -name "*.cpp")
TADAHMD_SRC := $(shell find $(TADAHMD) -maxdepth 1 -name "*.cpp")

CORE_OBJ := $(CORE_SRC:%.cpp=$(BUILD_DIR)/%.o)
MODELS_OBJ := $(MODELS_SRC:%.cpp=$(BUILD_DIR)/%.o)
TADAHMD_OBJ := $(TADAHMD_SRC:%.cpp=$(BUILD_DIR)/%.o)

TADAH_CONFIG_DIR:=$(ROOT_DIR)/TADAH.MD/CORE/config/
tadah_lnk:=-DTADAH_CONFIG_DIR=\"${TADAH_CONFIG_DIR}\"

LIB.so:=libtadah.md.so
LIB.a:=libtadah.md.a

dynamic: $(CORE_OBJ) $(MODELS_OBJ) $(TADAHMD_OBJ) $(LIB.so) move remove_dirs
static: $(CORE_OBJ) $(MODELS_OBJ) $(TADAHMD_OBJ) $(LIB.a) move remove_dirs

obj: $(CORE_OBJ) $(MODELS_OBJ) $(TADAHMD_OBJ) move remove_dirs

.DEFAULT_GOAL := static

$(CORE_OBJ): $(BUILD_DIR)/%.o : %.cpp 
	mkdir -p $(dir $@)
	$(CC) $(CXXFLAGS) $(tadah_lnk) -c $< -o $@

$(MODELS_OBJ): $(BUILD_DIR)/%.o : %.cpp
	mkdir -p $(dir $@)
	$(CC) $(CXXFLAGS) $(tadah_lnk) -c $< -o $@

$(TADAHMD_OBJ): $(BUILD_DIR)/%.o : %.cpp
	mkdir -p $(dir $@)
	$(CC) $(CXXFLAGS) $(tadah_lnk) -c $< -o $@

move: $(CORE_OBJ) $(MODELS_OBJ) $(TADAHMD_OBJ)
	$(foreach f, $(shell find $(BUILD_DIR) -mindepth 2 -name "*.o"), $(shell mv $(f) $(BUILD_DIR)))

remove_dirs: $(CORE_OBJ) $(MODELS_OBJ) $(TADAHMD_OBJ) move
	$(shell find $(BUILD_DIR) -type d -empty -delete)

$(LIB.so): $(CORE_OBJ) $(MODELS_OBJ) $(TADAHMD_OBJ) move
	@echo Building $(LIB.so)...
	$(CC) $(CXXFLAGS) -shared -o $(BUILD_DIR)/$@ $(BUILD_DIR)/*.o

$(LIB.a): $(CORE_OBJ) $(MODELS_OBJ) $(TADAHMD_OBJ) move
	@echo Building $(LIB.a)...
	ar rcs $(BUILD_DIR)/$@ $(BUILD_DIR)/*.o

clean:
	rm -f $(BUILD_DIR)/*.o
	rm -f $(BUILD_DIR)/$(LIB.so)
	rm -f $(BUILD_DIR)/$(LIB.a)

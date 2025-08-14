TARGET_NAME := chemeq
OPT=-O3
LDFLAGS=-lm
INCDIRS=./include

# tool macros
CC ?= gcc
CXX ?= # FILL: the compiler
CFLAGS := -Wall -Wextra $(foreach D,$(INCDIRS),-I$(D)) $(OPT) $(DEPFLAGS)
CXXFLAGS := # FILL: compile flags
DBGFLAGS := -g
COBJFLAGS := $(CFLAGS) -c

# path macros
BIN_PATH := bin
OBJ_PATH := obj
SRC_PATH := src
SHR_PATH := $(SRC_PATH)/shared
TRG_PATH := $(SRC_PATH)/targets
INC_PATH := include
DBG_PATH := debug

# compile macros
ifeq ($(OS),Windows_NT)
	TARGET_NAME := $(addsuffix .exe,$(TARGET_NAME))
endif
TARGET := $(BIN_PATH)/$(TARGET_NAME)
TARGET_DEBUG := $(DBG_PATH)/$(TARGET_NAME)

# src files & obj files
SRC_SHARED := $(foreach x, $(SHR_PATH), $(wildcard $(addprefix $(x)/*,.c)))
#SRC_TARGET := $(foreach x, $(TRG_PATH), $(wildcard $(addprefix $(x)/*,.c*)))
SRC_TARGET := $(TRG_PATH)/$(TARGET_NAME).c
SRC_ALL := $(SRC_TARGET) $(SRC_SHARED)
OBJ := $(addprefix $(OBJ_PATH)/, $(addsuffix .o, $(notdir $(basename $(SRC_ALL)))))
OBJ_DEBUG := $(addprefix $(DBG_PATH)/, $(addsuffix .o, $(notdir $(basename $(SRC_ALL)))))

# clean files list
DISTCLEAN_LIST := $(OBJ) \
                  $(OBJ_DEBUG)
CLEAN_LIST := $(TARGET) \
			  $(TARGET_DEBUG) \
			  $(DISTCLEAN_LIST)

# default rule
default: makedir all

# non-phony targets
$(TARGET): $(OBJ)
	$(CC) -o $@ $(OBJ) $(CFLAGS) $(LDFLAGS)

$(OBJ_PATH)/%.o: $(SHR_PATH)/%.c*
	$(CC) $(COBJFLAGS) -o $@ $<

$(OBJ_PATH)/$(TARGET_NAME).o: $(TRG_PATH)/$(TARGET_NAME).c*
	$(CC) $(COBJFLAGS) -o $@ $<

$(DBG_PATH)/%.o: $(TRG_PATH)/$(TARGET_NAME).c*
	$(CC) $(COBJFLAGS) $(DBGFLAGS) -o $@ $<

$(TARGET_DEBUG): $(OBJ_DEBUG)
	$(CC) $(CFLAGS) $(DBGFLAGS) $(OBJ_DEBUG) -o $@

# phony rules
.PHONY: makedir
makedir:
	@mkdir -p $(BIN_PATH) $(OBJ_PATH) $(DBG_PATH) $(INC_PATH)

.PHONY: all
all: $(TARGET)

.PHONY: debug
debug: $(TARGET_DEBUG)

.PHONY: clean
clean:
	@echo CLEAN $(CLEAN_LIST)
	@rm -f $(CLEAN_LIST)

.PHONY: cleanall
cleanall:
	rm --verbose ./bin/* ./obj/* ./debug/*

.PHONY: distclean
distclean:
	@echo CLEAN $(DISTCLEAN_LIST)
	@rm -f $(DISTCLEAN_LIST)

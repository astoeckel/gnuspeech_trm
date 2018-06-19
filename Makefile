# based on https://spin.atomicobject.com/2016/08/26/makefile-c-projects/

TARGET_LIB ?= libgnuspeech_trm.so
TARGET_PROG ?= gnuspeech_trm_main

BUILD_DIR ?= ./build
SRC_DIRS ?= ./gnuspeech_trm

SRCS := $(shell find $(SRC_DIRS) -name '*.cpp')
SRCS_LIB := $(shell find $(SRC_DIRS) -name '*.cpp' -and -not -name '*_main.cpp')

OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
OBJS_LIB := $(SRCS_LIB:%=$(BUILD_DIR)/%.o)

DEPS := $(OBJS:.o=.d)

INC_DIRS := $(shell find $(SRC_DIRS) -type d)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS ?= $(INC_FLAGS) -std=c++11 -g -O0 -Wall -Wextra -pedantic -fPIC -MMD -MP -fvisibility=hidden -fvisibility-inlines-hidden

all: $(BUILD_DIR)/$(TARGET_LIB) $(BUILD_DIR)/$(TARGET_PROG)

$(BUILD_DIR)/$(TARGET_LIB): $(OBJS_LIB)
	$(CXX) -fvisibility=hidden -fvisibility-inlines-hidden -fPIC -shared $(OBJS_LIB) -o $@ -lm $(LDFLAGS)

$(BUILD_DIR)/$(TARGET_PROG): $(OBJS)
	$(CXX) $(OBJS) -o $@ -lm $(LDFLAGS)

# c++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

.PHONY: clean

clean:
	$(RM) -r $(BUILD_DIR)

-include $(DEPS)

MKDIR_P ?= mkdir -p


SRC_DIR = ./src
GCC_OPTS = -O3 -std=c99

all: colortobase

colortobase: $(SRC_DIR)/colortobase.c
	gcc $(GCC_OPTS) $< -o $@

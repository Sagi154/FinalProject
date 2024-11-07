CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
SRC_FILE = symnmf.c

all: symnmf

symnmf: symnmf.o
	$(CC) $(CFLAGS) symnmf.o -o symnmf -lm

symnmf.o: $(SRC_FILE)
	$(CC) -c $(SRC_FILE) $(CFLAGS)

clean:
	rm -f symnmf symnmf.o
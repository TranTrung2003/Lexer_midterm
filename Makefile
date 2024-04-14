CC ?= gcc

scan: VC.o Lexer.o
	$(CC) $^ -o $@

.PHONY: clean
clean:
	rm -rf VC.o Lexer.o scan



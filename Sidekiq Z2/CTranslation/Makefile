COMP = gcc
EXEC = main_executable
RECV = receive
TRANS = transmit

$(EXEC) :
	$(COMP) main.c -o $(EXEC) -lm -lfftw3
$(RECV) :
	$(COMP) receive.c -o $(RECV) -lm -lfftw3
$(TRANS) :
	$(COMP) transmit.c -o $(TRANS) -lm -lfftw3
clean :
	rm $(EXEC)
	rm $(RECV)
	rm $(TRANS)
all :
	make $(TRANS)
	make $(RECV)
	make $(EXEC)
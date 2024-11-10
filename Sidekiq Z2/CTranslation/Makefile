COMP = gcc
EXEC = main_executable
RECV = receive
TRANS = transmit
INCLUDE =
LINK =
LIBS = -lm

$(EXEC) :
	$(COMP) main.c -o $(EXEC) -I $(INCLUDE) -L $(LINK) $(LIBS)
$(RECV) :
	$(COMP) receive.c -o $(RECV) -I $(INCLUDE) -L $(LINK) $(LIBS)
$(TRANS) :
	$(COMP) transmit.c -o $(TRANS) -I $(INCLUDE) -L $(LINK) $(LIBS)
clean :
	rm $(EXEC)
	rm $(RECV)
	rm $(TRANS)
all :
	make $(TRANS)
	make $(RECV)
	make $(EXEC)
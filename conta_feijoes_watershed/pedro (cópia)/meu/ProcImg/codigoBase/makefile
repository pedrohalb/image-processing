CC=gcc
#CCFLAGS=-Wall
CCFLAGS=-g
LDFLAGS=-lm 

# Fonte e objetos
SOURCES=$(wildcard *.c)
OBJECTS=$(SOURCES:.c=.o)
TARGET=contafeijao

# Verifica o sistema operacional para definir o caminho correto da biblioteca
ifeq ($(OS), Windows_NT)
   LIB_SRC = ..\\utils\\imagelib.c
   LIB_OBJ = ..\\utils\\imagelib.o
else
   LIB_SRC = ../utils/imagelib.c
   LIB_OBJ = ../utils/imagelib.o
endif

# Regra principal
all: $(TARGET)

# Regra para compilar o executável
$(TARGET): $(OBJECTS) $(LIB_OBJ)
	$(CC) -g -o $@ $^ $(LDFLAGS)

# Regra para compilar os objetos
%.o: %.c
	$(CC) $(CCFLAGS) -c $<

# Regra para compilar o objeto da biblioteca
$(LIB_OBJ): $(LIB_SRC)
	$(CC) $(CCFLAGS) -c $< -o $@

# Regra para limpar os arquivos de objeto e o executável
clean:
	rm -f *.o *-result.p*m $(TARGET) $(LIB_OBJ)
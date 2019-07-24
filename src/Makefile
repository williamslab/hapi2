CPPSRCS= main.cc cmdlineopts.cc phaser.cc analysis.cc swapphase.cc
CSRCS= memcpy.c
OBJS= $(patsubst %.cc,%.o,$(CPPSRCS)) $(patsubst %.c,%.o,$(CSRCS))
EXEC= hapi2

GPP = g++
GCC = gcc
DEFINES= 
CFLAGS = -Wall $(DEFINES)
CPPFLAGS = -std=c++11 $(CFLAGS)
ifdef DEBUG           # to use run `make DEBUG=1`
  CFLAGS += -g
else
  CFLAGS += -O2
endif

# profiling:
ifdef PROFILE       # to use run `make PROFILE=1
  CFLAGS += -pg
endif

LIBS = -static-libstdc++ -static-libgcc -lgenetio -lz -lhts -lpthread -Wl,--wrap=memcpy

# dependency variables / commands
DEPDIR = .deps
df = $(DEPDIR)/$(*F)

all: $(EXEC)

$(EXEC): $(OBJS) $(HEADERS)
	$(GPP) -o $(EXEC) $(OBJS) $(CFLAGS) $(LIBS)

# This way of building dependencies (per-file) described at
# http://make.paulandlesley.org/autodep.html

.c.o:
	@mkdir -p $(DEPDIR)
	$(GCC) -MMD $(CFLAGS) -o $@ -c $<
	@cp $*.d $(df).P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
	  rm -f $*.d

.cc.o:
	@mkdir -p $(DEPDIR)
	$(GPP) -MMD $(CPPFLAGS) -o $@ -c $<
	@cp $*.d $(df).P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
	  rm -f $*.d


# include the .P dependency files, but don't warn if they don't exist (the -)
-include $(CPPSRCS:%.cc=$(DEPDIR)/%.P)
-include $(CSRCS:%.c=$(DEPDIR)/%.P)
# The following applies if we don't use a dependency directory:
#-include $(SRCS:.cc=.P)

tags: $(CPPSRCS) $(CSRCS) *.h
	ctags --language-force=c++ --extra=+q --fields=+i --excmd=n *.cc *.h

clean:
	rm -f $(EXEC) $(OBJS)

clean-deps:
	rm -f $(DEPDIR)/*.P

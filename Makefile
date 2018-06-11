include version.mak
VERSION=${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_SUBMINOR}

VERSION_FLAGS= \
	-DVERSION_MAJOR="\"${VERSION_MAJOR}"\" \
	-DVERSION_MINOR="\"${VERSION_MINOR}"\" \
	-DVERSION_SUBMINOR="\"${VERSION_SUBMINOR}"\" \
	-DREVISION_DATE="\"${REVISION_DATE}"\" \
	-DSUBVERSION_REV="\"${SUBVERSION_REV}"\"

CFLAGS = -O3 ${VERSION_FLAGS} -Wall -Wextra -Werror
CC     = gcc
LDLIBS = -lm

incFiles = loop_aligner.h feed.h motifs.h seq_ops.h clock.h utilities.h \
           version.mak

default: noise_cancelling_repeat_finder
	mv noise_cancelling_repeat_finder NCRF

noise_cancelling_repeat_finder: noise_cancelling_repeat_finder.o loop_aligner.o feed.o motifs.o seq_ops.o clock.o utilities.o

%.o: %.c Makefile ${incFiles}
	${CC} -c ${CFLAGS} $< -o $@

clean: cleano
	rm -f NCRF noise_cancelling_repeat_finder

cleano:
	rm -f *.o

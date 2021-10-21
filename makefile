CC     = g++
CFLAGS = -Wall -Wextra -ggdb -I ../../../Klassen/ -I .

OliInSeq-v0.9.6: global-types-and-parameters.cpp remove_outlier_sequences.cpp global-types-and-parameters.h 
	$(CC) $(CFLAGS) -o OliInSeq-v0.9.6 global-types-and-parameters.cpp remove_outlier_sequences.cpp

clean:
	rm OliInSeq-v0.9.6




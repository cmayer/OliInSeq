CC     = g++
# CFLAGS = -Wall -Wextra -ggdb  -I .
CFLAGS = -Wall -Wextra -O2  -I .

OliInSeq-v0.9.6: global-types-and-parameters.cpp remove_outlier_sequences.cpp global-types-and-parameters.h 
	$(CC) $(CFLAGS) -o OliInSeq-v0.9.6 global-types-and-parameters.cpp remove_outlier_sequences.cpp

clean:
	rm OliInSeq-v0.9.6




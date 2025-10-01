# Change list
A list of changes made for the purpose of performance

1. Signature density calcualted at start instead of every new signature
2. Instead of that weird random allocation *compute_new_term_sig* uses memset and std::random_shuffle
3. Uses simd to add 32 int8_t's at a time when adding up the signature
4. straight up removed find_sig/hash_table though it might be worth adding back for larger files

# Potential future changes
(Find_sig needs further investigations)

1. partition can handle threads. it can be the sole file reader giving data and doc number to threads. doc number indicates where to write in a file so multiple threads can write.

- The file is currently read a line at a time so a meta line is read each time.
- multiple threads maybe able to read/write same file at same time. look into mmap (virtually mapping file into memory to access as pointers OR simply open multiple file descriptors)
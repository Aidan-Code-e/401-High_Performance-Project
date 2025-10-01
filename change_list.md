# Change list
A list of changes made for the purpose of performance

1. Signature density calcualted at start instead of every new signature
2. Instead of that weird random allocation *compute_new_term_sig* uses memset and std::random_shuffle
3. 

# Potential future changes
(Find_sig needs further investigations)

1. The loop in signature_add is NEEDS to be addressed

- The file is currently read a line at a time so a meta line is read each time.
- multiple threads maybe able to read/write same file at same time. look into mmap (virtually mapping file into memory to access as pointers OR simply open multiple file descriptors)
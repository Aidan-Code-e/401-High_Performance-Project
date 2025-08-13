# Change list
A list of changes made for the purpose of performance

1. Signature density calcualted at start instead of every new signature
2. Instead of that weird random allocation *compute_new_term_sig* uses memset and std::random_shuffle
3. 
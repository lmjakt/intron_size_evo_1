#include <R.h>
#include <Rinternals.h>
#include <string.h>

// Split a string strictly on a sincle split character.
// This differs from strsplit() in that it will return an empty
// string for a trailing split character

struct splitted {
  char **sub_strings;
  int n;
  int capacity;
};

void grow_capacity(struct splitted *sp){
  sp->capacity *= 2;
  char **buffer = malloc( sizeof(char*) * sp->capacity );
  memcpy( buffer, sp->sub_strings, sizeof(char*) * sp->n );
  free( sp->sub_strings );
  sp->sub_strings = buffer;
}

void push_splitted(struct splitted *sp, char *w){
  if(sp->n >= sp->capacity)
    grow_capacity(sp);
  sp->sub_strings[sp->n] = w;
  sp->n++;
}

struct splitted init_splitted(int capacity){
  struct splitted sp;
  sp.capacity = capacity;
  sp.sub_strings = malloc(sizeof(char*) * sp.capacity);
  sp.n = 0;
  return(sp);
}

void clear_splitted( struct splitted *sp ){
  sp->n = 0;
}

void split_string( char *string, char split_c, struct splitted* sp ){
  int i = 0;
  clear_splitted(sp);
  push_splitted( sp, string );
  while( string[i] != 0 ){
    if( string[i] == split_c ){
      string[i] = 0;
      push_splitted(sp, string + i + 1 );
    }
    i++;
  }
}

SEXP strict_split(SEXP strings, SEXP split){
  if(TYPEOF(strings) != STRSXP || length(strings) < 1)
    error("strings must be a character vector and have positive length");
  if(TYPEOF(split) != STRSXP || length(split) < 1)
      error("split must be a character vector and have positive length");

  char split_c = CHAR(STRING_ELT(split, 0))[0];
  if(!split_c)
    error("The split string must contain at least one character (the first of which will be used)");
  
  int n = length(strings);
  SEXP ret_data = PROTECT(allocVector( VECSXP, n ));
  // we reuse the same struct for each of the strings so that we do not need to malloc
  // and free stuff as we g along;
  struct splitted sp = init_splitted(10);

  // go through each string;
  for(int i=0; i < n; i++){
    const char *string = CHAR( STRING_ELT(strings, i) );
    // we need to copy this as the function modifies it putting 0s in
    size_t string_l = strlen( string );
    char *string_cpy = malloc(sizeof(char) * (string_l + 1));
    strncpy( string_cpy, string, string_l + 1 ); // actually no different to memcpy
    split_string( string_cpy, split_c, &sp );

    SET_VECTOR_ELT( ret_data, i, allocVector(STRSXP, sp.n) );
    for(int j=0; j < sp.n; ++j)
      SET_STRING_ELT( VECTOR_ELT( ret_data, i ), j, mkChar( sp.sub_strings[j] ));

    // clear the relevant data structures
    clear_splitted(&sp);
    free( string_cpy );
  }
  free( sp.sub_strings );

  UNPROTECT(1);
  return(ret_data);
}
  

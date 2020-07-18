#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <time.h>

#define BLOCKSIZE 10000000
#define TERMINATOR '@'

typedef struct {
  char   *sequence;
  char   *sequence2;  // dual indexed reads
  char   *sequenceRev; // paired reads
  int    original_pos;
} a_barcode;

typedef struct {
   char   *sequence;
   int    original_pos;
} a_hairpin;

typedef struct {
  long sequence_index;
} end_node;

typedef struct trie_node trie_node;
struct trie_node {
    char    base;
    long    count;
    // links is [@, A, C, G, T]
    trie_node *links[5];
    end_node *end;
};


a_barcode **barcodes;
a_hairpin **hairpins;
trie_node *hairpin_trie_head;
trie_node *barcode_single_trie_head;
trie_node *barcode_paired_trie_head;
trie_node *barcode_dualindex_trie_head;

long *barcode_positions;
int barcode_positions_size;
long *hairpin_positions;
int hairpin_positions_size;

int is_PairedReads; 
int is_DualIndexingReads;
int barcodesInHeader;
int num_barcode;
int num_hairpin;
long num_read;
long **summary;
int barcode_length;
int barcode2_length;
int barcode_length_rev;
int hairpin_length;

/*
int barcode_start;
int barcode_end;
int barcode2_start;
int barcode2_end;
int barcode_start_rev;
int barcode_end_rev;
int hairpin_start;
int hairpin_end;
int allow_shifting;
int shifting_n_base; 
int allow_shifted_mismatch;
*/

int allow_mismatch;
int barcode_n_mismatch;
int hairpin_n_mismatch;
int isverbose;
int barcodes_in_header;

long barcodecount;
long hairpincount;
long bchpcount;

long longest_read_length;

// Management functions for building and traversing a Trie.
int
Get_Links_Position(char base) {
  /*
  Determine the array position of the given base. 
  0, 1, 2, 3, 4 are @ A C G T respectively.
  base: the char to convert to int. Expects either @, A, C, G or T
  return: an int in the inclusive range 0 - 4
  */
  switch (base) {
      case 'A':
        return 1;
        break;
      case 'C':
        return 2;
        break;
      case 'G':
        return 3;
      case 'T':
        return 4;
      case TERMINATOR:
      default:
        return 0;
  }
}

trie_node*
Initialise_Node(char base){
  /*
  Initialise a trie node, which is a struct containing the base, the insertion count, 
  and an array to link the attached trie_nodes.
  base: the base rep to create a node of
  return: a pointer to the created node
  */
  trie_node* this_node = (trie_node *)malloc(sizeof(trie_node));
  this_node->base = base;
  this_node->count = 0;
  this_node->end = NULL;
  
  int i;
  for (i = 0; i < 5; i++) {
      this_node->links[i] = NULL;
  }
  
  return this_node;
}

trie_node*
Initialise_End_Node(char base, long sequence_index) {
  /* 
  The end node is a special trie node, containing the index of the sequence 
  which ends at this node in the corresponding barcode or hairpin array
  base: the base of this new node, which will always be the terminator character.
  sequence_index: the index in the associated array of sequence which ends at this node.
  */
  trie_node *this_node = Initialise_Node(base);
  end_node *end = (end_node *)malloc(sizeof(end_node));

  end->sequence_index = sequence_index;
  this_node->end = end;
  return this_node;
}

bool
Base_In_Node(trie_node* node, char base) {
  /*
  Check if the given node contains a link to the given base.
  Checks for the NULL pointer in the position signified by base.
  node: the node to check existance of base in
  base: the base to check for
  return: true if base exists in node->links
          false otherwise.
  */
  if (node->links[Get_Links_Position(base)] != NULL) {
    return true;
  }
  return false;
}

trie_node* 
Add_Node(trie_node *node, char base) {
  /* 
  Adds a trie node to the given node's internal list of nodes,
  returning a pointer to that new node
  Requires the node to not currently contain a link to a node in the array at index base.
  node: a pointer to the node to insert at
  base: the char rep of the node to create
  return: a pointer to the new node created
  */
  node->count++;
  trie_node *new_node = Initialise_Node(base);
  node->links[Get_Links_Position(base)] = new_node;
  return new_node;
}

trie_node*
Add_End_Node(trie_node *node, char base, long sequence_index) {
  /* 
  Adds an end node to the trie, which contains all of the 
  same data as a regular node, but with an additional array
  to store the index of the completed hairpin sequences
  node: a pointer to the node to insert at
  base: the char rep of the node to create
  sequence_index: the index of the sequence ending at this node
  return: a pointer to the new node created
  */
  node->count++;
  trie_node *new_node = Initialise_End_Node(base, sequence_index);
  node->links[Get_Links_Position(base)] = new_node;
  return new_node;
}

trie_node*
Build_Trie_Hairpins(void) {
  /*
  Build a trie using the global hairpins array.
  For every hairpin in the hairpins array, add it to the trie by:
    Starting at the head node, add a node for the current char if none exists, or follow to the relevant node
    Repeat for every char in the hairpin, finally inserting a terminator character at the end of the trie path
    Hairpins are assumed to be unique, which is checked in R wrapper function
  
  return: a pointer to the head of the trie, which contains the empty character
  */
  trie_node *head = Initialise_Node('\0');
  char *cur_seq;
  int hp_i, insert_i;
  char insert_base;
  trie_node *current_node;
  // for every hairpin in hairpins, insert it into the trie
  // This is done by following the trie until we reach a base not already in the trie,
  // inserting it and then continuing this process until all bases are linked together
  for (hp_i = 1; hp_i <= num_hairpin; hp_i++) {
    current_node = head;
    cur_seq = hairpins[hp_i]->sequence;

    for (insert_i = 0; insert_i < hairpin_length; insert_i++) {
      insert_base = cur_seq[insert_i];

      if (Base_In_Node(current_node, insert_base)) {
        // if the base is in the current node, simple increment the count
        // and move onto the linked node
        current_node->count++;
        current_node = current_node->links[Get_Links_Position(insert_base)];
      } else {
        current_node = Add_Node(current_node, insert_base);
      }
    }
    // insert the final TERMINATOR @
    current_node = Add_End_Node(current_node, TERMINATOR, hairpins[hp_i]->original_pos);

    // increment the last node in the sequence's count, before we insert the next string
    current_node->count++;
  }

  return head;
}

trie_node*
Build_Trie_Barcodes(bool is_paired, bool is_dualindex) {
  /*
  Build a trie using the barcodes array, with parameters allowing support for building paired read, or dual indexing barcodes.
  For every barcode in the barcodes array, add it to the trie by:
    Starting at the head node, add a node for the current char if none exists, or follow to the relevant node
    Repeat for every char in the barcode, finally inserting a terminator character at the end of the trie path
  
  is_paired: boolean value indicating if we should create a paired read trie, storing all sequences in barcodes->sequenceRev
  is_dualindex: boolean value indicating if we should insert barcodes->sequence2 into the created trie.
    Both of these booleans cannot be true. is_paired will be prioritised over is_dualindex
  return: a pointer to the head of the trie, which contains the empty character
  */
  trie_node *head = Initialise_Node('\0');
  trie_node *current_node;
  int length_test, bc_i, insert_i;
  char *cur_seq;
  char insert_base;
  // determine which length we need to use for inserting the barcodes
  if (is_paired) {
    length_test = barcode_length_rev;
  } else if (is_dualindex) {
    length_test = barcode2_length;   
  } else {
     length_test = barcode_length;
  }
  // For every barcode in the barcodes array, add it to the trie by iteratively following the trie
  // through each character in the barcode, adding the characters that don't exist yet.
  // The original position of the barcodes are recorded at the terminator character
  for (bc_i = 1; bc_i <= num_barcode; bc_i++) {
    current_node = head;

    // determine which barcode sequence we will insert
    if (is_paired) {
      cur_seq = barcodes[bc_i]->sequenceRev;
    } else if (is_dualindex) {
      cur_seq = barcodes[bc_i]->sequence2;
    } else {
      cur_seq = barcodes[bc_i] ->sequence;
    }

    // loop through each character in the barcode to insert it into the trie
    for (insert_i = 0; insert_i < length_test; insert_i++) {
      insert_base = cur_seq[insert_i];

      if (Base_In_Node(current_node, insert_base)) {
        // if the base is in the current node, simple increment the count
        // and move onto the linked node
        current_node->count++;
        current_node = current_node->links[Get_Links_Position(insert_base)];
      } else {
        current_node = Add_Node(current_node, insert_base);
      }
    }
    // insert the final TERMINATOR @
    // Barcodes are assumed to be unique, as uniqueness checking is performed in wrapper R functino
    current_node = Add_End_Node(current_node, TERMINATOR, bc_i);
    // increment the last node in the sequence's count, before we insert the next string
    current_node->count++;
  }

  return head;
}

// Management functions for recording barcode and hairpin found positions
long *
Initialise_Resize_Array(int size) {
  /*
  Intialise an array of length size, setting each element as 0.
  Used for a resizable array DS.

  size: the size of the new array
  return: a pointer to the new array.
  */
  long *new = (long *)malloc(sizeof(long) * size);
  
  int i;
  for (i = 0; i < size; i++) {
    new[i] = 0;
  }
  return new;
}

int
Expand_Resize_Array(long **resize_array, int size) {
  /*
  Expands the given array, by instantiating a new array of 
  length size * 2, and copying every element from the old array
  into the new one.

  resize_array: a double pointer to the array to expand. Must be double pointer
    in order to alter the given pointer to the array, so the reassignment of the
    array pointer can be done within this method, as opposed to after this method
    call.
  size: the size of the given array
  return: the size of the new array.
  post-cond: the double pointer given now references the new expanded array,
    and the old array is freed from memory. 
  */
  long *new = Initialise_Resize_Array(size * 2);

  int i;
  for (i = 0; i < size; i++) {
    new[i] = (*resize_array)[i];
  }

  // free the old array, and assign the new array to the given pointer
  free(*resize_array);
  *resize_array = new;
  return size * 2;
}

int
Increment_Resize_Array(long **array, int size, int position) {
  /*
  Increments the count at the given position in the given array.
  If the given position is outside of the bounds of the array,
  the array is expanded and then the position is incremented.

  array: a double pointer to the array to increment
  size: the size of the given array
  position: the index to increment
  return: the size of the new array, if altered, or the old size if not
  post-cond: the double pointer given now references the expanded array,
    with position incremented, or references the old array, with position i
    incremented.
  */
  int current_size = size;
  while (position >= current_size) {
    // This loop is only executed if position is outside the bounds of the array.
    current_size = Expand_Resize_Array(array, current_size);
  }
  (*array)[position]++;
  return current_size;
}

// Functions for reading barcodes and hairpins from input files, and storing in struct arrays
int Get_Lines_In_File(FILE* fin) {
  /*
  Iterates over the file and counts the number of lines in the file.
  Rewinds the file afterwards, for reading later.
  fin: the FILE object to determine the length of
  return: the length of the file.
  */
  int N=0, ch, last_ch='\n';
  while (1) { 
    ch=fgetc(fin);
    if (ch=='\n') { ++N; }
    else if (ch==EOF) { 
      if (last_ch!='\n') { ++N; } // Capture non-newline-terminated last line.
      break;
    }

    last_ch=ch;
  }
  rewind(fin);
  return N;
}

void
Read_In_Barcodes(char* filename){
  /*
  Read in the barcodes from a given textfile, which needs to contain newline char seperated barcode entries.
  All barcodes must be the same length.
  Barcodes get read into an a_barcode struct, along with other information,
  and the struct is then stored in the barcodes array.
  filename: a string of the file name to read from
  return: void.
  post-condition: barcodes contains structs of all the barcodes in the given file.
  */
  FILE *fin; /// a pointer to the file object
  char *line = NULL;
  size_t len = 1000;
  char *readline; /// a empty container for a pointer to readline

  fin = fopen(filename,"r");

  /// Getting number of lines in the file.
  num_barcode = Get_Lines_In_File(fin);
  barcodes=(a_barcode**)R_alloc(num_barcode+1, sizeof(a_barcode*));

  line = (char *)malloc(sizeof(char) * (len+1)); /// allocate space for the line and assign it to line
  a_barcode *new_barcode; /// create a new_barcode struct variable
  int count = 0; 
  char * token; 

  while ((readline = fgets(line, len, fin)) != NULL){
    count++;
    new_barcode = (a_barcode *)malloc(sizeof(a_barcode));
    new_barcode->sequence = (char *)malloc(barcode_length * sizeof(char));

    strncpy(new_barcode->sequence, line, barcode_length); // copy the barcode line up to the length of the barcode into the barcode struct
    new_barcode->original_pos = count;
    if (is_PairedReads > 0) {
      // strtok returns the first token from the string up to the first sep \t. Subsequent calls return the next token
      token = strtok(line, "\t");
      token = strtok(NULL, "\t");
      new_barcode->sequenceRev = (char *)malloc(barcode_length_rev * sizeof(char));
      strncpy(new_barcode->sequenceRev, token, barcode_length_rev);
    } else if (is_DualIndexingReads > 0) {
      token = strtok(line, "\t");
      token = strtok(NULL, "\t");
      new_barcode->sequence2 = (char *)malloc(barcode_length_rev * sizeof(char));
      strncpy(new_barcode->sequence2, token, barcode2_length);
    } else {
      new_barcode->sequenceRev = NULL;
      new_barcode->sequence2 = NULL;
    };
    barcodes[count] = new_barcode;
  }
  fclose(fin);
  free(line);

  Rprintf(" -- Number of Barcodes : %d\n", num_barcode);
}

void
Read_In_Hairpins(char *filename){
  /*
  Read in the hairpins from a given file, storing each hairpin in an a_hairpin struct
  which is then referenced in the hairpins array.
  filename: the string of the file to read
  return: void.
  */
  FILE *fin;
  char * line = NULL;
  size_t len = 1000;
  char *readline;

  fin = fopen(filename,"r");

  // Getting number of lines in the file.
  num_hairpin = Get_Lines_In_File(fin);
  // allocates a pointer to an array of pointers a a_hairpin structs
  hairpins=(a_hairpin**)R_alloc(num_hairpin+1, sizeof(a_hairpin*));

  // allocates len+1 bytes of memory as a char* 
  line = (char *)malloc(len+1);
  a_hairpin *new_hairpin;
  int count = 0;

  while ((readline = fgets(line, len, fin)) != NULL){ // fgets reads chars into 'line', with max length of 'len' from 'fin'
    count++;
    new_hairpin = (a_hairpin *)malloc(sizeof(a_hairpin));
    new_hairpin->sequence = (char *)malloc(hairpin_length * sizeof(char)); // allocated space for our char array 
    new_hairpin->original_pos = count;
    strncpy(new_hairpin->sequence, line, hairpin_length); // copy the data retrieved (stored in line) to our new_hairpin struct
    hairpins[count] = new_hairpin; // store the pointer location in our array of hairpin structs
  }
  fclose(fin);
  free(line);

  Rprintf(" -- Number of Hairpins : %d\n", num_hairpin);
}

void
Check_Hairpins(void){
  /*
  Checks all the read in hairpins to ensure that they contain only A,T,G,C characters
  */
  int p, q;
  char base;
  for(p = 1; p <= num_hairpin; p++){
    for(q = 0; q < hairpin_length; q++){
      base = hairpins[p]->sequence[q];
      if ((base != 'A') && (base != 'T') && (base != 'G') && (base != 'C')){
	      Rprintf("Hairpin no.%d: %s contains invalid base %c\n", p, hairpins[p]->sequence, base);
      }
    }
  }
}

// Functions for locating barcodes and hairpins using a Trie initially, and reverting to interative mismatch search if the Trie yields no results
int
locate_sequence_in_trie(trie_node *trie_head, char *read, int *found_position) {
  /*
  Search through this read until we locate a known sequence in the given trie. Return that sequences' index.
  Otherwise, if we reach the end of the read, return -1
  trie_head: the head of the trie to search, can be either a barcode trie or a hairpin trie.
  read: the fastq read to search through
  found_position: a pointer to an int which will be used as an out parameter, to store the position in the read the sequence was found at
  return: the index of the read in the hairpins or barcode array, or -1 if not found
  */
  long read_length = strlen(read);
  int i, j;
  char base;
  trie_node *current_node;
  end_node *end;
  for (i = 0; i < read_length; i++) {
    current_node = trie_head;
    // search from i until we find a TERMINATOR
    for (j = i; j < read_length; j++) {
      base = read[j];
      if (Base_In_Node(current_node, TERMINATOR)) {
        // IF the current node can be terminated, return the found hairpin.
        current_node = current_node->links[Get_Links_Position(TERMINATOR)];
        end = current_node->end;
        *found_position = i;
        return end->sequence_index;
      } else if (Base_In_Node(current_node, base)) {
        // If we can continue traversing the trie, move to the next node
        current_node = current_node->links[Get_Links_Position(base)];
      } else {
        // else, we have to start searching from the next position in the read.
        break;
      }
    }
    // last check if we've reached the end of the read, and the TERMINATOR node exists
    if (Base_In_Node(current_node, TERMINATOR)) {
      current_node = current_node->links[Get_Links_Position(TERMINATOR)];
      end = current_node->end;
      *found_position = i;
      return end->sequence_index;
    }
  }
  *found_position = -1;
  return -1;
}

bool
Valid_Match(char *sequence1, char *sequence2, int length, int threshold){
  /*
  Traverse each sequence and count the number of incorrect bases.
  Naive approach. Doesn't take into account insertions or deletions.
  Fast. 
  For barcodes and hairpins
  sequence1: the first sequence to compare
  sequence2: the second sequence to compare
  length: the length of the sequences to compare up to
  threshold: the number of mismatches allowed, after which we return false
  return: true if a valid match, false if not
  */
  int i_base;
  int mismatchbasecount = 0;
  for (i_base = 0; i_base < length; i_base++) {
    if (sequence1[i_base] != sequence2[i_base]) {
      mismatchbasecount++;
      if (mismatchbasecount > threshold) {
        return false;
      }		  
    }
  }
  if (mismatchbasecount <= threshold) {
    return true;
  } else {
    return false;
  }
}

bool
Same_Forward_Barcode(int barcode1_index, int barcode2_index) {
  /*
  Determines if the given barcode indices in the barcodes array have the same forward barcode

  barcode1_index: the index of the first barcode to compare, in barcodes array
  barcode2_index: the index of the second barcode to compare.
  return: true if the same forward barcode, false otherwise
  */  
  if (barcode1_index  < 1 || barcode2_index < 1) {
    return false;
  }
  if (strcmp(barcodes[barcode1_index]->sequence, barcodes[barcode2_index]->sequence) == 0) {
    return true;
  }
  return false;
}

int
locate_mismatch_barcode_single(char *a_barcode) {
  /*
  Finds if the given barcode appears in the barcode array, with allowance
  for mismatches with bases
  a_barcode: the barcode to find
  return: the index of the found barcode, or -1 if not found
  */
 int i;
 for (i = 1; i <= num_barcode; i++) {
   if (Valid_Match(a_barcode, barcodes[i]->sequence, barcode_length, barcode_n_mismatch)) {
     return barcodes[i]->original_pos;
   }
 }
 return -1;
}

int
locate_mismatch_barcode_paired(char *a_barcode, char *a_barcode_rev) {
  /*
  Mismatch searching for the input barcodes in the barcode array.
  Allows for slight variation in the bases appearing in each barcode.
  a_barcode: the forward barcode read
  a_barcode_rev: the reverse barcode read for paired end matching
  return: the index of the found barcode, or -1 if not found.
  */
  int i;
  for (i = 1; i <= num_barcode; i++){
    if ((Valid_Match(a_barcode, barcodes[i]->sequence, barcode_length, barcode_n_mismatch)) && 
        (Valid_Match(a_barcode_rev, barcodes[i]->sequenceRev, barcode_length_rev, barcode_n_mismatch))) {
        return barcodes[i]->original_pos;
    }
  }
  return -1;
}

int
binary_search_barcode_paired(char *a_barcode, char *a_barcode_rev) {
  /*
  Using the given forward and reverse read barcodes, binary search the sorted
  barcodes array to locate the original position (the position that barcode was in in the unsorted array)
  and return it.
  a_barcode: the forward barcode read
  a_barcode_rev: the reverse barcode read, matched to the second input read.
  return: the original position of the barcode forward and reverse combination, or -1 if not found.
  */
  int imin, imax, imid;
  imin = 1;
  imax = num_barcode;

  while (imax >= imin) {
    imid = (imax + imin) / 2;
    // compare forward barcode sequence
    if (strncmp(barcodes[imid]->sequence, a_barcode, barcode_length) < 0) {
      imin = imid + 1;
    } else if (strncmp(barcodes[imid]->sequence, a_barcode, barcode_length) > 0) {
      imax = imid - 1;
    } else {
      // same forward sequence, compare reverse barcode sequence
      if (strncmp(barcodes[imid]->sequenceRev, a_barcode_rev, barcode_length_rev) < 0) {
        imin = imid + 1;
      } else if (strncmp(barcodes[imid]->sequenceRev, a_barcode_rev, barcode_length_rev) > 0) {
        imax = imid - 1;
      } else {
        return barcodes[imid]->original_pos;     
      } 
    }    
  }

  return -1;
}

int
binary_search_barcode_dualindex(char *a_barcode, char *a_barcode2) {
  /*
  Using the given forward and 2nd  barcodes, binary search the sorted
  barcodes array to locate the original position (the position that barcode was in in the unsorted array)
  and return it.
  a_barcode: the forward barcode read
  a_barcode2: the 2nd barcode read
  return: the original position of the barcode forward and reverse combination, or -1 if not found.
  */
  int imin, imax, imid;
  imin = 1;
  imax = num_barcode;

  while (imax >= imin) {
    imid = (imax + imin) / 2;
    // compare forward barcode sequence
    if (strncmp(barcodes[imid]->sequence, a_barcode, barcode_length) < 0) {
      imin = imid + 1;
    } else if (strncmp(barcodes[imid]->sequence, a_barcode, barcode_length) > 0) {
      imax = imid - 1;
    } else {
      // same forward sequence, compare reverse barcode sequence
      if (strncmp(barcodes[imid]->sequence2, a_barcode2, barcode2_length) < 0) {
        imin = imid + 1;
      } else if (strncmp(barcodes[imid]->sequence2, a_barcode2, barcode2_length) > 0) {
        imax = imid - 1;
      } else {
        return barcodes[imid]->original_pos;     
      } 
    }    
  }

  return -1;
}

int
locate_barcode(char *read, int *found_position) {
  /*
  Using the pre-built barcode trie, locate a barcode in the read.
  Iteratively search the trie from each position in the read.
  If we don't find a barcode in the read, and mismatch is enabled,
  perform mismatch searching of every index in the read.
  read: the read string to search for a barcode
  return: the index of the barcode in the barcode array, or -1 if not found.
  */
  int barcode_index = locate_sequence_in_trie(barcode_single_trie_head, read, found_position);
  if (barcode_index > 0) {
    return barcode_index;
  }

  // search the read for a mismatched barcode
  if (allow_mismatch > 0) {
    int i;
    int read_length = strlen(read);
    char *a_barcode = (char *)malloc(barcode_length * sizeof(char));
    for (i = 0; i < read_length - barcode_length; i++) {
      strncpy(a_barcode, read + i, barcode_length);
      barcode_index = locate_mismatch_barcode_single(a_barcode);

      if (barcode_index > 0) {
        *found_position = i;
        return barcode_index;
      }
    }
  }
  // if not found, set found_position as -1
  *found_position = -1;
  return -1;
}

int
locate_barcode_paired(char *read, char *read_rev, int *found_position) {
  /*
  Locate a barcode pair in the given forward and reverse reads,
  initially using trie matching to find barcode strings in the reads,
  and the calling binary_search_barcode_paired to locate the original position
  of the barcode pairing (if it existed)
  read: the forward read to scan for a barcode
  read_rev: the reverse read
  found_position: a pointer to an int in which to store the position of the first barcode found in 'read'
  return: the original position of the barcode forward and reverse match if found, -1 otherwise.
  found_position is set as -1 if otherwise.
  */
  int found_rev_position = 0;
  int found_for_position = 0;
  int barcode1_index, barcode2_index, found;
  char *barcode1, *barcode2;

  // match each barcode using trie matching, and copy the barcode from the read to local variables,
  // to call binary search on.
  barcode1_index = locate_sequence_in_trie(barcode_single_trie_head, read, &found_for_position);
  if (barcode1_index < 1) {
    // if we didn't find a barcode, terminate
    *found_position = -1;
    return -1;
  }

  barcode2_index = locate_sequence_in_trie(barcode_paired_trie_head, read_rev, &found_rev_position);
  if (barcode2_index > 0) {
    // if we found a second barcode, and thus a first barcode, terminate returning required information.
    barcode1 = (char *)malloc(barcode_length * sizeof(char));
    strncpy(barcode1, read + found_for_position, barcode_length);
    barcode2 = (char *)malloc(barcode_length_rev * sizeof(char));
    strncpy(barcode2, read_rev + found_rev_position, barcode_length_rev);

    found = binary_search_barcode_paired(barcode1, barcode2);
    if (found > 0) {
        *found_position = found_for_position;
        return found;
    }
  }

  // very inefficient mismatching. Improve with wildcard trie matching
  if (allow_mismatch > 0) {
    int i, j;
    int read_length = strlen(read);
    int read_rev_length = strlen(read_rev);
    char *a_barcode = (char *)malloc(barcode_length * sizeof(char));
    char *a_barcode_rev = (char *)malloc(barcode_length_rev * sizeof(char));
    for (i = 0; i < read_length - barcode_length; i++) {
      strncpy(a_barcode, read + i, barcode_length);
      if (locate_mismatch_barcode_single(a_barcode) > 0) {
        for (j = 0; j < read_rev_length; j++) {
          strncpy(a_barcode_rev, read_rev + j, barcode_length_rev);
          found = locate_mismatch_barcode_paired(a_barcode, a_barcode_rev);

          if (found > 0) {
            *found_position = i;
            return found;
          }
        }
      }
    }
  }

  *found_position = -1;
  return -1;
}

int
locate_barcode_dualIndexing(char *read, int *found_position){
  /* 
  Using a similar procedure to paired matching for barcodes, search for the first barcode in the read
  And if we found one, search for the second in the read from the end of that barcode onwards.
  This cuts down our search space to improve performance.
  read: the sequence to search for barcodes in
  found_position: a pointer to the int to store the found position of the first barcode in
  return: the original position index of the found barcode pair, or -1 if none found.
  post-condition: found_position contains the index of the start of the first barcode in the read, or 0 if none found
  */
  int found_for_position = 0;
  int found_dual_position = 0;
  int barcode1_index, barcode2_index, found;
  char *barcode1, *barcode2;
  // locate the forward read barcode
  barcode1_index = locate_sequence_in_trie(barcode_single_trie_head, read, &found_for_position);
  if (barcode1_index < 1) {
    *found_position = -1;
    return -1;
  }

  // locate the second barcode in the forward read
  barcode2_index = locate_sequence_in_trie(barcode_dualindex_trie_head, read + *found_position, &found_dual_position);
  if (barcode2_index > 0) {
    // if we found a barcode2, and thus we've found a barcode1, binary search to find it in the barcodes array
    barcode1 = (char *)malloc(barcode_length * sizeof(char));
    strncpy(barcode1, read + found_for_position, barcode_length);
    barcode2 = (char *)malloc(barcode2_length * sizeof(char));
    strncpy(barcode2, read + found_dual_position, barcode2_length);

    found = binary_search_barcode_dualindex(barcode1, barcode2);
    if (found > 0) {
        *found_position = found_for_position;
        return found;
    }
  }

  if (allow_mismatch > 0) {
    int i, j;
    int read_length = strlen(read);
    char *a_barcode = (char *)malloc(barcode_length * sizeof(char));
    char *a_barcode2 = (char *)malloc(barcode2_length * sizeof(char));
    for (i = 0; i < read_length - barcode_length; i++) {
      strncpy(a_barcode, read + i, barcode_length);
      if (locate_mismatch_barcode_single(a_barcode) > 0) {
        for (j = 0; j < read_length; j++) {
          strncpy(a_barcode2, read + j, barcode2_length);
          found = locate_mismatch_barcode_paired(a_barcode, a_barcode2);

          if (found > 0) {
            *found_position = i;
            return found;
          }
        }
      }
    }
  }
  
  /* @TODO REMOVE THIS
  if (allow_mismatch > 0) {
    int i;
    for (i = 1; i <= num_barcode; i++){
      if ((Valid_Match(a_barcode, barcodes[i]->sequence, barcode_length, barcode_n_mismatch)) && 
	    (Valid_Match(a_barcode2, barcodes[i]->sequence2, barcode2_length, barcode_n_mismatch))) {
        return barcodes[i]->original_pos;
      }
    }
  }
  */
  *found_position = -1;
  return -1;
}


int
locate_mismatch_hairpin(char *a_hairpin){
  /*
  Using a given hairpin, compare to all hairpins in the hairpin array
  to find a match, using a given mismatch amount (the number of bases which can vary before
  we decide the hairpin is not a match)
  a_hairpin: the hairpin to match
  return: the index of the hairpin in the hairpins array, or -1 if not found
  */
  int i;
  for (i = 1; i <= num_hairpin; i++){
    if (Valid_Match(a_hairpin, hairpins[i]->sequence, hairpin_length, hairpin_n_mismatch)) {
      return hairpins[i]->original_pos;
    }
  }
  return -1;
}

int
locate_hairpin(char *read, int *barcode_found_position, int *hairpin_found_position) {
  /*
  This function is a handler function for finding an exact match hairpin, and if that fails
  and mismatch hairpin is enabled, calling the OLD mismatch hairpin function on every position in the read
  until we find a hairpin, or reach the end of the read.
  read: the fastq read to find a hairpin in
  found_position: a pointer to the position that the barcode of this read was found in
  return: the index of the read in the sorted hairpins array, or -1 if no hairpin was found
  */
  int barcode_start = *barcode_found_position;
  if (barcode_start == -1) {
    barcode_start = 1 - barcode_length;
  }
  // using pointer arithmetic, search the hairpin trie starting the read at the end of the found barcode (or the start of the read if barcode not found)
  int hairpin_index = locate_sequence_in_trie(hairpin_trie_head, read + barcode_start + barcode_length - 1, hairpin_found_position);
  if (hairpin_index > 0) {
    return hairpin_index;
  }

  // search the read for a mismatched hairpin
  // @TODO test the allow_mismatch functionality of this function
  if (allow_mismatch > 0) {
    int i;
    int read_length = strlen(read);
    char *a_hairpin = (char *)malloc(hairpin_length * sizeof(char));
    for (i = 0; i < read_length - hairpin_length; i++) {
      // copy the possible hairpin at this i position and check if a mismatch hairpin exists
      strncpy(a_hairpin, read + i, hairpin_length);
      hairpin_index = locate_mismatch_hairpin(a_hairpin);

      if (hairpin_index > 0) {
        *hairpin_found_position = i;
        return hairpin_index;
      }
    }
  }
  *barcode_found_position = -1;
  *hairpin_found_position = -1;
  return -1;
}


// Barcode sorting
int 
barcode_compare(a_barcode *barcode1, a_barcode *barcode2){
  /*
  Compares two barcodes based on their sequence, returning 
  0 for matches,
  negative if barcode1 < barcode2
  positive if barcode2 < barcode1

  barcode1: The first barcode to compare
  barcode2: The second barcode to compare
  return: an integer. 0 indicates matches, negative: barcode1 < barcode2, positive: barcode1 > barcode2
  */
  int ans;
  // strncmp returns 0 if identical, negative if arg1 is < arg2, and positive if arg2 < arg1
  ans = strncmp(barcode1->sequence, barcode2->sequence, barcode_length);
  if (ans == 0) {
    if (is_PairedReads > 0){  
      ans = strncmp(barcode1->sequenceRev, barcode2->sequenceRev, barcode_length_rev);
    } else if (is_DualIndexingReads > 0){
      ans = strncmp(barcode1->sequence2, barcode2->sequence2, barcode2_length);
    }
  }

  return ans;
}

void
Sort_Barcodes(void){
  /* 
  Sorts the barcodes in lexographical order
  Implements bubble sort to do so.
  return: void
  */
  int i, j;
  a_barcode *temp;
  for(i = 1; i < num_barcode; i++){
    for(j = i+1; j <= num_barcode; j++){
      if (barcode_compare(barcodes[i], barcodes[j]) > 0) {
        // if the barcodes[i] is > barcodes[j], swap them
        temp = barcodes[i];
        barcodes[i] = barcodes[j];
        barcodes[j] = temp;
      }
    }
  }
}

int
Base_to_Int(char* base) {
  /*
  Determine the position of a base in the array used for count sort
  base: the base to convert to an int
  return: an int in the range 0-4 inclusive
  */
  switch (*base) {
        case 'A':
          return 1;
          // positions[1]++ usingthis above returns the pre-incremented value of positions
          break;
        case 'C':
          return 2;
          break;
        case 'G':
          return 3;
          break;
        case 'T':
          return 4;
          break;
        case TERMINATOR:
        default:
          return 0;
      }
}

// Hairpin sorting
void
Count_Sort_Hairpins(long index, a_hairpin** input_hairpins, a_hairpin** sorted_hairpins) {
  /* 
  Implements Count Sort, which stable sorts the input_hairpins, making use of the intermediate
  sorted_hairpins array given, to store the sorted hairpins as we find their positions.
  
  index: the index of the hairpin to sort based on
  input_hairpins: a pointer to the unsorted array of hairpins
  sorted_hairpins: a pointer to a allocated array of a_hairpin structs, which may or may not contain pointers.
                  These will be overridden
  */
  long counts[5]; // counts for empty, A, C, G, T
  long positions[5]; // positions for empty A, C, G, T
  int arr_pos;
  // intialise the counts array to 0
  for (arr_pos = 0; arr_pos < 5; arr_pos++) {
    counts[arr_pos] = 0;
  }

  // Count each appearance of the bases in the array
  char base;
  for (arr_pos = 1; arr_pos <= num_hairpin; arr_pos++) {
      // count the number of occurances of each base
      base = input_hairpins[arr_pos]->sequence[index];
      counts[Base_to_Int(&base)]++;
  }

  // positions holds [whitespace, A, C, G, T]
  positions[0] = 1;
  // determine the positions of the first appearance of each base in the array.
  // position[i] = position[i-1] + count[i-1]
  for (arr_pos = 1; arr_pos < 5; arr_pos++) {
    positions[arr_pos] = positions[arr_pos-1] + counts[arr_pos-1];
  }

  // construct the sorted hairpins array from our stored positions
  for (arr_pos = 1; arr_pos <= num_hairpin; arr_pos++) {
      base = input_hairpins[arr_pos]->sequence[index];
      // using the position of this base, insert this hairpin into the sorted array, and increment that position
      sorted_hairpins[positions[Base_to_Int(&base)]++] = input_hairpins[arr_pos];
  }

  // now, sorted_hairpins should now contain all values from input_hairpins but sorted
  // make input_hairpins reflect the values of sorted_hairpins
  int j;
  for (j= 1; j <= num_hairpin; j++) {
    input_hairpins[j] = sorted_hairpins[j];
  }
}

void 
Sort_Hairpins(void) {
  /* 
  Implements Radix Sort, which performs count sort on an array
  of hairpins, on each subsequent base from right to left.
  At end, hairpins array contains the same structs but sorted lexographically
  Makes use of an intermediate sort temporary array, which holds references to existing structs.
  These should not be freed on function termination, as doing so will destroy our hairpin data.
    We only need to free the intermediate array pointers.
  */
  // Create our storage for the intermediate steps of count sort.
  a_hairpin** temporary_hairpins = (a_hairpin **)malloc((num_hairpin+1) * sizeof(a_hairpin*));

  long i;
  // run radix_sort on this local array of the hairpins, freeing memory as we go
  for (i = hairpin_length; i >= 0; i--) {
    // run count sort, which will sort based on the given index
    Count_Sort_Hairpins(i, hairpins, temporary_hairpins);
    // the hairpins array is now stable sorted based on the ith index
  }
  // free the array we created.
  free(temporary_hairpins);
}


void
Process_Hairpin_Reads(char *filename, char *filename2){
  /*
  Analyse the given files, searching for barcodes and hairpins and incrementing
  the related counts in the summary table, as well as recording the positions in each
  read the forward barcodes and hairpins were found in.
  filename: the forward read file to anaylse
  filename2: the optional reverse read file to analyse, only used if is_PairedReads is true.
  */
  FILE *fin = NULL; // the sample file
  FILE *finRev = NULL; // the dual read
  char *line = NULL;
  char *line2 = NULL;
  size_t len = 1000;
  char *readline;
  char *readline2;
  long num_read_thisfile = 0;

  // opens the given files and allocate space for the line
  line = (char *)malloc(sizeof(char) * (len+1));
  fin = fopen(filename,"r");
  if (is_PairedReads > 0) {
    finRev = fopen(filename2, "r");
    line2 = (char *)malloc(len+1 * sizeof(char));
  }

  if (isverbose > 0){
    if (is_PairedReads > 0) {
      Rprintf("Processing reads in %s and %s.\n", filename, filename2);
    } else {
      Rprintf("Processing reads in %s.\n", filename);
    }
  }
  long line_count = 0;

  int barcode_index = -1;
  int hairpin_index;
  int barcode_start_position = 0;
  int hairpin_start_position = 0;
  long this_read_length = 0;

  // analyze each line of the file given
  while ((readline = fgets(line, len, fin)) != NULL){
    if (is_PairedReads > 0) {
      readline2 = fgets(line2, len, finRev);
	    if (readline2 == NULL) {
	      break;
	    }
    }

    // fastq files are 4 lines per read, with the sequence being on the second line.
    line_count++;  

    if ((line_count % 4) != 2) {
      // if the barcodes are in the header of each fastq group, we need to search line_count % 4 == 1 for the barcode
      if ((line_count % 4 ) == 1) {
        if (barcodesInHeader > 0) {
          // search the header line for the barcode. New Implementations of barcode searching using header line should go here
          if (is_PairedReads > 0) {
            barcode_index = locate_barcode_paired(line, line2, &barcode_start_position);
            barcode_start_position = -1;
          } else {
            barcode_index = locate_barcode(line, &barcode_start_position);
            barcode_start_position = -1;
          }
        }
      }
      continue;
    }

    // store the length of this line if greater than previous lines
    this_read_length = strlen(line);
    if (this_read_length > longest_read_length) {
      longest_read_length = this_read_length;
    }

    // Print out the current number of reads to screen, mod 10 million
    if ((isverbose > 0) && (num_read_thisfile % BLOCKSIZE == 0)) {
      Rprintf(" -- Processing %d million reads\n", (num_read_thisfile / BLOCKSIZE + 1) * 10);
    }
    num_read++;
    num_read_thisfile++;
    
    // Match the barcodes based on the input arguments for the type of barcode matching
    if (is_PairedReads > 0){    
      // Using trie matching, find a matching barcode forward and reverse sequence in the two lines
      barcode_index = locate_barcode_paired(line, line2, &barcode_start_position);
    } else if (is_DualIndexingReads > 0) {
      // Uising trie matching, find a matching barcode forward and 2nd sequence in the single line
      barcode_index = locate_barcode_dualIndexing(line, &barcode_start_position);
    } else if (barcodesInHeader <= 0) { 
      // Nothing needs to be done if the barcodes are in the header, as they are already found
      // The main barcode matching for single read sequences. Uses trie matching to find the 'original index' of the barcode in the barcodes array
      // This index is used to increment a counter in the summary table
      barcode_index = locate_barcode(line, &barcode_start_position);
    }

    if (barcode_index > 0) {
      // Record the position this barcode was found in, in the read.
      barcodecount++;
      if (barcodesInHeader <= 0) {
        // We don't care about the position of the barcodes found if the barcodes are found in the header line
        barcode_positions_size = Increment_Resize_Array(&barcode_positions, barcode_positions_size, barcode_start_position); 
      }
    }

    // Using trie matching, find a valid hairpin within the read.
    hairpin_index = locate_hairpin(line, &barcode_start_position, &hairpin_start_position); 
    if (hairpin_index > 0){
      hairpincount++;
      hairpin_positions_size = Increment_Resize_Array(&hairpin_positions, hairpin_positions_size, hairpin_start_position);
    }
       
    // Increment the count for the specific barcode and hairpin
    if ((barcode_index > 0) && (hairpin_index > 0)) {
      summary[hairpin_index][barcode_index]++;
      bchpcount++;
    }

  }

  if (isverbose > 0) {
    if (is_PairedReads > 0) {
      Rprintf("Number of reads in file %s and %s: %ld\n", filename, filename2, num_read_thisfile);  
    } else {
      Rprintf("Number of reads in file %s : %ld\n", filename, num_read_thisfile);  
    }
  }
  
  fclose(fin); 
  free(line);

  if (is_PairedReads > 0){  
    fclose(finRev);
    free(line2);
  }
}


void 
Initialise(int IsPaired, int IsDualIndexing, 
           int barcodeLength, int barcode2Length, int barcodeLengthRev,
           int hairpinLength,
           int allowMismatch, int barcodemismatch, int hairpinmismatch, 
           int verbose,
           int BarcodesInHeader){
	/* 
  Initiliases all local variables with given values
  isPaired: determines whether two reads are given, for forward and reverse reads
  isDualIndexing: determines whether 2 barcode sequences should be found in each read
  barcodeLength: the length of each barcode
  barcode2Length: the length of the dualIndex barcodes
  barcodeLengthRev: the length of the reverse read barcodes
  hairpinLength: the haripin length
  allowMismatch: if sequences with a pre-determined base mismatch is allowed
  barcodemismatch: the number of mismatch bases allowed before the sequence does not match
  hairpinmismatch: same as barcodemismatch, for hairpins
  verbose: if the user expects extra text printing during execution
  BarcodesInHeader: determines if barcodes should be search for in the header line of each read 
    that is, the line immediately above the read line in a fastq file (the linecount % 4 == 1, if the first line is line 1)
  */
  num_barcode = 0;
  num_hairpin = 0;

  is_PairedReads = IsPaired;
  is_DualIndexingReads = IsDualIndexing;
  barcodesInHeader = BarcodesInHeader;
  barcode_length = barcodeLength;
  barcode2_length = barcode2Length;
  barcode_length_rev = barcodeLengthRev;
  hairpin_length = hairpinLength;

  allow_mismatch = allowMismatch;
  barcode_n_mismatch = barcodemismatch;
  hairpin_n_mismatch = hairpinmismatch;
  isverbose = verbose;

  num_read = 0;
  barcodecount = 0;
  hairpincount = 0;
  bchpcount = 0;

  long read_length = 50;
  longest_read_length = 0;
  barcode_positions = Initialise_Resize_Array(read_length);
  barcode_positions_size = read_length;
  hairpin_positions = Initialise_Resize_Array(read_length);
  hairpin_positions_size = read_length;
}

void
Output_Summary_Table(char *output){
  /* 
  Writes to the file given the data present in our out summary table
  Includes the header for that column

  output: the file to write to
  */
  Rprintf("outputting the summary table to file %s\n", output);
  int i, j;
  FILE *fout;
  fout = fopen(output, "w");
  for(i = 1; i <= num_hairpin; i++) {
    fprintf(fout, "%ld", summary[i][1]);
    for(j = 2; j <= num_barcode; j++) {
      fprintf(fout, "\t%ld", summary[i][j]);
    }
    fprintf(fout, "\n");
  }
  fclose(fout);
  Rprintf("Finished outputting the summary table\n");
}

void
Output_Sequence_Locations(char *output, long *arr, int size) {
  /* 
  Writes to the file given the data present in the resize_array indicating
  either hairpin read position finds or barcode read position finds.
  Includes the header for that column

  output: the file to write to
  */
  Rprintf("Outputting Sequence locs to file: %s\n", output);
  
  int j;
  long max_size;
  if (size < longest_read_length) {
    max_size = size;
  } else {
    max_size = longest_read_length;
  }

  FILE *fout;
  fout = fopen(output, "w");
  fprintf(fout, "%ld", arr[0]);
  for(j = 1; j < max_size; j++) {
    fprintf(fout, "\n%ld", arr[j]);
  }
  fprintf(fout, "\n");
  fclose(fout);
  Rprintf("Finished Outputting sequences locs\n");
}

/*
long
Calculate_Average_Position_Array(long *arr, int size) {
  
  Calculates the average value of an array, by summing and dividing by size.

  arr: the array to average
  size: the length of the array.
  return: the average value of the array.
  
  long sum = 0;
  long sum_weights = 0;
  int i;
  for (i = 0; i < size; i++) {
    sum += (arr[i] * i);
    sum_weights += arr[i];
  }
  return sum/sum_weights;
}
*/
void
Clear_Trie(trie_node *node) {
  /* 
  Recursive function to clear a trie. 
  Calls this function on each linked node, then frees this node

  node: a pointer to the current node to free

  */
    int i;
    if (node->end != NULL) {
        free(node->end);
    }
    for (i = 0; i < 5; i++) {
        if (node->links[i] != NULL) {
            //char link_base = node->links[i]->base;
            Clear_Trie(node->links[i]);
        }
    }
    free(node);
}

void
Clean_Up(void){
  /*
  Deallocate all space for arrays created
  */
  Rprintf("Cleaning up the data\n");
  int index;
  // free the barcode array
  for (index = 1; index <= num_barcode; index++){
    free(barcodes[index]->sequence);
    if (is_PairedReads > 0){
      free(barcodes[index]->sequenceRev);
    } 
    if (is_DualIndexingReads > 0){
      free(barcodes[index]->sequence2);
    }
    free(barcodes[index]);
  }

  // free the hairpin array
  for (index = 1; index <= num_hairpin; index++){
    free(hairpins[index]->sequence);
    free(hairpins[index]);
  }

  // free the summary table
  for (index = 0; index <= num_hairpin; index++){
    free(summary[index]);
  }
  free(summary);

  Rprintf("Clearing the Tries\n");
  //free the hairpin & barcode tries
  Clear_Trie(barcode_single_trie_head);
  if (is_PairedReads) {
    Clear_Trie(barcode_paired_trie_head);
  } else if (is_DualIndexingReads) {
    Clear_Trie(barcode_dualindex_trie_head);
  }
  Clear_Trie(hairpin_trie_head);
  
  Rprintf("Clearing the resize arrays\n");
  free(barcode_positions);
  free(hairpin_positions);
}

void
Allocate_Summary_Table(void){
  /*
  Allocate a 2D array of long, which are to store our barcode/hairpin read counts.
  summary[i][j] = the count of reads with hairpin i and barcode j
  */
  int i, j;
  
  summary = (long **)malloc((num_hairpin+1) * sizeof(long *));
  for (i=0; i <= num_hairpin; i++){
    summary[i] = (long *)malloc((num_barcode+1) * sizeof(long));
  }
 
  for (i = 0; i <= num_hairpin; i++){
    for (j = 0; j <= num_barcode; j++){ 
      summary[i][j] = 0;  
	}
  }
}

void 
processHairpinReads(int *isPairedReads, int *isDualIndexingReads, 
                    char **file, char **file2, int *filecount, 
                    char**barcodeseqs, char**hairpinseqs,
                    int *barcodeLength, int *barcode2Length, int *barcodeLengthRev,
                    int *hairpinLength,
                    int *allowMismatch, int *barcodemismatch, int *hairpinmismatch,
                    char **output, int *verbose, int *barcodesInHeader,
                    char **barcodePosFile, char **hairpinPosFile)
{  
  /* 
  The entry point for the processAmplicons function.
  This function reads in all barcode and hairpin data, and for every file given
  searches for matching hairpins and barcodes and records their ID, and the position
  found in each read.
  This data is made available to the master R function through the given output files.
  All input data is given as pointers which must be dereferenced, corresponding to the
  way R provides C functions with arguments.

  isPairedReads: denotes if 2 read files are provided, and 2 sets of barcodes
  isDualIndexingReads: denotes if 2 sets of barcodes are provided, both of which should be located in each read
  file: the main array of files to read.
  file2: the secondary array of files to read. (reverse reads)
  filecount: the number of files given
  barcodeseqs: the file containing all the barcode information
  hairpinseqs: the file containing all hairpin information
  barcodeLength, barcode2Length, barcodeLengthRev: the lengths of the corresponding barcode sequences (if provided)
  haripinLength: the length of the hairpin sequences
  allowMismatch: denotes if mismatch barcode and hairpin bases are allowed
  barcodemismatch, hairpinmismatch: denotes the number of bases allowed to mismatch before considered incorrect
  output: the main output file, to which the barcode and hairpin match counts will be recorded
  verbose: denotes if extra text output should be provided upon execution of the function
  barcodesInHeader: denotes if barcodes should be matched in the header of each read
  barcodePosFile, hairpinPosFile: the files to which the positions of barcode and hairpin matches should be recorded.
  */
  // retrieves all our pointer data and stores it as local variables
  Initialise(*isPairedReads, *isDualIndexingReads,
             *barcodeLength, *barcode2Length, *barcodeLengthRev, 
             *hairpinLength,
             *allowMismatch, *barcodemismatch, *hairpinmismatch, 
             *verbose, *barcodesInHeader);

  Read_In_Barcodes(*barcodeseqs); 

  // build our barcode trie based on paired reads or dual indexing
  if (is_PairedReads > 0) {
    barcode_paired_trie_head = Build_Trie_Barcodes(true, false);
  } else if (is_DualIndexingReads > 0) {
    barcode_dualindex_trie_head = Build_Trie_Barcodes(false, true);
  }
  // Always build the single read trie no matter the index method is, as locating the forward read barcode is always required
  barcode_single_trie_head = Build_Trie_Barcodes(false, false);
  Sort_Barcodes(); // bubble sort. Is there a better sort??

  Read_In_Hairpins(*hairpinseqs);
  Sort_Hairpins();  // radix sort
  // Checks if all hairpins only contain ATGC bases
  Check_Hairpins();
  hairpin_trie_head = Build_Trie_Hairpins(); 

  Allocate_Summary_Table();
  
  int i_file;

  // For each file given, run Process_Hairpin_Reads
  for (i_file = 0; i_file < *filecount; i_file++){
    Process_Hairpin_Reads(file[i_file], file2[i_file]);
  }
  
  Rprintf("\nThe input run parameters are: \n");
  Rprintf(" -- Barcode: length %d\n", barcode_length);  
  if (is_DualIndexingReads){
    Rprintf(" -- Second Barcode in forward read: length %d\n", barcode2_length); 
  }
  if (is_PairedReads){
    Rprintf(" -- Barcode in reverse read: length %d\n", barcode_length_rev); 
  }
  Rprintf(" -- Hairpin: length %d\n", hairpin_length); 

  if (allow_mismatch > 0) {
    Rprintf(" -- Allow sequence mismatch, <= %d base in barcode sequence and <= %d base in hairpin sequence. \n", barcode_n_mismatch, hairpin_n_mismatch );
  } else {
    Rprintf(" -- Mismatch in barcode/hairpin sequences not allowed. \n");
  } 

  Rprintf("\nTotal number of read is %ld \n", num_read);
  Rprintf("There are %ld reads (%.4f percent) with barcode matches\n", barcodecount, 100.0*barcodecount/num_read);
  Rprintf("There are %ld reads (%.4f percent) with hairpin matches\n", hairpincount, 100.0*hairpincount/num_read);
  Rprintf("There are %ld reads (%.4f percent) with both barcode and hairpin matches\n", bchpcount, 100.0*bchpcount/num_read);

  //Rprintf("The average position of the barcode matches was %ld\n", Calculate_Average_Position_Array(barcode_positions, barcode_positions_size));
  //Rprintf("The average position of the hairpin matches was %ld\n", Calculate_Average_Position_Array(hairpin_positions, hairpin_positions_size));
  Output_Summary_Table(*output);

  Output_Sequence_Locations(*barcodePosFile, barcode_positions, barcode_positions_size);
  Output_Sequence_Locations(*hairpinPosFile, hairpin_positions, hairpin_positions_size);

  Clean_Up();
}

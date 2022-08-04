// A1


/* DNA sequence mapper:
 *
 * Skeleton code written by Jianzhong Qi, April 2022
 * Edited by: HENRY ROUTSON 1307261, April 29 2022
 * 
 */

// algorithms are fun

// When using .0
// https://docs.microsoft.com/en-us/cpp/c-language/c-floating-point-constants

#include <stdio.h>
#include <math.h>
#include <string.h>

#define STAGE_NUM_ONE 1						  /* stage numbers */ 
#define STAGE_NUM_TWO 2
#define STAGE_NUM_THREE 3
#define STAGE_NUM_FOUR 4
#define STAGE_NUM_FIVE 5
#define STAGE_HEADER "Stage %d\n==========\n" /* stage header format string */

#define MAX_READ_ID_LENGTH 100				  /* maximum read ID length */
#define MAX_READ_LENGTH 100					  /* maximum read length */
#define MAX_NUM_READS 100					  /* maximum number of reads */
#define MAX_REF_LENGTH 1000					  /* maximum reference DNA length */

#define TRUE 1                                /* define booleans */
#define FALSE 0                               /* define booleans */
#define END_CHAR '\0'                         /* char at end of strings */
#define DEFAULT_INDEX 0                       /* for finding min max or other */

// STAGE 3
#define HIGHEST_ACCEPTABLE_ERROR_CHAR '*'     /* char to define error bases */
#define HIGH_ERROR_MASK '*'                   /* char to mask error bases */

// STAGE 4
#define BASES "ACGT"                          /* DNA bases */
#define NBASES (int) strlen(BASES)            /* number of DNA bases */

// STAGE 5
#define QUALITY_ASCII_LOWER_BOUND 33.0        /* assesment specific constants */
#define ERR_PROB_PAR 10.0                     /* assesment specific constants */   
#define MATCH_EQL_SCALAR log2(ERR_PROB_PAR)/ERR_PROB_PAR  /* precompute */
#define LOG2_POINT_25 -2                      /* precompute */
#define LOG2_1 0                              /* precompute */


typedef char read_id_t[MAX_READ_ID_LENGTH+1]; /* a read ID */
typedef char read_t[MAX_READ_LENGTH+1];		  /* a read */
typedef char score_t[MAX_READ_LENGTH+1];	  /* quality scores of a read */
typedef char ref_t[MAX_REF_LENGTH+1];		  /* a reference DNA sequence */


/****************************************************************/

/* function prototypes */
int take_one_read(read_t one_read, score_t scores_of_one_read);
void print_stage_header(int stage_num);

void stage_one(read_t one_read, score_t scores_of_one_read);
void stage_two(read_t reads[], score_t scores[], int *num_reads);
void stage_three(read_t reads[], score_t scores[], int num_reads);
void stage_four(ref_t ref_sequence);
void stage_five(read_t reads[], score_t scores[], int num_reads, 
	ref_t ref_sequence);

/****************************************************************/

/* own function prototypes */
int string_min_index(char string[]);
int index_of_base_with_smallest_quality_score(char string[]);
int count_chars(char string[], char character);
double string_average(char string[]);
double match_score_i(char R[], char S[], char S_scores[], int i);
double match_score(char R[], char S[], char S_scores[]);
void print_closest_subsequence(char sequence[], char sub[], char sub_scores[]);

/****************************************************************/


/* main function controls all the action, do NOT modify this function */
int
main(int argc, char *argv[]) {
	/* to hold all input reads and quality scores */
	read_t reads[MAX_NUM_READS];	
	score_t scores[MAX_NUM_READS];
	/* to hold the number of input reads */
	int num_reads = 0;	
	/* to hold the input reference sequence */
	ref_t ref_sequence;
	
	/* stage 1: process one read */
	stage_one(reads[0], scores[0]); 
	num_reads++;
	
	/* stage 2: process all reads */
	stage_two(reads, scores, &num_reads);
	
	/* stage 3: mask bases with high error probability */ 
	stage_three(reads, scores, num_reads);
	
	/* stage 4: process reference sequence */
	stage_four(ref_sequence);
	
	/* stage 5: map reads to the reference sequence */
	stage_five(reads, scores, num_reads, ref_sequence);
	
	/* all done; take some rest */
	return 0;
}

/* print stage header given stage number */
void 
print_stage_header(int stage_num) {
	printf(STAGE_HEADER, stage_num);
}


/****************************************************************/
/* own functions */

int string_min_index(char string[]) {
    /* return the index of the lowest integer value character */

    int string_len = strlen(string); 

    // return the smallest index among the tied values
    int min_i = 0;
    for (int cur_i = 1; cur_i < string_len; cur_i++) { 
    
        if (string[ min_i ] > string[ cur_i ]) {
            min_i = cur_i; 
        }
    }
    
    return min_i;
}


int index_of_base_with_smallest_quality_score(char string[]) {
    /* return the index of the base with the lowest quality score
       where the quality score is determined by character integer value */
    
    return string_min_index(string);
}


int count_chars(char string[], char character) {
    /* count how many times character occurs in the string A */

    int strlen_string = strlen(string);
    int count = 0;
    for (int i = 0; i < strlen_string; i++) {
        count += (string[i] == character); // add 1 to sum if char
    }
    return count;
}


double string_average(char string[]) {
    /* return the average value of the integer representations 
       in the given string  */
    
    int strlen_string = strlen(string);  
    double sum_string = 0;
    for (int i = 0; i < strlen_string; i++) { 
        sum_string += string[i]; 
    }
    return sum_string / strlen_string;
}


double log2_match_score_i(char R[], char S[], char S_scores[], int i) { 
    /* return log2 of the match score of two characters
       R[i] and S[i]
       as defined for the assignment   */
    
    /*    Precompute -----------
          for R[i] == S[i]
          to avoid floating point errors 
          and to improve performance

    (Using log power rule)
    
    log2( 1.0 / pow( 10.0, ( (S_scores[i] - 33.0) / 10.0 ) ) ) =
    -log2( pow( 10.0, ( (S_scores[i] - 33.0) / 10.0 ) ) ) =
    ( (S_scores[i] - 33.0) / 10.0 )  * -log2(10.0)  =
    (33.0 - S_scores[i])  *  log2(10.0) / 10.0                     */

    // R is a read  
    // S is a string tring to match a Substring of R
    // S_scores is the quality scores of the bases in S

    if (R[i] == S[i]) { 
        return (QUALITY_ASCII_LOWER_BOUND - S_scores[i]) * MATCH_EQL_SCALAR; 
    } 
    if (R[i] ==  '*') { return LOG2_POINT_25; }
    return LOG2_1;
}


double match_score(char R[], char S[], char S_scores[]) { 
    /* return the match score of two strings
       as defined for the assignment   */
       
    // R is a read
    // S is a string tring to match a Substring of R
    // S_scores is the quality scores of the bases in S

    int length_S = strlen(S); // S is always the correct length, R is not
    double sum = 0;
    for (int i = 0; i < length_S; i++) {
        sum -= log2_match_score_i(R, S, S_scores, i);
    }
    return sum;
}


void print_closest_subsequence(char sequence[], char sub[], char sub_scores[]) {
    /* print subsequence found in the sequence
       with the highest match score to sub */
    
    int strlen_sub = strlen(sub); 
    
    double cur_score, max_score = -INFINITY;
    int max_score_index = DEFAULT_INDEX; 
    
    // Create sliding window to find max match score
    
    // + 1 because sub starts on i
    int i_max_plus1 =  strlen(sequence) - strlen_sub + 1; 
    for (int i = 0; i < i_max_plus1; i++) {

        cur_score = match_score(sequence + i, sub, sub_scores);
        if (max_score < cur_score) { 
            max_score = cur_score;
            max_score_index = i;
        }
    }
    
    for (int i = 0; i < strlen_sub; i++) {
        putchar(*(sequence + max_score_index + i));
    }
}


/****************************************************************/



/* process a read record */
int 
take_one_read(read_t one_read, score_t scores_of_one_read) {
    /* update one_read and scores_of_one_read values
       returns 1 if there is a read, 0 if not               */

	read_id_t id;
    scanf("%s", id);
    if (id[0] == '#') { return FALSE; }
	scanf("%s", one_read);
	getchar();
	getchar(); // skip '+' and '\n'
	scanf("%s", scores_of_one_read);
    return TRUE;
}


/* stage 1: process one read */
void 
stage_one(read_t one_read, score_t scores_of_one_read) {
    /* print base with smallest quality score and it's index */

	print_stage_header(STAGE_NUM_ONE);
    
    // load values into one_read and score_of_one_read
    take_one_read(one_read, scores_of_one_read); 
    
    int indx = index_of_base_with_smallest_quality_score(scores_of_one_read);
    
    printf("Base with the smallest quality score: %c\n", one_read[indx]);
    printf("Index: %i\n", indx);
    printf("\n");
}


/* stage 2: process all reads */
void 
stage_two(read_t reads[], score_t scores[], int *num_reads) {
    /* print number of reads,
       smallest average read quality score 
       and the read with that average         */
    
    print_stage_header(STAGE_NUM_TWO);

    double average, smallest_average = INFINITY; 
    int smallest_average_index = DEFAULT_INDEX;
    
    // (*num_reads is index)
    // take_one_read returns 1 if there are reads, else 0
    while (take_one_read(reads[*num_reads], scores[*num_reads])) { 

        average = string_average(scores[*num_reads]);        
        if (average < smallest_average) {
            smallest_average = average;
            smallest_average_index = *num_reads;
        }
        
        ++*num_reads;
    }
    
    printf("Total number of reads: %i\n", *num_reads);
    printf("Smallest average quality score: %.2lf\n", smallest_average);
    printf("Read with the smallest average quality score:\n");
    printf("%s\n", reads[smallest_average_index]); 
    printf("\n");
}


/* stage 3: mask bases with high error probability */ 
void 
stage_three(read_t reads[], score_t scores[], int num_reads) {
    /* print reads,
       with high error probability bases masked with the '*' char    */

    print_stage_header(STAGE_NUM_THREE);    
    
        
    /*    Precompute -----------

    Q = scores[r][c]
    P(Q) > 0.15

    1 / 10 ^ ((Q - 33) / 10) > 0.15 = 15 / 100
    10 ^ ((Q - 33) / 10) < 100 / 15 = 20 / 3
    (Q - 33) / 10 < log10(20 / 3)
    Q < 10 * log10(20 / 3) + 33
    Q < 41.23908740944318757...

    Q is an integer
    Hence,
    Q < 42
    Q < '*'                                 */
    
    int c; // c is column index
    for (int r = 0; r < num_reads; r++) { // r is read / row index

        c = 0;
        while (scores[r][c] != END_CHAR) { 

            if (scores[r][c] < HIGHEST_ACCEPTABLE_ERROR_CHAR) {
                reads[r][c] = HIGH_ERROR_MASK;
            }
            c++;
        }
        printf("%s\n", reads[r]);
    }
    
    printf("\n");
}


/* stage 4: process reference sequence */
void stage_four(ref_t ref_sequence) {
    /* print number of bases in reference sequence
       for each base, print base and count            */
    
    print_stage_header(STAGE_NUM_FOUR);
    
    scanf("%s", ref_sequence);
    printf("Length of the reference sequence: ");
    printf("%i\n", (int) strlen(ref_sequence));

    // print count for each base
    for (int b = 0; b < NBASES; b++) {
        printf("Number of %c bases: ", BASES[b]);
        printf("%i\n", count_chars(ref_sequence, BASES[b]));
    }
    
    printf("\n");
}


/* stage 5: map reads to the reference sequence */
void 
stage_five(read_t reads[], score_t scores[], int num_reads, 
           ref_t ref_sequence) {
    /* print reads and their closest match */

    print_stage_header(STAGE_NUM_FIVE);
    
    // for each read
    for (int r = 0; r < num_reads; r++) {
        printf("Read:  ");
        printf("%s\n", reads[r]);
        printf("Match: ");
        print_closest_subsequence(ref_sequence, reads[r], scores[r]);
        printf("\n");
    }
}





















// A2


/* Restaurant recommender:
 *
 * Skeleton code written by Jianzhong Qi, April 2022
 * Edited by: HENRY ROUTSON 1307261
 *
 */

// algorithms are fun
// algorithms are awesome

// #define NDEBUG // uncomment to turn off assert statements

// structs 
#define MAX_CUISINE_NAME_LEN 20 // max number of characters in a cuisine name
#define MAX_REST_NAME_LEN 100 // max number of characters in a restaurant name
#define CHARS_IN_DATETIME 19 // number of characters in the datetime
#define MAX_RESTS 99 // max number of restaurants
#define MAX_CUSTOMER_ID_LEN 6 // number characters in customer ID

// stage 4
#define K_MOST_SIMILAR 2 // reccommends based on the k most similar customers

#include <stdio.h>
#include <string.h>
#include <math.h>

/* stage heading */
#define HEADING "==========================Stage %d==========================\n"

/* structs */

typedef struct {
  double x, y;
} coordinates_t;

typedef struct {
  int ID;
  int avg_price; // (per head)
  coordinates_t coordinates;
  char cuisine[MAX_CUISINE_NAME_LEN];
  char name[MAX_REST_NAME_LEN];
} restaurant_t;

typedef struct {
  char date_time[CHARS_IN_DATETIME];
  char customerID[MAX_CUSTOMER_ID_LEN];
  int restaurantID;
  int amount;
} transaction_t;

typedef struct {
  char ID[MAX_CUSTOMER_ID_LEN];
  int restaurants_visits[MAX_RESTS]; // (implicit counter)
} customer_t;

/* files */

typedef customer_t data_t;
#include "listops.c"

/* function prototypes */

void print_stage_header(int stage);

// stage 1
int write_restaurant_line_to(restaurant_t * restaurant_ad);

// stage 2
int write_transaction_line_to(transaction_t * transaction_ad);
void integrate_transaction(list_t * p_list, transaction_t * p_transaction,
                           restaurant_t * rests, int num_rests);
void print_customer(customer_t * p_rest, int num_rests);
void print_customers_list(node_t * p_node, int num_rests);

// stage 3
double distance(coordinates_t * p_r1, coordinates_t * p_r2);
int similar_rest(restaurant_t * p_rest1, restaurant_t * p_rest2);
void recommendations_3(int * rest_visits, restaurant_t * rests, int num_rests);

// stage 4
int visiting_similarity(int visits1[], int visits2[], int num_rests);
void recommendations_4(node_t * p_nd_Cust, node_t * p_nd_Othr, int num_rests);

/****************************************************************/
// own main

int main(int argc, char * argv[]) {

  // STAGE 1 -----------------------------
  print_stage_header(1);

  int num_rests = 0;
  int cheapest_index = 0;
  restaurant_t rests[MAX_RESTS];

  while (write_restaurant_line_to( & rests[num_rests])) {
  
    int cur_is_cheapest = 
        rests[num_rests].avg_price <  rests[cheapest_index].avg_price ||
        // OR
       (rests[num_rests].avg_price == rests[cheapest_index].avg_price &&
        rests[num_rests].ID < rests[cheapest_index].ID); // tie breaker

    if (cur_is_cheapest) {
      cheapest_index = num_rests; // update cheapest value
    }
    
    num_rests++;
  }

  printf("Number of restaurants: %i\n", num_rests);
  printf("Restaurant with the smallest average price: ");
  printf("%s\n\n", rests[cheapest_index].name);

  scanf("#####");
  

  // STAGE 2 -----------------------------
  print_stage_header(2);

  list_t * p_customer_list = make_empty_list();

  transaction_t transaction;
  while (write_transaction_line_to( & transaction)) {
    integrate_transaction(p_customer_list, & transaction, rests, num_rests);
  }

  node_t * p_head = p_customer_list -> head;
  print_customers_list(p_head, num_rests);
  printf("\n");
  

  // STAGE 3 -----------------------------
  print_stage_header(3);

  node_t * p_node = p_head;
  while (p_node) {
    recommendations_3(p_node->data.restaurants_visits, rests, num_rests);
    p_node = p_node -> next;
  }
  print_customers_list(p_head, num_rests);
  printf("\n");


  // STAGE 4 -----------------------------
  print_stage_header(4);

  p_node = p_head;
  while (p_node) {  // "for each customer C"
    recommendations_4(p_node, p_head, num_rests);
    p_node = p_node -> next;
  }
  print_customers_list(p_head, num_rests);
    
    
  /*
  STAGE 4 COMPLEXITY ----
  
      O(c^2 * r * k) run
      O(k) mem
  
  */

  // -------------------

  return 0;
}

/* print stage header given stage number */
void
print_stage_header(int stage_num) {
  printf(HEADING, stage_num);
}

/****************************************************************/
// own functions

// stage 1 ---

int write_restaurant_line_to(restaurant_t * p_restaurant) {
  // reads line and writes to address of restaurant struct
  // 699229 19.3 55.4 30 indian curry_palace

  if (scanf("%6i", & (p_restaurant -> ID)) == 0) { return 0; }
  scanf("%lf", & (p_restaurant -> coordinates.x));
  scanf("%lf", & (p_restaurant -> coordinates.y));
  scanf("%i", & (p_restaurant -> avg_price));
  scanf("%s", (char * ) & (p_restaurant -> cuisine));
  scanf("%s", (char * ) & (p_restaurant -> name));

  return 1;
}

// stage 2 ---

int write_transaction_line_to(transaction_t * p_transaction) {
  // reads line and writes to address of transaction struct
  // 2022:04:05:17:05:55 jxrpfj 190947 42

  if (getchar() == -1) { return 0; } // getchar returns '\n' or -1 if EOF
  scanf("%s ", (char * ) & (p_transaction -> date_time));
  scanf("%6s ", (char * ) & (p_transaction -> customerID));
  scanf("%6i ", & (p_transaction -> restaurantID));
  scanf("%i", & (p_transaction -> amount));

  return 1;
}

void integrate_transaction(list_t * p_list, transaction_t * p_transaction,
                           restaurant_t * rests, int num_rests) {
  /* Add transaction data to customer list
  Args:
      p_list  Pointer to customer list
      p_transaction  Pointer to transaction 
      rests  Array of restaurant structs
      num_rests  Number of restaurants    */

    
  // find customer with transaction customer ID 
  node_t * p_node = p_list -> head;
  while (1) {

    if (p_node == '\0') {
      // or add customer if not in linked list

      customer_t customer;
      strcpy(customer.ID, p_transaction -> customerID);
      memset(customer.restaurants_visits, 0, sizeof(int) * MAX_RESTS);

      insert_at_foot(p_list, customer);
      p_node = p_list -> foot;
      break;
    }

    int string_dif = strcmp(p_node -> data.ID, p_transaction -> customerID);
    if (string_dif == 0) { break; }

    p_node = p_node -> next;
    
  }

  assert(strcmp(p_node -> data.ID, p_transaction -> customerID) == 0); 

  // find restaurant with transaction restaurant ID 
  int i = 0;
  while (p_transaction -> restaurantID != rests[i].ID) {
    i++; 
    assert(i < num_rests);
  }
  
  // update number of restaurant visits
  p_node -> data.restaurants_visits[i] += 1;

}

void print_customer(customer_t * p_customer, int num_rests) {
  // print a customer struct, with their ID and restaurants_visits

  printf("%s:", p_customer -> ID);

  for (int i = 0; i < num_rests; i++) {
    printf("  ");
    int to_print = p_customer -> restaurants_visits[i];
    switch (to_print) {
      case -1: printf("-"); break;
      case -2: printf("+"); break;
      case -3: printf("*"); break;
      default: printf("%i", to_print);
    }

  }
  printf("\n");

}

void print_customers_list(node_t * p_node, int num_rests) {
  /* Print a linked list of customers
  Args:
      p_node  pointer to head node in customer list
      num_rests  Number of restaurants (to print restaurants_visits)    */ 


  while (p_node) {
    print_customer( & p_node -> data, num_rests);
    p_node = p_node -> next;
  }

}

// stage 3 ---

double distance(coordinates_t * p_r1, coordinates_t * p_r2) {
  // calculate distance between coordinates with Pythagorean theorem 
    
  double x_dif = (p_r1 -> x) - (p_r2 -> x);
  double y_dif = (p_r1 -> y) - (p_r2 -> y);
  return pow(pow(x_dif, 2) + pow(y_dif, 2), 0.5);
}

int similar_rest(restaurant_t * p_rest1, restaurant_t * p_rest2) {
  /* Return Boolean for if restaurant are similar 
  Args:
      p_rest1  Restaurant struct
      p_rest2  
  Returns:
      1 if restaurants are similar, 0 otherwise
      under conditions in assignment    */

  int is_similar =
    abs((p_rest1 -> avg_price) - (p_rest2 -> avg_price)) <= 20 || // OR
    strcmp(p_rest1 -> cuisine, p_rest2 -> cuisine) == 0    || // OR
    distance( & p_rest1 -> coordinates, & p_rest2 -> coordinates) <= 30;

  return is_similar;
}

void recommendations_3(int * rest_visits, restaurant_t * rests, int num_rests) {
  /* Updates a customers restaurants_visits to indicate recommendations
  using restaurant similarity
  Args:
      visits1  Array of integer number of visits for each restaurant
      visits2  
      num_rests  Number of restaurants / elements in visits arrays    */

  int r;
  // for each visited restaurant 
  for (int R = 0; R < num_rests; R++) {
    if (rest_visits[R] > 0) {

      // for each unvisited restaurant similar to the visited one
      for (r = 0; r < num_rests; r++) {
        if (r != R && !rest_visits[r] && similar_rest(rests + R, rests + r)) {

          // change unvisited value of 0 to - to recommend
          assert(rest_visits[r] == 0);
          rest_visits[r] = -1;

        }
      }
    }
  }
}

// stage 4 ---

int visiting_similarity(int visits1[], int visits2[], int num_rests) {
  /* 
  Args:
      visits1  Array of integer number of visits for each restaurant
      visits2  
      num_rests  Number of restaurants / elements in visits arrays
  Returns: visiting similarity as defined in assignment    */
    
  int sum = 0;
  for (int i = 0; i < num_rests; i++) {
    if (visits1[i] > 0 && visits2[i] > 0) {
      sum += visits1[i] * visits2[i];
    }
  }
  return sum;
}

int min_index(int A[], int n) {
  /* 
  Args:
      A  An array
      n  Number of elements in A
  Returns: integer index in A where A[index] is min in A    */
  
  int index = --n;
  while (n--) {
    if (A[index] > A[n]) { index = n; }
  }
  return index;
}


/* 
                    Cannot find the error in stage 4
*/

void recommendations_4(node_t * p_nd_Cust, node_t * p_nd_Othr, int num_rests) {
  /* Updates a customers restaurants_visits to indicate recommendations
  using customer similarity
  Args:
      p_nd_Cust  Current customer node
      p_nd_Othr  Customer linked list head node
      num_rests  Number of restaurants / elements in visits arrays    */
    

  int *Cust_visits = p_nd_Cust -> data.restaurants_visits;

  int *highest_sims_visits[K_MOST_SIMILAR];
  memset(highest_sims_visits, '\0', sizeof(int *) * K_MOST_SIMILAR);

  int  highest_sims[K_MOST_SIMILAR];
  memset(highest_sims, 0, sizeof(int) * K_MOST_SIMILAR);

  // Find most similar customers

  while (p_nd_Othr) {
    if (p_nd_Othr != p_nd_Cust) {
      int *Othr_visits = p_nd_Othr -> data.restaurants_visits;

      int sim = visiting_similarity(Cust_visits, Othr_visits, num_rests);
      if (sim > 0) {

        int i = min_index(highest_sims, K_MOST_SIMILAR);
        if (highest_sims[i] < sim) {
          highest_sims[i] = sim;
          highest_sims_visits[i] = Othr_visits;
        }
      }
    }

    p_nd_Othr = p_nd_Othr -> next;
  }

  // Make recommendations

  for (int k = 0; k < K_MOST_SIMILAR; k++) {

    int *Otr_visits2 = highest_sims_visits[k];
    assert(Otr_visits2 != Cust_visits); // (pointers)
      
    if (Otr_visits2 != '\0') {  
      assert(highest_sims[k] > 0);
    
      // asserts
      for (int k2 = 0; k2 < K_MOST_SIMILAR; k2++) { 
        if (k != k2) { 
          assert(highest_sims_visits[k] != highest_sims_visits[k2]); 
        }
      }

      for (int r = 0; r < num_rests; r++) {

        // if restaurant is unvisited, but visited by other similar
        if (Cust_visits[r] <= 0 &&  Otr_visits2[r] > 0) {
          Cust_visits[r] -= 1; // recommend
        }
      }
        
    }
  }
}


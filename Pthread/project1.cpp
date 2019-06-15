#include <limits.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

typedef struct Node {
	struct Node* next;
	int value;			// process working time
}Node;

// aggregate variables
long sum = 0;
long odd = 0;
long min = INT_MAX;
long max = INT_MIN;
bool done = false;

// gloval variable
int idle_thread_cnt = 0;

// mutex
pthread_mutex_t lock;
pthread_mutex_t glock;
pthread_mutex_t sync_lock;
pthread_cond_t sync_cond;
pthread_cond_t done_cond;
pthread_cond_t arrive_cond;

// task queue
Node* task_head = NULL;

// function prototypes
void *update(void*);
void init();
void inqueue(Node **, int);
int dequeue(Node **);
void display_list(Node *);
void head_init(Node **);

/*
 * initiate objects
 */
void init(){
	head_init(&task_head);
	pthread_mutex_init(&lock, NULL);
	pthread_mutex_init(&glock, NULL);
	pthread_mutex_init(&sync_lock, NULL);
	pthread_cond_init(&sync_cond, NULL);
	pthread_cond_init(&done_cond, NULL);
	pthread_cond_init(&arrive_cond, NULL);
}

/*
 * clean up objects
 */
void cleanup(){
	pthread_mutex_destroy(&lock);
	pthread_mutex_destroy(&glock);
	pthread_mutex_destroy(&sync_lock);
	pthread_cond_destroy(&sync_cond);
	pthread_cond_destroy(&arrive_cond);
}


/*
 * Queue operation
 * inqueue, dequeue, display_list, head_init
 */
void inqueue(Node** head, int value){
	Node* new_node = (Node *)malloc(sizeof(Node));
	new_node->value = value;
	new_node->next=NULL;
	Node* current_node = (*head);
	while(current_node){
		if(current_node->next == NULL){
			current_node->next = new_node;
			break;
		}
		current_node = current_node->next;
	}
	return;
}

int dequeue(Node** head){
	Node* temp_node = NULL;
	Node* current_node = (*head);

	if((*head)->next == NULL) printf("fuck you\n");
	temp_node = current_node->next;
	int retval = current_node->next->value;
	current_node->next = current_node->next->next;
	if((temp_node)!=NULL)
		free(temp_node);

	return retval;
}

void display_list(Node* head) {
	Node* current_node = head->next;
	while (current_node){
		printf("[%d]-->", current_node->value);
		current_node = current_node->next;
	}
	printf("NULL\n");
}

void head_init(Node** head){
	(*head) = (Node *)malloc(sizeof(Node));
	(*head)->next = NULL;
}

int is_queue_empty(Node** head){
	if(((*head)->next) == NULL)
		return 1;
	return 0;
}


/*
 * update global aggregate variables given a number
 * parallelized version.
 */
void *update(void*)
{
	int number;
	pthread_t tid = pthread_self();

	// for pthread creation sync
	pthread_mutex_lock(&sync_lock);
	pthread_cond_signal(&sync_cond);
	pthread_mutex_unlock(&sync_lock);

	// wait for signal arrive!
	pthread_mutex_lock(&lock);
	while(is_queue_empty(&task_head) == 1 && !done){
		pthread_cond_wait(&arrive_cond, &lock);
	}

	if(done)
		goto here;

	do{			
		// work while master thread work

		// get a work from the task queue	
		idle_thread_cnt--;				// decrease idle thread count
		number = dequeue(&task_head);
		pthread_mutex_unlock(&lock);

    	// simulate computation
		sleep(number);

   		// update aggregate variables
		pthread_mutex_lock(&glock);
		sum += number;
		if (number % 2 == 1) {
			odd++;
		}
		if (number < min) {
			min = number;
		}
		if (number > max) {
			max = number;
		}
		pthread_mutex_unlock(&glock);

		pthread_mutex_lock(&lock);
		// wait for signal arrive!
		idle_thread_cnt++;
		while(is_queue_empty(&task_head) == 1 && !done){
			pthread_cond_wait(&arrive_cond, &lock);
		}
	}while(!done);
here :
	pthread_mutex_unlock(&lock);
	return NULL;
}

int main(int argc, char* argv[])
{
    // check and parse command line options
    if (argc != 3) {
        printf("Usage: sum <infile> <num_threads>\n");
        exit(EXIT_FAILURE);
    }
    char *fn = argv[1];
	int num_threads = atoi(argv[2]);

	// initialize mutex, cond, gval
	init();
	idle_thread_cnt = num_threads;		// when the num_threads equals 0, it means all threads is working.

    // load numbers and add them to the queue
    FILE* fin = fopen(fn, "r");

	pthread_t worker_threads[num_threads];
	int i;

	for(i = 0 ; i < num_threads ; i++){
			pthread_mutex_lock(&sync_lock);
			if(pthread_create(&worker_threads[i], NULL, update, NULL) != 0){
				printf("PTHREAD_CREATE FAILURE\n");
            	exit(EXIT_FAILURE);
			}
			pthread_cond_wait(&sync_cond,&sync_lock);
			pthread_mutex_unlock(&sync_lock);
	}

    char action;
    long num;

	// make task queue first.
    while (fscanf(fin, "%c %ld\n", &action, &num) == 2) {
        if (action == 'p') {            // process
			pthread_mutex_lock(&lock);

			// inqueue() to the queue
			inqueue(&task_head, num);
			pthread_cond_signal(&arrive_cond);
			pthread_mutex_unlock(&lock);
        } else if (action == 'w') {     // wait
            sleep(num);
        } else {
            printf("ERROR: Unrecognized action: '%c'\n", action);
            exit(EXIT_FAILURE);
        }
    }
    fclose(fin);
	

	/* 
	 * master thread work done.
	 * Broadcast signal to finish job. 
	 * Before broadcasting we should consider whether all of the worker thread work done or not.
	 */
	while(1){
		if(idle_thread_cnt == num_threads && is_queue_empty(&task_head)){
			pthread_mutex_lock(&lock);
			done = true;
			pthread_cond_broadcast(&arrive_cond);
			pthread_mutex_unlock(&lock);
			break;
		}
	}

	for (i = 0; i < num_threads ; i++){
		pthread_join(worker_threads[i], NULL);
	}

    // print results
    printf("%ld %ld %ld %ld\n", sum, odd, min, max);

    // clean up and return
	cleanup();
    return (EXIT_SUCCESS);
}


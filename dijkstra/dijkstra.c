//Inclue header file==============================================
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
//================================================================

//Define struct===================================================
//Link---------------------------
typedef struct _link
{
    int link_id;
    int from_node_id;
    int to_node_id;
    float FFT;
    float capacity;
    int next_inlink_id;
    int next_outlink_id;
} Link_t;

//Node---------------------------
typedef struct _node
{
    int node_id;
    int first_inlink_id;
    int first_outlink_id;
    int heap_position; //heap
    float spl;         //Shortest path length
    int search_flag;   //Use in dijkstra
} Node_t;

//Function prottype===============================================
//-------------------------------------------
int input_parameter_data(char *folder_path_, int *Num_Link_, int *Num_Node_);
int input_netwok_data(char *folder_path_, Link_t *Link_table_, int *Num_Link_, Node_t *Node_table_, int *Num_Node_);

//-------------------------------------------
int dijkstra(Link_t *Link_table_, int *Num_Link_, Node_t *Node_table_, int *Num_Node_);//, int *heap_tree_, int* last_position_);
int pickup_heap(int *heap_tree_, int *last_position_, Node_t *Node_table_);
int add_heap(int *heap_tree_, int *last_position_, Node_t *Node_table_, int add_node_id_);
int up_heap(int *heap_tree_, Node_t *Node_table_, int tr_position_);
int down_heap(int *heap_tree_, int *last_position_, Node_t *Node_table_, int tr_position_);

//main============================================================
int main()
{
    //Setting folder path=========================================
    char *folder_path = malloc(50 * sizeof(char));
    printf("Input folder path\n");
    scanf("%s", folder_path);

    // Input data=================================================
    // Define parameter-------------------------
    int *Num_Node = malloc(sizeof(int));
    int *Num_Link = malloc(sizeof(int));

    //Input parameter data-----------------------
    input_parameter_data(folder_path, Num_Link, Num_Node);

    printf("dbg1\n");
    //Construct Node_t array and Link_t array-----
    int node_table_len = *Num_Node + 1;
    int link_table_len = *Num_Link + 1;
    Node_t *Node_table = malloc(node_table_len * sizeof(Node_t));
    Link_t *Link_table = malloc(link_table_len * sizeof(Link_t));

    //Input network data------------------------
    printf("dbg2\n");
    input_netwok_data(folder_path, Link_table, Num_Link, Node_table, Num_Node);
    printf("dbg3\n");

    printf("Link:%i\n", *Num_Link);
    printf("Node:%i\n", *Num_Node);


    dijkstra(Link_table, Num_Link, Node_table, Num_Node);

    for (int tr_link_id = 1; tr_link_id <= *Num_Link; ++tr_link_id)
    {
        printf("link-id:%i\tcapa:%f\n", tr_link_id, Link_table[tr_link_id].capacity);
    }

    for (int tr_node_id = 0; tr_node_id < *Num_Node; ++tr_node_id)
    {
        printf("node-id:%i\tspl:%f\n", tr_node_id, Node_table[tr_node_id].spl);
    }

    free(Num_Link);
    free(Num_Node);
    free(Link_table);
    free(Node_table);
    return 0;
}

//dijkstra=======================================================
int dijkstra(Link_t *Link_table_, int *Num_Link_, Node_t *Node_table_, int *Num_Node_)
{
    printf("dbg4\n");
    int pivot_node_id = 0;
    Node_table_[pivot_node_id].spl = 0;


    int *heap_tree = malloc((*Num_Node_ + 1) * sizeof(int));
    heap_tree[1] = pivot_node_id; //Set Origin Node ID
    
    int *last_position = malloc(sizeof(int));
    *last_position = 1;
    


    while (1)
    {
        printf("==================================");
        pivot_node_id = pickup_heap(heap_tree, last_position, Node_table_);

        printf("pivot_node_id:%i\n", pivot_node_id);
        printf("last_position:%i\n", *last_position);

        int tr_link_id = Node_table_[pivot_node_id].first_outlink_id;
        while (tr_link_id != -1)
        {
            printf("tr_link_id:%i\n", tr_link_id);
            int tr_node_id = Link_table_[tr_link_id].to_node_id;
            if (Node_table_[tr_node_id].search_flag != 2)
            {
                Node_table_[tr_node_id].search_flag = 1;
                if (Node_table_[tr_node_id].spl > Node_table_[pivot_node_id].spl + Link_table_[tr_link_id].FFT)
                {
                    Node_table_[tr_node_id].spl = Node_table_[pivot_node_id].spl + Link_table_[tr_link_id].FFT;
                    printf("Add node %i\n", tr_node_id);
                    add_heap(heap_tree, last_position, Node_table_, tr_node_id);
                }
            }
            tr_link_id = Link_table_[tr_link_id].next_outlink_id;
        }
        //--------------------------

        Node_table_[pivot_node_id].search_flag = 2;

        if (*last_position == 0)
        {
            break;
        }
    }
    free(last_position);
    free(heap_tree);
    return 0;
}

//==============================================================
int add_heap(int *heap_tree_, int *last_position_, Node_t *Node_table_, int add_node_id_)
{
    if (Node_table_[add_node_id_].heap_position == -1)
    {
        //new node
        int add_position = *last_position_ + 1;

        heap_tree_[add_position] = add_node_id_;
        Node_table_[add_node_id_].heap_position = add_position;
        *last_position_ = add_position;
        
        printf("last_posi:%i\n", *last_position_);
        getchar();

        up_heap(heap_tree_, Node_table_, add_position);
    }
    else
    {
        //already exist
        int add_position = Node_table_[add_node_id_].heap_position;
        up_heap(heap_tree_, Node_table_, add_position);
    }
    return 0;
}

//===============================================================
int pickup_heap(int *heap_tree_, int *last_position_, Node_t *Node_table_)
{
    int top_node_id = heap_tree_[1];
    int last_position = *last_position_;

    heap_tree_[1] = heap_tree_[last_position];
    Node_table_[heap_tree_[last_position]].heap_position = 1;

    *last_position_ = last_position - 1;

    down_heap(heap_tree_, last_position_, Node_table_, 1);
    return top_node_id;
}

//===============================================================
int down_heap(int *heap_tree_, int *last_position_, Node_t *Node_table_, int tr_position_)
{
    int tr_position = tr_position_;
    int small_child_position;

    while (tr_position * 2 <= *last_position_)
    {
        int child1_position = tr_position * 2;
        int child2_position = tr_position * 2 + 1;

        if (child2_position > *last_position_)
        {
            small_child_position = child1_position;
        }
        else if (Node_table_[heap_tree_[child1_position]].spl < Node_table_[heap_tree_[child2_position]].spl)
        {
            small_child_position = child1_position;
        }
        else
        {
            small_child_position = child2_position;
        }

        if (Node_table_[heap_tree_[tr_position]].spl > Node_table_[heap_tree_[small_child_position]].spl)
        {
            int tmp_smallest_node_id = heap_tree_[small_child_position];

            heap_tree_[small_child_position] = heap_tree_[tr_position];
            Node_table_[heap_tree_[tr_position]].heap_position = small_child_position;

            heap_tree_[tr_position] = tmp_smallest_node_id;
            Node_table_[tmp_smallest_node_id].heap_position = tr_position;

            tr_position = small_child_position;
        }
        else
        {
            //---
            break;
        }
    }
    return 0;
}

//================================================================
int up_heap(int *heap_tree_, Node_t *Node_table_, int tr_position_)
{
    int tr_position = tr_position_;
    int parent_position;

    while (tr_position > 1)
    {
        parent_position = tr_position / 2;
        if (Node_table_[heap_tree_[parent_position]].spl < Node_table_[heap_tree_[tr_position]].spl)
        {
            int tmp_smallest_node_id = heap_tree_[parent_position];

            heap_tree_[parent_position] = heap_tree_[tr_position];
            Node_table_[heap_tree_[tr_position]].heap_position = parent_position;

            heap_tree_[tr_position] = tmp_smallest_node_id;
            Node_table_[tmp_smallest_node_id].heap_position = tr_position;

            tr_position = parent_position;
        }
        else
        {
            break;
        }
    }
    return 0;
}

// Input paremeter data function===================================
int input_parameter_data(char *folder_path_, int *Num_Link_, int *Num_Node_)
{
    char tmp_folder_path[50];
    strcpy(tmp_folder_path, folder_path_);
    FILE *input_file = fopen(strcat(tmp_folder_path, "\\parameter.dat"), "r");
    if (input_file == NULL)
    {
        fprintf(stderr, "File Open Failed\n");
        exit(1);
    }
    // int node;
    fscanf(input_file, "%i%i\n", Num_Node_, Num_Link_);
    fclose(input_file);
    return 0;
}
//--------------------------------------------

// Input networks data function====================================
int input_netwok_data(char *folder_path_, Link_t *Link_table_, int *Num_Link_, Node_t *Node_table_, int *Num_Node_)
{
    //Init Node_table----------------------------
    for (int i = 0; i <= *Num_Node_; ++i)
    {
        Node_table_[i].node_id = i;
        Node_table_[i].first_inlink_id = -1;
        Node_table_[i].first_outlink_id = -1;
        Node_table_[i].heap_position = -1; //heap
        Node_table_[i].spl = 999999;         //Shortest path length
        Node_table_[i].search_flag = 0;   //Use in dijkstra
    }

    //------------------------------------------
    char tmp_folder_path[50];
    strcpy(tmp_folder_path, folder_path_);
    FILE *input_file = fopen(strcat(tmp_folder_path, "\\network.dat"), "r");
    if (input_file == NULL)
    {
        fprintf(stderr, "File Open Failed\n");
        exit(1);
    }
    int Link_id = 1;
    while (Link_id <= *Num_Link_)
    {
        int from_node_id;
        int to_node_id;
        float FFT;
        float capacity;
        fscanf(input_file, "%i%i%f%f\n", &from_node_id, &to_node_id, &FFT, &capacity);
        Link_table_[Link_id].from_node_id = from_node_id;
        Link_table_[Link_id].to_node_id = to_node_id;
        Link_table_[Link_id].FFT = FFT;
        Link_table_[Link_id].capacity = capacity;
        Link_table_[Link_id].next_inlink_id = -1;
        Link_table_[Link_id].next_outlink_id = -1;
        //Add inlink and outlink information to Node_table-------------
        //from_node--------------------
        Node_table_[from_node_id].node_id = from_node_id;
        if (Node_table_[from_node_id].first_outlink_id == -1)
        {
            Node_table_[from_node_id].first_outlink_id = Link_id;
        }
        else
        {
            int tr_link_id;
            tr_link_id = Node_table_[from_node_id].first_outlink_id;
            while (1)
            {
                if (Link_table_[tr_link_id].next_outlink_id == -1)
                {
                    Link_table_[tr_link_id].next_outlink_id = Link_id;
                    Link_table_[Link_id].next_outlink_id = -1;
                    break;
                }
                tr_link_id = Link_table_[tr_link_id].next_outlink_id;
            }
        }
        //-------------------------------------------
        //to_node------------------------------------
        Node_table_[to_node_id].node_id = to_node_id;
        if (Node_table_[to_node_id].first_inlink_id == -1)
        {
            Node_table_[to_node_id].first_inlink_id = Link_id;
            Link_table_[Link_id].next_inlink_id = -1;
        }
        else
        {
            int tr_link_id;
            tr_link_id = Node_table_[to_node_id].first_inlink_id;
            while (1)
            {
                if (Link_table_[tr_link_id].next_inlink_id == -1)
                {
                    Link_table_[tr_link_id].next_inlink_id = Link_id;
                    Link_table_[Link_id].next_inlink_id = -1;
                    break;
                }
                tr_link_id = Link_table_[tr_link_id].next_inlink_id;
            }
        }
        //----------------------------------------------
        Link_id = Link_id + 1;
    }
    fclose(input_file);
    return 0;
}
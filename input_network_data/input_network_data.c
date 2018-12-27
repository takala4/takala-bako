//Inclue header file==============================================
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
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
} Node_t;

//Function prottype===============================================
int input_parameter_data(char *folder_path_, int *Num_Node_, int *Num_Link_);
int input_netwok_data(char *folder_path_, Link_t *Link_table_, Node_t *Node_table_, int *Num_Node_);

//main============================================================
int main(){
    //Setting folder path=========================================
    char *folder_path = malloc(50*sizeof(char));
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
    input_netwok_data(folder_path, Link_table, Node_table, Num_Node);
    printf("dbg3\n");

    for(int tr_link_id=1; tr_link_id <= *Num_Link; ++tr_link_id)
    {
        printf("link-id:%i\tcapa:%f\n",tr_link_id, Link_table[tr_link_id].capacity);
    }

    return 0;
}

// Input paremeter data function===================================
int input_parameter_data(char *folder_path_, int *Num_Node_, int *Num_Link_)
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
int input_netwok_data(char *folder_path_, Link_t *Link_table_, Node_t *Node_table_, int *Num_Node_)
{
    //Init Node_table----------------------------
    for (int i = 0; i <= *Num_Node_; ++i)
    {
        Node_table_[i].node_id = i;
        Node_table_[i].first_inlink_id = -1;
        Node_table_[i].first_outlink_id = -1;
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
    while (!feof(input_file))
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
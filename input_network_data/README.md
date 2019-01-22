# Input network data

A C function to read a simple road network data 

# Road Network

* Each node id are integer
* 

# Input Data format

* file name : `network.dat`

|from node id|to node id|free flow travel time|capacity|
|---|---|---|---|
|1  |2  |10 |500|
|1  |3  |13 |300|
|2  |3  |12 |100|
|2  |4  |15 |200|
|3  |4  |17 |600|

![Alt text](https://g.gravizo.com/svg?
  digraph G {
    aize ="4,4";
    main [shape=box];
    main -> parse [weight=8];
    parse -> execute;
    main -> init [style=dotted];
    main -> cleanup;
    execute -> { make_string; printf}
    init -> make_string;
    edge [color=red];
    main -> printf [style=bold,label="100 times"];
    make_string [label="make a string"];
    node [shape=box,style=filled,color=".7 .3 1.0"];
    execute -> compare;
  }
)

# Usage

#


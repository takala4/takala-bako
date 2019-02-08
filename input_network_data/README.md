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

![Alt text](https://g.gravizo.com/svg?%20digraph%20G%20%7B%0A%20%20%20%201%20%3B%0A%20%20%20%202%20%3B%0A%20%20%20%203%20%3B%0A%20%20%20%204%20%3B%0A%20%20%20%201%20-%3E%202%3B%0A%20%20%20%201%20-%3E%203%3B%0A%20%20%20%202%20-%3E%203%3B%0A%20%20%20%202%20-%3E%204%3B%0A%20%20%20%203%20-%3E%204%3B%0A%20%20%20%20%7Brank%20%3D%20same%3B%202%3B%7D%0A%20%20%20%20%7Brank%20%3D%20same%3B%201%3B%204%3B%7D%0A%20%20%20%20%7Brank%20%3D%20same%3B%203%3B%7D%0A%7D)

```graphviz
digraph hierarchy {

                nodesep=1.0 // increases the separation between nodes

                node [color=Red,fontname=Courier,shape=box] //All nodes will this shape and colour
                edge [color=Blue, style=dashed] //All the lines look like this

                Headteacher->{Deputy1 Deputy2 BusinessManager}
                Deputy1->{Teacher1 Teacher2}
                BusinessManager->ITManager
                {rank=same;ITManager Teacher1 Teacher2}  // Put them on the same level
}
```

# Usage

#


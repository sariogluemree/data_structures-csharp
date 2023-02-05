using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DataProje4
{
    public class Node
    {
        public int key, height;
        public Node left, right;

        public Node(int d)
        {
            key = d;
            height = 1;
        }
    }
    public class AVLTree //https://www.geeksforgeeks.org/avl-tree-set-1-insertion/
    {

        public Node root;

        // A utility function to get 
        // the height of the tree  
        public int height(Node N)
        {
            if (N == null)
            {
                return 0;
            }
            return N.height;
        }

        // A utility function to get 
        // maximum of two integers  
        public int max(int a, int b)
        {
            return (a > b) ? a : b;
        }

        // A utility function to right  
        // rotate subtree rooted with y  
        // See the diagram given above.  
        public Node rightRotate(Node y)
        {
            Node x = y.left;
            Node T2 = x.right;

            // Perform rotation  
            x.right = y;
            y.left = T2;

            // Update heights  
            y.height = max(height(y.left),
                        height(y.right)) + 1;
            x.height = max(height(x.left),
                        height(x.right)) + 1;

            // Return new root  
            return x;
        }

        // A utility function to left 
        // rotate subtree rooted with x  
        // See the diagram given above.  
        public Node leftRotate(Node x)
        {
            Node y = x.right;
            Node T2 = y.left;

            // Perform rotation  
            y.left = x;
            x.right = T2;

            // Update heights  
            x.height = max(height(x.left),
                        height(x.right)) + 1;
            y.height = max(height(y.left),
                        height(y.right)) + 1;

            // Return new root  
            return y;
        }

        // Get Balance factor of node N  
        public int getBalance(Node N)
        {
            if (N == null)
                return 0;

            return height(N.left) - height(N.right);
        }

        public Node insert(Node node, int key)
        {

            /* 1. Perform the normal BST insertion */
            if (node == null)
                return (new Node(key));

            if (key < node.key)
                node.left = insert(node.left, key);
            else if (key > node.key)
                node.right = insert(node.right, key);
            else // Duplicate keys not allowed  
                return node;

            /* 2. Update height of this ancestor node */
            node.height = 1 + max(height(node.left),
                                height(node.right));

            /* 3. Get the balance factor of this ancestor  
                node to check whether this node became  
                unbalanced */
            int balance = getBalance(node);

            // If this node becomes unbalanced, then there  
            // are 4 cases Left Left Case  
            if (balance > 1 && key < node.left.key)
                return rightRotate(node);

            // Right Right Case  
            if (balance < -1 && key > node.right.key)
                return leftRotate(node);

            // Left Right Case  
            if (balance > 1 && key > node.left.key)
            {
                node.left = leftRotate(node.left);
                return rightRotate(node);
            }

            // Right Left Case  
            if (balance < -1 && key < node.right.key)
            {
                node.right = rightRotate(node.right);
                return leftRotate(node);
            }

            /* return the (unchanged) node pointer */
            return node;
        }
        public void PreOrder(Node node)
        {
            if (node != null)
            {
                Console.Write(node.key + " ");
                PreOrder(node.left);
                PreOrder(node.right);
            }
        }
    }

    public class MST
    {
        // Number of vertices in the graph
        static int V = 5;

        // A utility function to find
        // the vertex with minimum key
        // value, from the set of vertices
        // not yet included in MST
        public static int minKey(int[] key, bool[] mstSet)
        {
            // Initialize min value
            int min = int.MaxValue, min_index = -1;

            for (int v = 0; v < V; v++)
                if (mstSet[v] == false && key[v] < min)
                {
                    min = key[v];
                    min_index = v;
                }
            return min_index;
        }

        // A utility function to print
        // the constructed MST stored in
        // parent[]
        public static void printMST(int[] parent, int[,] graph)
        {
            Console.WriteLine("Edge \tWeight");
            for (int i = 1; i < V; i++)
                Console.WriteLine(parent[i] + " - " + i + "\t" + graph[i, parent[i]]);
        }

        // Function to construct and
        // print MST for a graph represented
        // using adjacency matrix representation
        public static void primMST(int[,] graph)
        {
            // Array to store constructed MST
            int[] parent = new int[V];

            // Key values used to pick minimum weight edge in cut
            int[] key = new int[V];

            // To represent set of vertices included in MST
            bool[] mstSet = new bool[V];

            // Initialize all keys as INFINITE
            for (int i = 0; i < V; i++)
            {
                key[i] = int.MaxValue;
                mstSet[i] = false;
            }

            // Always include first 1st vertex in MST.
            // Make key 0 so that this vertex is
            // picked as first vertex
            // First node is always root of MST
            key[0] = 0;
            parent[0] = -1;

            // The MST will have V vertices
            for (int count = 0; count < V - 1; count++)
            {
                // Pick thd minimum key vertex
                // from the set of vertices
                // not yet included in MST
                int u = minKey(key, mstSet);

                // Add the picked vertex
                // to the MST Set
                mstSet[u] = true;

                // Update key value and parent
                // index of the adjacent vertices
                // of the picked vertex. Consider
                // only those vertices which are
                // not yet included in MST
                for (int v = 0; v < V; v++)

                    // graph[u][v] is non zero only
                    // for adjacent vertices of m
                    // mstSet[v] is false for vertices
                    // not yet included in MST Update
                    // the key only if graph[u][v] is
                    // smaller than key[v]
                    if (graph[u, v] != 0 && mstSet[v] == false
                        && graph[u, v] < key[v])
                    {
                        parent[v] = u;
                        key[v] = graph[u, v];
                    }
            }
            printMST(parent, graph);  // print the constructed MST
        }

        class StackX
        {
            private const int SIZE = 20;
            private int[] st;
            private int top;
            public StackX()   // constructor
            {
                st = new int[SIZE];   // make array
                top = -1;
            }
            public void push(int j)   // put item on stack
            { st[++top] = j; }
            public int pop()   // take item off stack
            { return st[top--]; }

            public int peek()   // peek at top of stack
            {
                return st[top];
            }
            public Boolean isEmpty()  // true if nothing on stack
            { return (top == -1); }
        }

        class Vertex
        {
            public char label;   // label (e.g. ‘A’)
            public Boolean wasVisited;
            public Vertex(char lab)   // constructor
            { label = lab; wasVisited = false; }
        }

        class Graph
        {
            private const int MAX_VERTS = 20;
            private Vertex[] vertexList;   // list of vertices
            private int[,] adjMat;   // adjacency matrix
            private int nVerts;    // current number of vertices
            private StackX theStack;

            public Graph()   // constructor
            {
                vertexList = new Vertex[MAX_VERTS];
                adjMat = new int[MAX_VERTS, MAX_VERTS];
                nVerts = 0; 
                for (int j = 0; j < MAX_VERTS; j++)    // set adjacency
                    for (int k = 0; k < MAX_VERTS; k++)    // matrix to 0
                        adjMat[j, k] = 0; theStack = new StackX();
            }
            public void addVertex(char lab)
            {
                vertexList[nVerts++] = new Vertex(lab);
            }
            public void addEdge(int start, int end)
            {
                adjMat[start, end] = 1;
                adjMat[end, start] = 1;
            }

            public void displayVertex(int v)
            {
                Console.WriteLine(vertexList[v].label);
            }
            public void dfs()  // depth-first search
            {                                           // begin at vertex 0
                vertexList[0].wasVisited = true;        // mark it
                displayVertex(0);                       // display it
                theStack.push(0);                       // push it
                while (!theStack.isEmpty())    // until stack empty
                {
                    // get an unvisited vertex adjacent to stack top
                    int v = getAdjUnvisitedVertex(theStack.peek());
                    if (v == -1)
                        theStack.pop();
                    else
                    {
                        vertexList[v].wasVisited = true;   // mark it
                        displayVertex(v);   // display it
                        theStack.push(v);   // push it
                    }
                }
                for (int j = 0; j < nVerts; j++)   // reset flags
                    vertexList[j].wasVisited = false;
            }
            // returns an unvisited vertex adjacent to v
            public int getAdjUnvisitedVertex(int v)
            {
                for (int j = 0; j < nVerts; j++)
                    if (adjMat[v, j] == 1 && vertexList[j].wasVisited == false)
                        return j;  // return first such vertex
                return -1;   // no such vertices
            }
        }

        class Program
        {
            public static int INFINITY = 1000;
            static void Main(string[] args)
            {
                AVLTree tree = new AVLTree();  // Create AVL tree

                /* Constructing tree given in the above figure */
                tree.root = tree.insert(tree.root, 10);
                tree.root = tree.insert(tree.root, 20);
                tree.root = tree.insert(tree.root, 30);
                tree.root = tree.insert(tree.root, 40);
                tree.root = tree.insert(tree.root, 50);
                tree.root = tree.insert(tree.root, 25);

                Console.Write("Preorder traversal" + " of constructed tree is : ");
                tree.PreOrder(tree.root);
                Console.ReadKey();
                Console.WriteLine();

                //b.Prim’s MST (Minimum Spanning Tree) 
                int[,] graph = new int[,] { { 0, 2, 0, 6, 0 },
                                      { 2, 0, 3, 8, 5 },
                                      { 0, 3, 0, 0, 7 },
                                      { 6, 8, 0, 0, 9 },
                                      { 0, 5, 7, 9, 0 } };
                primMST(graph);  // Print the solution
                Console.ReadKey();
                Console.WriteLine();

                //c.DFT (Depth-First Traverse)
                Graph theGraph = new Graph();
                theGraph.addVertex('A');    // 0 (start for dfs)
                theGraph.addVertex('B');    // 1
                theGraph.addVertex('C');    // 2
                theGraph.addVertex('D');    // 3
                theGraph.addVertex('E');    // 4

                theGraph.addEdge(0, 1);     // AB
                theGraph.addEdge(1, 2);     // BC
                theGraph.addEdge(0, 3);     // AD
                theGraph.addEdge(3, 4);     // DE
                Console.WriteLine("Visits: ");
                theGraph.dfs();    // depth-first search
                Console.WriteLine(" ");
                Console.ReadKey();


                //a.Dijkstra’s Shortest Path 
                int N = 5;
                int SRC = 0;

                // int[][] cost	= new int[N][N];
                int[,] cost = {
                               { INFINITY,    5,    3, INFINITY,    2},
                               { INFINITY, INFINITY,    2,    6, INFINITY},
                               { INFINITY,    1, INFINITY,    2, INFINITY},
                               { INFINITY, INFINITY, INFINITY, INFINITY, INFINITY},
                               { INFINITY,    6,   10,    4,    INFINITY}};

                int[] distances = new int[N];

                int[] previous = Distance(N, cost, distances, SRC);

                Console.WriteLine("0 numaralı düğümden tüm düğümlere en kısa yol uzunlukları ");
                Console.WriteLine("----------------------------------------------------------");
                for (int i = 0; i < distances.Length; ++i)
                    if (distances[i] != INFINITY)
                        Console.WriteLine(distances[i]);
                    else Console.WriteLine("INFINITY");

                int DEST = 1;
                Console.WriteLine("\n Shortest path from " + SRC + " to " + DEST + " (straight):");
                printShortestPathStraight(DEST, previous);
                Console.WriteLine("\n\n Shortest path from " + SRC + " to " + DEST + " (reverse) :");
                printShortestPathReverse(DEST, previous);
                Console.ReadLine();
            }

            public static int[] Distance(int N, int[,] cost, int[] D, int src)
            {

                int w, v, min;

                bool[] visited = new bool[N];

                int[] previous = new int[N]; //for tracking shortest paths (güzergah)

                //initialization of D[], visited[] and previous[] arrays according to src node
                for (v = 0; v < N; v++)
                {
                    if (v != src)
                    {
                        visited[v] = false;
                        D[v] = cost[src, v];
                        if (D[v] != INFINITY) //there is a connection between src and v
                        {
                            previous[v] = src;
                        }
                        else //no path from source
                        {
                            previous[v] = -1;
                        }
                    }
                    else
                    {
                        visited[v] = true;
                        D[v] = 0;
                        previous[v] = -1;
                    }

                }

                // Searching for shortest paths
                for (int i = 0; i < N; ++i)
                {
                    min = INFINITY;
                    for (w = 0; w < N; w++)
                        if (!visited[w])
                            if (D[w] < min)
                            {
                                v = w;
                                min = D[w];
                            }

                    visited[v] = true;

                    for (w = 0; w < N; w++)
                        if (!visited[w])
                            if (min + cost[v, w] < D[w])
                            {
                                D[w] = min + cost[v, w];
                                previous[w] = v; //update the path info
                            }
                }

                return previous;
            }

            //  Print Shortest Path Straight
            public static void printShortestPathStraight(int dest, int[] previous)
            {
                Stack<int> pathStack = new Stack<int>();

                int current = dest;
                pathStack.Push(current);

                while (previous[current] != -1)
                {
                    current = previous[current];
                    pathStack.Push(current);
                }

                if (pathStack.Count == 1)
                {
                    Console.Write(" NO PATH");
                    return;
                }

                while (pathStack.Count > 0)
                {
                    Console.Write(" -> " + pathStack.Pop());
                }
            }

            // Print Shortest Path Reverse
            public static void printShortestPathReverse(int dest, int[] previous)
            {
                int current = dest;
                Console.Write(dest + " <- ");

                while (previous[current] != -1)
                {
                    current = previous[current];
                    Console.Write(current + " <- ");
                }

            }
        }
    }
}


/*
This code should be run With the following command line arguments:
programName txtWithAllTheNodesInTheNetwork txtWithTheScoresForAllNodes txtWithTheNetworkConnections txtWithUpstreamInput txtWithDownstreamInput NameOfTheOutputFile
txtWithAllTheNodesInTheNetwork should be formatted as one node name per line
txtWithTheScoresForAllNodes should be formatted as NodeNAme \t NodeScore per line, the score in line n will be used for node in line n
txtWithTheNetworkConnections should have 2 columns separated by a \t, the edge will come from the first column node an go to the second column node
txtWithUpstreamInput & txtWithDownstreamInput should be just one column with the nodes that should be considered sources/sinks (up=source/down=sink)
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <limits.h>
#include <string.h>
#include <queue>
#include <iomanip>

using namespace std;
//ofstream myfile ("/home/user/Desktop/AT Progs/NetTestResults/TestResultsNetFUNKY");

/* Returns true if there is a path from source 's' to sink 't' in
residual graph. Also fills parent[] to store the path */
bool bfs(int **rGraph, int s, int t, int parent[],int V)
{
    // Create a visited array and mark all vertices as not visited
    bool visited[V];
    memset(visited, 0, sizeof(visited));

    // Create a queue, enqueue source vertex and mark source vertex
    // as visited
    queue <int> q;
    q.push(s);
    visited[s] = true;
    parent[s] = -1;

    // Standard BFS (Breadth First Traversal) Loop
    while (!q.empty())
    {
        int u = q.front();
        q.pop();

        for (int v=0; v<V; v++)
        {
            if (visited[v]==false && *(*(rGraph + u) + v) > 0)
            {
                q.push(v);
                parent[v] = u;
                visited[v] = true;
            }
        }
    }

    // If we reached sink in BFS starting from source, then return
    // true, else false
    return (visited[t] == true);
}

// Returns the maximum flow from s to t in the given graph
int fordFulkersonModified(int **graph, int **fGraph, int s, int t,int V)
{
    int u, v;

    // Create a residual graph and fill the residual graph with
    // given capacities in the original graph as residual capacities
    // in residual graph
    int **rGraph = new int*[V];
    for (size_t indice = 0; indice < V; indice++) {
      rGraph[indice] = new int[V];
    }
    for (int indice=0; indice<V; indice++)
    {
        for(int j=0; j<V; j++)
        {
            *(*(rGraph + indice) + j)=0;
        }

    }
    // Residual graph where rGraph[i][j] indicates
    // residual capacity of edge from i to j (if there
    // is an edge. If rGraph[i][j] is 0, then there is not)
    for (u = 0; u < V; u++)
        for (v = 0; v < V; v++)
             *(*(rGraph + u) + v) = *(*(graph + u) + v);

    int parent[V]; // This array is filled by BFS to store path

    int max_flow = 0; // There is no flow initially

    // Augment the flow while there is path from source to sink
    while (bfs(rGraph, s, t, parent,V)) // while there is a path from s to t
    {
        // Find minimum residual capacity of the edges along the
        // path filled by BFS. Or we can say find the maximum flow
        // through the path found.
        int path_flow = INT_MAX;
        for (v=t; v!=s; v=parent[v]) // starts at sink and walks the path backward until the source
        {
            u = parent[v];
            path_flow = min(path_flow, rGraph[u][v]);
        }

        // update residual capacities of the edges and reverse edges
        // along the path
        for (v=t; v != s; v=parent[v])
        {
            u = parent[v];
            *(*(rGraph + v) + u) += path_flow;
            *(*(rGraph + u) + v) -= path_flow;

            *(*(fGraph + u) + v) += path_flow;// flow graph fGraph[i][j] = all the flow going from node i to node j
            if( *(*(graph + u) + v)<path_flow  )//if the graph edge is smaller than the flow, than the flow came through a residual edge
            {
                *(*(fGraph + u) + v)-= path_flow; // in case it is a backward/residual edge, removes what was just added,
                *(*(fGraph + v) + u)-= path_flow; //and removes flow from v to u, because it was redirected
            }


        }

        // Add path flow to overall flow
        max_flow += path_flow;
    }
    for (size_t indice = V; indice > 0; )
    {
      delete[] rGraph[--indice];
    }
    delete[]rGraph;

    return max_flow;// Return the overall flow
}

void GlobalSourcesAndSinks(int **graph, int sourceArray[], int sinkArray[],int V)
{
    sourceArray[0] = 0;
    sinkArray[V] = 0;
    for(int i=0; i<V; i++)
    {
        if(i==0)
        {
            cout << "node 0 is the global source"<< endl;
            //myfile<< "node 0 is the global source"<< endl;
        }
        else if (i==V-1)
        {
            cout << "node "<< V-1 << " is the global sink"<< endl;
            //myfile<< "node "<< V-1 << " is the global sink"<< endl;
        }
        else
        {
            if(sourceArray[i]==0 && sinkArray[i]==0)
            {
                sourceArray[i]=1;
                sinkArray[i]=1;
                cout<< "node "<< i << " is a isolated node"<< endl;
                //myfile<< "node "<< i << " is a isolated node"<< endl;
            }
            if (sourceArray[i] == 0)
            {
                *(*(graph + 0) + i) = INT_MAX;
                cout << "node "<< i << " is a source"<< endl;
                //myfile << "node "<< i << " is a source"<< endl;
            }
            if (sinkArray[i] == 0)
            {
                *(*(graph + i) + V-1) = INT_MAX;
                cout << "node "<< i << " is a sink"<< endl;
                //myfile << "node "<< i << " is a sink"<< endl;
            }

        }

    }
    cout << endl;

}

void FlowThroughEachNode (int **fGraph, int fArray[],int V)
{

    for (int i=0; i<V; i++)
    {
        fArray[i]=0;
    }

    for (int i = 0; i<V; i++)
    {
        for(int j = 0; j<V; j++)
            fArray[i] = fArray[i]+*(*(fGraph + i) + j);
        if(fArray[i]!=0)
            cout << "The flow coming through node " << i << " is: " <<fArray[i]<<endl;
    }
    cout << endl;
}

void ScoreEachNode(int **graph, int fArray[], float scoreArray[], int sinkArray[], int sourceArray[], int ScoreNumber[],int V,string NameNumber[], ofstream &myfile)//ScoreNumber[V-2], string NameNumber V-2
{
    float aux;
    int Aux[15];
    for(int i=0; i<15; i++)
    {
        Aux[i]=0;
    }
    for (int i=0; i<V; i++)
    {
        aux=0;
        if(sinkArray[i]!=0 && sourceArray[i]!=0)
        {
            for(int j=1; j<V-1; j++)
            {
                if(*(*(graph + i) + j) != 0)
                    aux = aux+ (float)ScoreNumber[j-1]/1000;
            }
            scoreArray[i]=fArray[i] * (1+aux);
        }
        else
            scoreArray[i] = 0; //Zero's the scor`e for sources and sinks;
        for(int j=0; j<15; j++)
        {
            if (scoreArray[i] > scoreArray[Aux[j]])
            {
                for(int k=14;k>j;k--)
                {
                    Aux[k]=Aux[k-1];
                }
                Aux[j]=i;
                break;
            }
        }
    }
    cout << "The 10 top scoring nodes are:";
    for(int i=0;i<15;i++)
    {
        cout << "\nNode "<< Aux[i] << " \t("<< NameNumber[Aux[i]-1] << ")\tScore: " << scoreArray[Aux[i]] << "\tFlow: " << fArray[Aux[i]];
        //myfile << "\nNode "<< Aux[i] << " \t("<< NameNumber[Aux[i]-1] << ")\tScore: " << scoreArray[Aux[i]] << "\tFlow: " << fArray[Aux[i]];
    }

}

int NumberOFNodesOnMainPath(int fArray[], int V)
{
    int NumberOfNodes = 2;
    for (int i=1;i<V-1;i++)
    {
        if (fArray[i] > 0)
        {
            NumberOfNodes ++;
        }
    }
    return NumberOfNodes;
}
void NodesOnMainPath(int fArray[], string NameNumber[],  string NodesOnReducedPath[], int V, ofstream &myfile)
{
    int j=0;
    for (int i=1;i<V-1;i++)
    {
        if (fArray[i] > 0)
        {
            myfile<< NameNumber[i-1]<< endl;
            NodesOnReducedPath[j]=NameNumber[i-1];
            j++;
        }
    }
}

void Disruption(int i,int **graph,int **fGraph, int s, int t, int fArray[], int initialFlow, string NameNumber[],int V, ofstream &myfile)//string NameNumber[V-2]
{
    int **auxGraph = new int*[V];
    for (size_t indice = 0; indice < V; indice++) {
      auxGraph[indice] = new int[V];
    }
    for (int indice=0; indice<V; indice++)
    {
        for(int j=0; j<V; j++)
        {
            *(*(auxGraph + indice) + j)=0;
        }

    }
    for (int j=0; j<V; j++)
    {
        //copies the graph to auxGraph, except for node i, which will be removed.
        if(j!=i)
            for(int k=0; k<V; k++)
            {
                *(*(auxGraph + j) + k)=*(*(graph + j) + k);
            }
        else
            for(int k=0; k<V; k++)
            {
                *(*(auxGraph + i) + k)=0;
            }

        //Re-initializes the flowGraph
        for(int k=0; k<V; k++)
        {
            *(*(fGraph + j) + k)=0;
        }
    }
    float flowWithoutNodeI = fordFulkersonModified(auxGraph,fGraph, 0, V-1,V);
    float disruption = 100-100*flowWithoutNodeI/initialFlow;
    float AbsDisrupt = initialFlow-flowWithoutNodeI;
    cout<< fArray[i] << "\t|\t" << fixed << setprecision(2) <<flowWithoutNodeI
    << "\t|\t" <<disruption << "%\t|\t" << AbsDisrupt<< "\t|\t" << NameNumber[i-1]<<endl;

    myfile<< NameNumber[i-1] << "\t" << fArray[i] <<"\t"<< fixed << setprecision(2) << disruption<< "\t" <<AbsDisrupt <<endl;
    for (size_t indice = V; indice > 0; )
    {
      delete[] auxGraph[--indice];
    }
    delete[]auxGraph;

}


// Driver program to test above functions
int main(int argc, char *argv[])
{

    int V=2;
    string aux;
    int aux2;
    if ( argc != 7 ) // argc should be 7 considering a correct call
    {
        cout<<"incorrect call, too many or too few arguments";
        return 0;
    }
    ifstream InFile (argv[1]);
    ifstream InFile1 (argv[1]); //AllNodes
    ifstream InFile2(argv[2]); //AllScores
    ifstream InFile22(argv[2]);
    ifstream InFile3(argv[3]); //AllConnections
    ifstream InFile4(argv[3]);
    ifstream InFileUps(argv[4]);
    ifstream InFileDown(argv[5]);
    ifstream InFileUps2(argv[4]);
    ifstream InFileDown2(argv[5]);
    ofstream myfile (argv[6]);
    for(int k=0 ; !InFile.eof(); k++) //finds the size of the file
    {
        getline (InFile, aux);
        V++;
    }
    string NameNumber[V-2];
    for(int k=0 ; !InFile1.eof(); k++) //reads the file, and stores which name corresponds to each node
    {
        getline (InFile1, NameNumber[k]);
    }


    int ScoreNumber[V-2];
    for(int i =0; i<V-2; i++)
    {
        ScoreNumber[i]=50;
    }
    for(int i=0 ; !InFile2.eof(); i++)//reads the file, and stores the score of each node
    {
        getline (InFile2, aux,'\t');
        for (int j=0; j<V-2; j++)
        {
            if(aux==NameNumber[j])
            {
                aux2=j;
                break;
            }
        }
        getline (InFile2, aux);
        ScoreNumber[aux2] = ScoreNumber[aux2]+atoi(aux.c_str());
    }


    int sourceArray[V];
    //sourceArray[i] = 0 -> mean that node i is a source
    for (int i=0;i<V;i++)
    {
        sourceArray[i]=1;
    }
    for (int i=0;!InFileUps.eof();i++)
    {
        getline (InFileUps, aux);
        for(int j=0;j<V-2;j++)
        {
            if(aux==NameNumber[j])
            {
                sourceArray[j+1]=0;
            }
        }
    }

    int sinkArray[V];
    //sinkArray[i] = 0 -> mean that node i is a sink
    for (int i=0;i<V;i++)
    {
        sinkArray[i]=1;
    }
    for (int i=0;!InFileDown.eof();i++)
    {
        getline (InFileDown, aux);
        for(int j=0;j<V-2;j++)
        {
            if(aux==NameNumber[j])
            {
                sinkArray[j+1]=0;
            }
        }
    }






    int **graph = new int*[V];
    int **fGraph= new int*[V];//*(*(graph + i) + j)
    for (size_t i = 0; i < V; i++) {
      graph[i] = new int[V];
      fGraph[i]= new int[V];
    }
    for (int i=0; i<V; i++)
    {
        for(int j=0; j<V; j++)
        {
            *(*(graph + i) + j)=0;
            *(*(fGraph + i) + j)=0;
        }

    }



    bool NameFound=false;
    for(int i=0; !InFile3.eof(); i++)
    {
        getline(InFile3, aux,'\t');
        for(int j=0; j<V-2; j++)
        {
            NameFound=false;
            if(NameNumber[j]==aux)
            {
                aux2=j;
                NameFound=true;
                break;
            }
        }
        if(NameFound==false)
        {
            cout<<  aux << " 1is present in the Connections file, but is not present in the NodeNames"<< endl;
        }
        getline(InFile3, aux);
        for(int j=0; j<V-2; j++)
        {
            //the graph is (n+2)x(n+2) where n is the number of nodes in the network node 0 is the global source and node n+1 = is the global sink
            NameFound=false;
            if(NameNumber[j]==aux)
            {
                *(*(graph + aux2+1) + j+1)= (ScoreNumber[aux2]+ScoreNumber[j])/2;//ScoreNumber[aux2]=NonAvrgedScores
                //(ScoreNumber[aux2]+ScoreNumber[j])/2 = AveragedScores

                NameFound=true;
                break;
            }//the +0.5 is in order to make any non int number be rounded to and int number instead of just ignoring the non integer part

        }
        if(NameFound==false)
        {
            cout<<  aux << " 2is present in the Connections file, but is not present in the NodeNames"<<NameNumber[aux2]<< endl;
        }
    }



    int fArray[V];
    //Flow array fArray[n] = all the flow coming from node n
    //float scoreArray[V];
    //scoreArray[i] = the score of node i+1 ## the size is V-2 because the global sink and global source wont be assigned scores

    GlobalSourcesAndSinks(graph, sourceArray, sinkArray, V);
    //assigns 0 and 1 to the sourceArray and sinkArray, accordingly to the characteristics of the nodes on graph

    cout << "First Iteration:\n";
    myfile<< "First Iteration:\n";
    int flowFirstIteraction = fordFulkersonModified(graph,fGraph, 0, V-1,V);
    cout << "The maximum possible flow is " << flowFirstIteraction << endl << endl;

    FlowThroughEachNode(fGraph, fArray,V);
    //prints the flow through each node and stores it in the flow array getting the information from the flow graph, provided by the fordFulkersonModified

    //ScoreEachNode(graph, fArray, scoreArray, sinkArray, sourceArray,ScoreNumber,V,NameNumber, myfile);
    //calculates the score of each node (except for sources and sinks, which are assigned a score of 0) and prints the top 3 scoring nodes

    cout << endl<<endl;

    cout << "Second Iteration:\n";
    myfile << "Second Iteration:\n";

    int NumberOfNodes = NumberOFNodesOnMainPath(fArray,V);
    cout<< "1\n";
    string NameNumberReduced[NumberOfNodes-2];
    int ScoreNumberReduced[NumberOfNodes-2];
    NodesOnMainPath(fArray, NameNumber, NameNumberReduced, V, myfile);

    bool Found = false;
    for(int i =0; i<NumberOfNodes-2; i++)
    {
        ScoreNumberReduced[i]=50;
    }
    for(int i=0 ; !InFile22.eof(); i++)//reads the file, and stores the score of each node
    {
        getline (InFile22, aux,'\t');
        for (int j=0; j<NumberOfNodes-2; j++)
        {
            if(aux==NameNumberReduced[j])
            {
                aux2=j;
                getline (InFile22, aux);
                ScoreNumberReduced[aux2] = ScoreNumberReduced[aux2]+atoi(aux.c_str());
                Found=true;
                break;
            }
        }
        if (Found == false)
        {
        getline (InFile22, aux);

        }
        Found=false;
    }


    for (int i=0;i< NumberOfNodes-2;i++)
    {
        cout << NameNumberReduced[i] << " score: "<< ScoreNumberReduced[i]<<endl;
    }
    system("pause");
    for (size_t i = V; i > 0; )
    {
      delete[] fGraph[--i];
      delete[] graph[i];
    }
    delete[]fGraph;
    delete[] graph;

    aux2=0;
    NameFound=false;
    int **rgraph = new int*[NumberOfNodes];
    int **rfGraph= new int*[NumberOfNodes];//*(*(graph + i) + j)
    for (size_t i = 0; i < NumberOfNodes; i++) {
      rgraph[i] = new int[NumberOfNodes];
      rfGraph[i]= new int[NumberOfNodes];
    }
    for (int i=0; i<NumberOfNodes; i++)
    {
        for(int j=0; j<NumberOfNodes; j++)
        {
            *(*(rgraph + i) + j)=0;
            *(*(rfGraph + i) + j)=0;
        }

    }




    for(int i=0; !InFile4.eof(); i++)
    {
        getline(InFile4, aux,'\t');
        for(int j=0; j<NumberOfNodes-2; j++)
        {
            if(NameNumberReduced[j]==aux)
            {
                aux2=j;
                NameFound=true;
                break;
            }
        }
        getline(InFile4, aux);
        if (NameFound==true);
        for(int j=0; j<NumberOfNodes-2; j++)
        {
            //the graph is (n+2)x(n+2) where n is the number of nodes in the network node 0 is the global source and node n+1 = is the global sink
            if(NameNumberReduced[j]==aux)
            {
                *(*(rgraph + aux2+1) + j+1)= (ScoreNumberReduced[aux2]+ScoreNumberReduced[j])/2;//ScoreNumber[aux2]=NonAvrgedScores
                //(ScoreNumber[aux2]+ScoreNumber[j])/2 = AveragedScores

                break;
            }//the +0.5 is in order to make any non int number be rounded to and int number instead of just ignoring the non integer part

        }
        NameFound=false;
        aux2=NumberOfNodes-2;
    }




    int sourceArrayReduced[NumberOfNodes];
    //sourceArray[i] = 0 -> mean that node i is a source
    for (int i=0;i<NumberOfNodes;i++)
    {
        sourceArrayReduced[i]=1;
    }
    for (int i=0;!InFileUps2.eof();i++)
    {
        getline (InFileUps2, aux);
        for(int j=0;j<NumberOfNodes-2;j++)
        {
            if(aux==NameNumberReduced[j])
            {
                sourceArrayReduced[j+1]=0;
            }
        }
    }


    int sinkArrayReduced[NumberOfNodes];
    //sinkArray[i] = 0 -> mean that node i is a sink
    for (int i=0;i<NumberOfNodes;i++)
    {
        sinkArrayReduced[i]=1;
    }

    for (int i=0;!InFileDown2.eof();i++)
    {
        getline (InFileDown2, aux);
        for(int j=0;j<NumberOfNodes-2;j++)
        {
            if(aux==NameNumberReduced[j])
            {
                sinkArrayReduced[j+1]=0;
            }
        }
    }


    GlobalSourcesAndSinks(rgraph, sourceArrayReduced, sinkArrayReduced, NumberOfNodes);
    cout << "first time: " << flowFirstIteraction<<endl;
    flowFirstIteraction = fordFulkersonModified(rgraph,rfGraph, 0, NumberOfNodes-1,NumberOfNodes);
    cout << "The maximum possible flow is " << flowFirstIteraction << endl << endl;

    int fArrayReduced[NumberOfNodes];
    FlowThroughEachNode(rfGraph, fArrayReduced,NumberOfNodes);

    cout<< "For Nodes that have flow going through them, and are neither sources nor sinks:\n";
    cout << "Flow\t|\tMFlow w/o Node\t|\tScore\t|\tDisrupt\t|\tAbsDisrupt\t|\tNode\n";
    //im sending FlowRankArray[i] so that the table will be printed in order by the Flow ranking of each node;
    for (int i=0; i<NumberOfNodes; i++)
    {
        if(fArrayReduced[i]>0 && sourceArrayReduced[i]!=0 && sinkArrayReduced[i]!=0)
        {
            Disruption(i,rgraph,rfGraph,0,NumberOfNodes-1, fArrayReduced, flowFirstIteraction, NameNumberReduced,NumberOfNodes, myfile);
        }
    }


    //cout << "From the first to second iteraction the flow has been reduced by: "<< 100-100*flowSecondIteraction/flowFirstIteraction<<"%\n";
    return 0;

}


//########
// source: http://kgm-par.googlecode.com/svn/trunk/ kgm-par-read-only
// modified: stoero
//########

#include <iostream>
#include <fstream>
#include <limits>
#include <list>

#include "mpi.h"
#include <vector>
#include <string>
#include <stdlib.h>
#include <map>
#include <set>

using namespace std;

//##############################################################################
//#include "Edge.h"
//##############################################################################

class Edge {
public:
    int FromNode;
    int ToNode;
    bool IsVisited;

    Edge(int fromNode, int toNode, bool isVisited);
    int operator==(Edge e);
};

//######################

Edge::Edge(int fromNode, int toNode, bool isVisited) {
    this->FromNode = fromNode;
    this->ToNode = toNode;
    this->IsVisited = isVisited;
}

int Edge::operator==(Edge e) {
    return FromNode==e.FromNode && ToNode==e.ToNode;
}

//##############################################################################
//#include "Graph.h"
//##############################################################################

class Graph {

public:
    int NodesCount;
    vector<vector<int> > AdjacencyList;

    Graph(ifstream& ifs);
    ~Graph();
};

//######################

Graph::Graph(ifstream& ifs)	{
	string line;

	getline(ifs, line);
	NodesCount = atoi(line.c_str());


	int lineLength;
	int i;

	int lineIndex = 0;
		
	while( getline( ifs, line ) )
	{
		lineLength = line.length();
		if(lineLength == NodesCount)
		{
			AdjacencyList.push_back(vector<int>());
			for (i=0; i < lineLength; i++)
			{
				char item = line[i];
				int number = atoi(&item);
				if(number == 1) {
					AdjacencyList[lineIndex].push_back(i);
				}
			}
			lineIndex++;
		}
	}	
}

Graph::~Graph() {
	
}

//##############################################################################
//#include "SpanningTree.h"
//##############################################################################

class SpanningTree {
	int startNode;
	int maxNodesCount;
	int degree;
	int nodesWithDegree;
	map<int, int> nodesDegree;

public:
	set<int> Nodes;
	vector<Edge*> Edges;
	
	SpanningTree(int nodesCount, int startNode);
	void AddEdge(Edge* edge);
	void RemoveEdge(Edge* edge);
	bool IsSpanningTree();
	int GetDegree();
	void Reset();
};

//######################

SpanningTree::SpanningTree(int nodesCount, int startNode) {
	this->maxNodesCount = nodesCount;
	this->startNode = startNode;
	Nodes.insert(startNode);
	nodesDegree[startNode] = 0;
	this->degree = 0;
    this->nodesWithDegree = 1;
}

void SpanningTree::Reset() {
	Nodes.clear();
	Nodes.insert(startNode);
	Edges.clear();
	nodesDegree.clear();
	nodesDegree[this->startNode] = 0;
	this->degree = 0;
	this->nodesWithDegree = 1;
}

void SpanningTree::AddEdge(Edge* edge) {

	Nodes.insert(edge->ToNode);
	nodesDegree[edge->ToNode] = 1;

	nodesDegree[edge->FromNode]++;
	int fromNodeDegree = nodesDegree[edge->FromNode];
	if(fromNodeDegree == degree)
	{
		nodesWithDegree++;
	}
	if(fromNodeDegree > degree)
	{
		degree = fromNodeDegree;
		nodesWithDegree = 1;
	}

	Edges.push_back(edge);
}

void SpanningTree::RemoveEdge(Edge* edge) {
	
	Nodes.erase(edge->ToNode);
	nodesDegree.erase(edge->ToNode);

	bool isFromNodeMaxDegree = nodesDegree[edge->FromNode] == degree;
	nodesDegree[edge->FromNode]--;

	if(isFromNodeMaxDegree)
	{
		if(nodesWithDegree > 1) 
		{
			nodesWithDegree--;
		}
		else
		{
			int maxDegree = 0;
            int nodesWithMaxDegree = 0;
			set<int>::iterator it;
			for (it=Nodes.begin(); it!=Nodes.end(); ++it)
			{
				int node = *it;
				int nodeDegree = nodesDegree[node];
				if(nodeDegree == maxDegree)
				{
					nodesWithMaxDegree++;
				}
				if(nodeDegree > maxDegree)
				{
					maxDegree = nodeDegree;
					nodesWithMaxDegree = 1;
				}
			}
			degree = maxDegree;
			nodesWithDegree = nodesWithMaxDegree;
		}
	}
	
	Edges.pop_back();
}

bool SpanningTree::IsSpanningTree() {
	return (Nodes.size() == maxNodesCount && Edges.size() 
                         == (maxNodesCount - 1));
}

int SpanningTree::GetDegree() {
	return degree;
}


//##############################################################################
//##############################################################################

#define CHECK_MSG_AMOUNT  100
#define MSG_WORK_REQUEST 1000
#define MSG_WORK_SENT    1001
#define MSG_WORK_NOWORK  1002
#define MSG_TOKEN        1003
#define MSG_FINISH       1004
#define MSG_NEW_MINIMUM  1005

#define STATUS_WORKING  2001
#define STATUS_IDLE     2002
#define STATUS_FINISHED 2999

#define WHITE_PROCESS 3001
#define BLACK_PROCESS 3002

#define WHITE_TOKEN 3003
#define BLACK_TOKEN 3004


Graph* LoadData(char* fileName);
void PrintGraph(Graph* g);
void PrintSpanningTree(SpanningTree* s);
bool sendWork(int target, MPI_Request* request, int* sendData, 
              int maxSendDataLenght, list<Edge*> &stack, SpanningTree* tree);
void sendFirstEdge(int target, int* sendData, Edge* edge);
void receiveWork(MPI_Status* status, int* sendData, int maxSendDataLenght, 
                 list<Edge*> &stack, SpanningTree* tree);
void receiveFirstEdge(MPI_Status* status, int* sendData, list<Edge*> &stack);

void sendNewMin(int minDegree, int my_rank, int processorCount);

//##############################################################################
// main
//##############################################################################

int main(int argc, char **argv) {
    if (argc != 2) {
        cout << "Incorrect program parameters" << endl;;
        return 1;
    }

    bool wasRequestForWork = false;
    int my_rank;
    int processorCount;
    unsigned int checkCounter;
    unsigned int processToAskForWork;
    int processStatus = 0;
    int token = 0;
    int color = WHITE_PROCESS;
    bool vectorSent = false;
    bool isResult = false;
    bool wasWorking = false;
	double time1, time2 = 0;

    MPI_Status status;
    MPI_Request request;

    int flag;
    const int minBound = 2;
    int i;
	
	Graph* graph = LoadData(argv[1]);
	if(graph == NULL)
	{
		cout << "Data load failed" << endl;
		return 2;
	}

    int maxSendDataLength = (graph->NodesCount * 2 + 1) * 2 + 2;
    int* sendData = new int[maxSendDataLength];

    //PrintGraph(graph);

    int startNode = 0;
    int minDegree = numeric_limits<int>::max();
    vector<Edge> minSpanningTree;

    SpanningTree* spanningTree = new SpanningTree(graph->NodesCount, startNode);
    list<Edge*> myStack;
    Edge* currentEdge;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processorCount);
	MPI_Barrier(MPI_COMM_WORLD);
	cout << "MPI_Barrier START wit MPI_Rank nr.: " << my_rank << endl;

	if(my_rank == 0) {
		time1 = MPI_Wtime();
	}
		
	if (my_rank == 0) {
		int startAdjacencySize = graph->AdjacencyList[startNode].size();
		if(startAdjacencySize < 1){
			cout << "Discontinuous Graf, startNode 0 neighbor" << endl;
			return 3;
		}

		cout << "starting to distribute" << endl;

		int lastSentAdjacent = 1;

		for(i = 1; i < processorCount; i++) 
		{
			if(i < startAdjacencySize)
			{
				int adjacent = graph->AdjacencyList[startNode][i];
				Edge* e = new Edge(startNode, adjacent, false);
				sendFirstEdge(i, sendData, e);
				lastSentAdjacent++;
				cout << "sent work to process " << i << endl;
			}
			else
			{
				MPI_Send(0, 0, MPI_INT, i, MSG_WORK_NOWORK, MPI_COMM_WORLD);
				cout << "sent NO_WORK to process " << i << endl;
			}
		}

		int adjacentForMe = graph->AdjacencyList[startNode][0];
		myStack.push_back(new Edge(startNode, adjacentForMe, false));
		cout << "Pridelil jsem si prvniho souseda" << endl;
		while(lastSentAdjacent < startAdjacencySize){
			adjacentForMe = graph->AdjacencyList[startNode][i];
			myStack.push_back(new Edge(startNode, adjacentForMe, false));
			cout << "Pridelil jsem si " << lastSentAdjacent 
                 << " souseda" << endl;
			lastSentAdjacent++;
		}

		processStatus = STATUS_WORKING;

	} else {
		while (processStatus == 0) {
			MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, 
                       &flag, &status);
			if (flag) {

				switch (status.MPI_TAG) {
					case MSG_WORK_NOWORK:
						MPI_Recv(0, 0, MPI_INT, 0, MSG_WORK_NOWORK, 
                                 MPI_COMM_WORLD, &status);
						processStatus = STATUS_IDLE;
						break;
					case MSG_WORK_SENT:
						receiveFirstEdge(&status, sendData, myStack);
						processStatus = STATUS_WORKING;
						break;
				}
			}
		}

	}
	
	processToAskForWork = (my_rank + 1) % processorCount;

    MPI_Barrier(MPI_COMM_WORLD);
	cout << "MPI_Barrier INIT WORK SENT" << endl;
    
    while (processStatus != STATUS_FINISHED) {

        if ((my_rank == 0) && (processorCount == 1) && 
            (processStatus == STATUS_IDLE)) {
            processStatus = STATUS_FINISHED;
        }
        if (processStatus == STATUS_IDLE) {
            
            if ((my_rank == 0) && (!vectorSent)) {
                color = WHITE_PROCESS;
                token = WHITE_TOKEN;
                MPI_Send(&token, 1, MPI_INT, 1, MSG_TOKEN, MPI_COMM_WORLD);
                vectorSent = true;
	
            } else if ((token != 0) && (my_rank != 0)) {
    			int nextProcesRank = (my_rank + 1) % processorCount;
	            MPI_Send(&token, 1, MPI_INT, nextProcesRank, MSG_TOKEN, 
                         MPI_COMM_WORLD);
				color = WHITE_PROCESS;
                token = 0;
            }
            
            int target = processToAskForWork % processorCount;
            if (!wasRequestForWork) {
                if (target != my_rank) {
                    MPI_Send(0, 0, MPI_INT, target, MSG_WORK_REQUEST, 
                             MPI_COMM_WORLD);
                    wasRequestForWork = true;
                }
            }

            processToAskForWork++;
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, 
                       &status);

            if (flag) {
                switch (status.MPI_TAG) {
                    case MSG_WORK_REQUEST:
                        MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, 
                                 MPI_COMM_WORLD, &status);
                        MPI_Send(0, 0, MPI_INT, status.MPI_SOURCE, 
                                 MSG_WORK_NOWORK, MPI_COMM_WORLD);
                        break;
                    case MSG_WORK_SENT:
                        receiveWork(&status, sendData, maxSendDataLength, 
                                    myStack, spanningTree);
                        wasRequestForWork = false;
                        processStatus = STATUS_WORKING;
                        break;
                    case MSG_WORK_NOWORK:
                        MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, 
                                 MPI_COMM_WORLD, &status);
                        wasRequestForWork = false;
                        processStatus = STATUS_IDLE;
                        break;
                    case MSG_TOKEN:
                        MPI_Recv(&token, 1, MPI_INT, status.MPI_SOURCE, 
                                 MSG_TOKEN, MPI_COMM_WORLD, &status);
                        if ((my_rank == 0) && (token == BLACK_TOKEN)) {
                            cout << "MASTER PRIJAL BLACK TOKEN" << endl;
                            vectorSent = false;
                        }
                        if ((my_rank == 0) && (token == WHITE_TOKEN)) {
                            processStatus = STATUS_FINISHED;
                            cout << "MAIN PROCESS STATUS_FINISHED" << endl;
                        } else {
                            if ((color == BLACK_PROCESS)) {
                                token = BLACK_TOKEN;
                            }
                        }
                        break;
                    case MSG_NEW_MINIMUM:
                        int message;
                        MPI_Recv(&message, 1, MPI_INT, MPI_ANY_SOURCE, 
                                 MSG_NEW_MINIMUM, MPI_COMM_WORLD, &status);
                        if(message == minBound){
                            isResult = false;
			    cout << "Prijat minBound procesem " << my_rank << endl;
                            processStatus = STATUS_FINISHED;
                        }
                        else if (message < minDegree) {
                            minDegree = message;
                            isResult = false;
                        }
                        break;
                    case MSG_FINISH:
                        MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, MSG_FINISH, 
                                 MPI_COMM_WORLD, &status);
                        processStatus = STATUS_FINISHED;
                        break;
                    default:;
                        break;
                }
            }
        }

        if (processStatus == STATUS_WORKING) {
            while (!myStack.empty() && processStatus != STATUS_FINISHED) {
                wasWorking = true;
				if(processorCount > 1){
					checkCounter++;
					if ((checkCounter % CHECK_MSG_AMOUNT) == 0) {
						MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, 
                                   &flag, &status);
						if (flag) {
							bool workSended;
							switch (status.MPI_TAG) {
								case MSG_WORK_REQUEST:
									MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, 
                                             MSG_WORK_REQUEST, MPI_COMM_WORLD, 
                                             &status);
									workSended = sendWork(status.MPI_SOURCE, 
                                                          &request, sendData, 
                                                          maxSendDataLength, 
                                                          myStack, 
                                                          spanningTree);
										if ((workSended) && 
                                            (my_rank > status.MPI_SOURCE)) {
											color = BLACK_PROCESS;
											if (token != 0) {
												token = BLACK_TOKEN;
											}
										}
									break;
								case MSG_NEW_MINIMUM:
									int message;
									MPI_Recv(&message, 1, MPI_INT, 
                                             MPI_ANY_SOURCE, MSG_NEW_MINIMUM, 
                                             MPI_COMM_WORLD, &status);
									if(message == minBound){
										isResult = false;
										cout << "Prijat minBound procesem " 
                                             << my_rank << endl;
										processStatus = STATUS_FINISHED;
									}
									else if (message < minDegree) {
										minDegree = message;
										isResult = false;
									}
									break;
								case MSG_TOKEN:
									MPI_Recv(&token, 1, MPI_INT, MPI_ANY_SOURCE, 
                                             MSG_TOKEN, MPI_COMM_WORLD, 
                                             &status);
									if (my_rank == 0) {
										vectorSent = false;
									}
									break;
								default:;
									break;
							}
						}
					}
				}

                currentEdge = myStack.back();
                if (!currentEdge->IsVisited) {
                    currentEdge->IsVisited = true;
                    spanningTree->AddEdge(currentEdge);

                    int spanningTreeDegree = spanningTree->GetDegree();
                    if (spanningTreeDegree >= minDegree) {
                        continue;
                    }
                    if (spanningTree->IsSpanningTree()) {
                        if (spanningTreeDegree < minDegree) {
                            minDegree = spanningTreeDegree;
                            sendNewMin(minDegree, my_rank, processorCount);
							isResult = true;

                            minSpanningTree.clear();
                            int i;
                            for (i = 0; i < spanningTree->Edges.size(); i++) {
                                minSpanningTree.push_back(
                                  Edge(spanningTree->Edges[i]->FromNode, 
                                       spanningTree->Edges[i]->ToNode, false));
                            }
                        }
                        if (spanningTreeDegree == minBound) {
                            processStatus = STATUS_FINISHED;
                            break;
                        }
                    }

                    set<int>::iterator it;
                    int i;
                    for (it = spanningTree->Nodes.begin(); 
                         it != spanningTree->Nodes.end(); ++it) {

                        int treeNode = *it;
                        for (i = 0; 
                             i < graph->AdjacencyList[treeNode].size(); i++) {
                            int adjacent = graph->AdjacencyList[treeNode][i];
                            if (spanningTree->Nodes.find(adjacent) == 
                                spanningTree->Nodes.end()) {
                                myStack.push_back(
                                  new Edge(treeNode, adjacent, false));
                            }
                        }
                    }

                } else {
                    myStack.pop_back();
                    spanningTree->RemoveEdge(currentEdge);
                    delete currentEdge;
                }
            }
            
            if(processStatus != STATUS_FINISHED) {
                processStatus = STATUS_IDLE;
            }
        }
    }

    delete graph;
    delete spanningTree;
    delete[] sendData;
	
    if (my_rank == 0 && minDegree != minBound) {
        for (i = 1; i < processorCount; i++) {
            MPI_Send(0, 0, MPI_INT, i, MSG_FINISH, MPI_COMM_WORLD);
            cout << "[MPI_Send] " << "Sender: " << my_rank << " Target " << i 
                 << " Message: " << MSG_FINISH << endl;
        }
    }

	MPI_Barrier(MPI_COMM_WORLD);
	cout << "MPI_Barrier END" << endl;

	if(my_rank == 0) {
		time2 = MPI_Wtime();

		double totalTime = time2 - time1;

		cout << "Total time is: " << totalTime << endl;
		printf ("Total time is %f.\n", totalTime);
	}
    if (isResult) {
        cout << "Minimum degree: " << minDegree << endl;
        cout << "Spanning Tree:" << endl;
        for (i = 0; i < minSpanningTree.size(); i++) {
            cout << minSpanningTree[i].FromNode << "->" 
                 << minSpanningTree[i].ToNode << endl;
        }
    }
    if (wasWorking){
        cout << "PROCESOR: " <<my_rank<<" pracoval." <<endl;
    }

    MPI_Finalize();
    return (0);

}

//##############################################################################
// sendNewMin
//##############################################################################
void sendNewMin(int minDegree, int resultRank, int processorCount) {
    for (int i = 0; i < processorCount; i++) {
        if (i != resultRank) {
            MPI_Send(&minDegree, 1, MPI_INT, i, MSG_NEW_MINIMUM, 
                     MPI_COMM_WORLD);
			
        }
    }
}

//##############################################################################
// sendWork
//##############################################################################
bool sendWork(int target, MPI_Request* request, int* sendData, 
              int maxSendDataLenght, list<Edge*> &stack, SpanningTree* tree) {
    list<Edge*>::reverse_iterator it;

    int edgesToSendCount = 0;

    it = stack.rbegin();
    while (it != stack.rend()) {
        if ((*it)->IsVisited) {
            break;
        }
        edgesToSendCount++;
        it++;
    }

    if (edgesToSendCount < 2) {
        MPI_Send(0, 0, MPI_INT, target, MSG_WORK_NOWORK, MPI_COMM_WORLD);
        return false;
    }

    edgesToSendCount = edgesToSendCount / 2;

    int dataLenght = (edgesToSendCount * 2) + 2 + (tree->Edges.size() * 2);

    it = stack.rbegin();

    int iSend = 0;
    for (int i = 0; i < edgesToSendCount; i++) {
        Edge* e = *it;

        sendData[iSend] = e->FromNode;
        sendData[iSend + 1] = e->ToNode;
        stack.remove(e);
        delete e;

        iSend += 2;
        it++;
    }
    sendData[iSend] = -1;
    iSend++;
    int treeSize = tree->Edges.size();
    for (int i = 0; i < treeSize; i++) {
        sendData[iSend] = tree->Edges[i]->FromNode;
        sendData[iSend + 1] = tree->Edges[i]->ToNode;
        iSend += 2;
    }
    sendData[iSend] = -1;

    MPI_Send(sendData, dataLenght, MPI_INT, target, MSG_WORK_SENT, 
             MPI_COMM_WORLD);
    return true;
}

//##############################################################################
// receiveWork
//##############################################################################
void receiveWork(MPI_Status* status, int* sendData, int maxSendDataLenght, 
                 list<Edge*> &stack, SpanningTree* tree) {
    tree->Reset();

    MPI_Recv(sendData, maxSendDataLenght, MPI_INT, MPI_ANY_SOURCE, 
             MSG_WORK_SENT, MPI_COMM_WORLD, status);
    int dataLenght = maxSendDataLenght;
    MPI_Get_count(status, MPI_INT, &dataLenght);

    bool wasMinusOne = false;
    for (int i = 2; i < dataLenght; i += 2) {
        if (wasMinusOne) {
            tree->AddEdge(new Edge(sendData[i - 2], sendData[i - 1], false));
        } else {
            stack.push_back(new Edge(sendData[i - 2], sendData[i - 1], false));
        }

        if (sendData[i] == -1) {
            if (i == dataLenght - 1) {
                break;
            }
            i++;
            if (sendData[i] == -1) {
                break;
            }
            wasMinusOne = true;
        }
    }
}

//##############################################################################
// sendFirstEdge
//##############################################################################
void sendFirstEdge(int target, int* sendData, Edge* edge) {
    sendData[0] = edge->FromNode;
    sendData[1] = edge->ToNode;
    MPI_Send(sendData, 2, MPI_INT, target, MSG_WORK_SENT, MPI_COMM_WORLD);
    delete edge;
}

//##############################################################################
// receiveFirstEdge
//##############################################################################
void receiveFirstEdge(MPI_Status* status, int* sendData, list<Edge*> &stack) {
    MPI_Recv(sendData, 2, MPI_INT, MPI_ANY_SOURCE, MSG_WORK_SENT, 
             MPI_COMM_WORLD, status);
    stack.push_back(new Edge(sendData[0], sendData[1], false));
}

//##############################################################################
// PrintGraph
//##############################################################################
void PrintGraph(Graph* g) {
    cout << "Nodes Count: " << g->NodesCount << endl;
    int i, j;
    for (i = 0; i < g->AdjacencyList.size(); i++) {
        for (j = 0; j < g->AdjacencyList[i].size(); j++) {
            cout << g->AdjacencyList[i][j];
        }
        cout << endl;
    }
}

//##############################################################################
// PrintSpanningTree
//##############################################################################
void PrintSpanningTree(SpanningTree* s) {
    cout << "Edges in spanning tree:" << endl;
    int i;
    for (i = 0; i < s->Edges.size(); i++) {
        cout << s->Edges[i]->FromNode << "->" << s->Edges[i]->ToNode << endl;
    }
}

//##############################################################################
// LoadData
//##############################################################################
Graph* LoadData(char* fileName) {
	ifstream infile;
	infile.open(fileName);

	if(!infile) { 
		cout << "Cannot open file.\n"; 
		return NULL; 
	} 

	Graph* g = new Graph(infile); 
	
	infile.close();

	return g;
}

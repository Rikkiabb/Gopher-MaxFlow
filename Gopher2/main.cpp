#include <algorithm>
#include <bitset>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <utility>
#include <vector>

using namespace std;


#include <cstdlib>
#include <vector>
#include <stack>
#include <queue>
using namespace std;

const double EPS = 1e-9;

struct flow_network {
    
    struct edge {
        int u, v, cap;
        edge *rev;
        bool forward;
        edge(int _u, int _v, int _cap, bool forw)
        : u(_u), v(_v), cap(_cap), forward(forw) { }
    };
    
    int n;
    vector<vector<edge*> > adj;
    flow_network(int _n) : n(_n), adj(_n) { }
    
    void add_edge(int u, int v, int cap) {
        edge *e = new edge(u, v, cap, true);
        edge *rev = new edge(v, u, 0, false);
        e->rev = rev;
        rev->rev = e;
        adj[u].push_back(e);
        adj[v].push_back(rev);
    }
    
    int augment(int s, int t) {
        vector<edge*> back(n, NULL);
        queue<int> Q;
        Q.push(s);
        back[s] = (edge*)1;
        while (!Q.empty()) {
            int u = Q.front(); Q.pop();
            for (int i = 0; i < adj[u].size(); i++) {
                int v = adj[u][i]->v;
                if (back[v] == NULL && adj[u][i]->cap > 0) {
                    back[v] = adj[u][i];
                    Q.push(v);
                }
            }
        }
        
        if (back[t] == NULL)
            return 0;
        
        stack<edge*> S;
        S.push(back[t]);
        int bneck = back[t]->cap;
        while (S.top()->u != s) {
            S.push(back[S.top()->u]);
            bneck = min(bneck, S.top()->cap);
        }
        
        while (!S.empty()) {
            S.top()->cap -= bneck;
            S.top()->rev->cap += bneck;
            S.pop();
        }
        
        return bneck;
    }
    
    int max_flow(int source, int sink) {
        int flow = 0;
        while (true) {
            int f = augment(source, sink);
            if (f == 0) {
                break;
            }
            
            flow += f;
        }
        
        return flow;
    }
};

double coordsLength(double x, double y, double x2, double y2){

    return sqrt(pow((x - x2), 2) + (pow((y - y2), 2)));

}

struct coords{
    
    double x;
    double y;
    coords(double _x, double _y){
    
        x = _x;
        y = _y;
    
    }

};


int main() {

    
    int n = 0, m = 0, s = 0, v = 0; // n = gopher, m = gopher holes, s = seconds, v = gophers run at velocity v
    while(cin >> n >> m >> s >> v){
    
    int SOURCE = 0,
    SINK = 1,
    LEFT = 2,
    RIGHT = LEFT + n,
    CNT = RIGHT + m;

    flow_network g(CNT);
    double maxLength = s * v;
    vector<coords>gophers;
    vector<coords>holes;

    for(int i = 0; i < n; i++){ // Coordinates of gophers
    
        double x = 0, y = 0;
        scanf("%lf %lf", &x, &y);
        coords newc(x, y);
        gophers.push_back(newc);
    
    }
    
    for(int i = 0; i < n; i++){ // Connect gophers to source
    
        g.add_edge(SOURCE, LEFT + i, 1);
    
    }
    
    for(int i = 0; i < m; i++){ // Coordinates of gopher holes
    
        double x = 0, y = 0;
        scanf("%lf %lf", &x, &y);
        coords newc(x, y);
        holes.push_back(newc);
        
    }
    
    for(int i = 0; i < m; i++){ // Connect holes to sink
    
        g.add_edge(RIGHT + i, SINK, 1);
    
    }
    
    for(int i = 0; i < n; i++){ // Connect gophers to holes
    
        for(int j = 0; j < m; j++){
        
            double length = coordsLength(gophers[i].x, gophers[i].y, holes[j].x, holes[j].y);
            if(length <= maxLength + EPS){
            
                    g.add_edge(LEFT + i, RIGHT + j, 1);
            
            }
        
        }
    
    }

    printf("%d\n", n - g.max_flow(SOURCE, SINK));
    }
    
    return 0;
}



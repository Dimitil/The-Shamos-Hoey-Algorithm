#include <iostream>
#include <list>
#include <QPointF>
#include <cmath>
#include <vector>
#include <set>
#include <QDebug>
#include <QPolygonF>
#include <QtGlobal>
#define FALSE   0
#define TRUE    1
#define LEFT    0
#define RIGHT   1
const double EPS = 0.000001;
// xyorder(): determines the xy lexicographical order of two points
//      returns: (+1) if p1 > p2; (-1) if p1 < p2; and  0 if equal
int xyorder( const QPointF* p1, const QPointF* p2 )
{
    // test the x-coord first
    if (p1->x() - p2->x() > EPS) return 1;
    if (p2->x() - p1->x() > EPS) return (-1);
    // and test the y-coord second
    if (p1->y() - p2->y() > EPS ) return 1;
    if (p2->y() - p1->y() > EPS ) return (-1);
    // when you exclude all other possibilities, what remains  is->->->
    return 0;  // they are the same point
}

// isLeft(): tests if point P2 is Left|On|Right of the line P0 to P1.
//      returns: >0 for left, 0 for on, and <0 for  right of the line.
inline double isLeft( QPointF P0, QPointF P1, QPointF P2 )
{
    return (P1.x() - P0.x())*(P2.y() - P0.y()) - (P2.x() - P0.x())*(P1.y() -  P0.y());
}

typedef struct _event Event;
struct _event {
    int      edge;          // polygon edge i is V[i] to V[i+1]
    int      type;          // event type: LEFT or RIGHT vertex
    QPointF*   eV;            // event vertex
};

int E_compare( const void* v1, const void* v2 ) // qsort compare two events
{
    Event**    pe1 = (Event**)v1;
    Event**    pe2 = (Event**)v2;

    return xyorder( (*pe1)->eV, (*pe2)->eV );
}
// the EventQueue is a presorted array (no insertions needed)
class EventQueue {
    int      ne;                // total number of events in array
    int      ix;                // index of next event on queue
    Event*   Edata;             // array of all events
    Event**  Eq;                // sorted list of event pointers
public:
              EventQueue(QPolygonF P);     // constructor
             ~EventQueue(void)           // destructor
                  { delete[] Eq; delete[] Edata;}

    Event*   next();                     // next event on queue
};
// EventQueue Routines
EventQueue::EventQueue(QPolygonF P )
{

    ix = 0;
    ne = 2 * (P.size()-1) ;           // vertex events for each point
    Edata = new Event[ne];
    Eq =  new Event*[ne];
    for (int i=0; i < ne; i++)           // init Eq array pointers
        Eq[i] = &Edata[i];

    // Initialize event queue with edge segment endpoints
    for (int i=0, j = 0; j < P.size() - 1; ) {        // init data for edge i
           Eq[i]->eV   = &(P[j]);
           i++;
           j++;
           Eq[i]->eV = &(P[j]);
           i++;
    }
    for ( int i = 0, j = 0; i < ne - 1; i +=2 ){
           Eq[i]->edge = j;
           Eq[i+1]->edge = j;
           j++;
           if (xyorder( Eq[i]->eV, Eq[i+1]->eV) < 0)  { // determine type
               Eq[i]->type    = LEFT;
                Eq[i+1]->type = RIGHT;
           }
           else {
               Eq[i]->type    = RIGHT;
                Eq[i+1]->type = LEFT;
           }
       }
    // Sort Eq[] by increasing x and y
    qsort(Eq, ne, sizeof(Event**), E_compare);
}
Event* EventQueue::next()
{
    if (ix >= ne)
        return nullptr;
    else
        return Eq[ix++];
}


// SweepLine Class

// SweepLine segment data struct
typedef struct _SL_segment SLseg;
struct _SL_segment {
    int      edge;          // polygon edge i is V[i to V[i+1]
    QPointF    lP;            // leftmost vertex pointi
    QPointF    rP;            // rightmost vertex point
    SLseg*   above;         // segment above this one
    SLseg*   below;         // segment below this one
};

struct cmp{
    bool operator()(SLseg* u, SLseg* d)
    {
        double uy1 = u->lP.y();
        double uy2 = u->rP.y();
        double dy1 = d->lP.y();
        double dy2 = d->rP.y();
        double minUY = std::min(uy1, uy2);
        double minDU = std::min(dy1, dy2);
        return minUY > minDU;
    }
};
// the Sweep Line itself
class SweepLine {
    int      nv;            // number of vertices in polygon
    QPolygonF* Pn;            // initial Polygon
    std::multiset<SLseg*, cmp>     Tree;          // balanced binary tree
public:
              SweepLine(QPolygonF P)            // constructor
                  { nv = P.size(); Pn = &P; }
             ~SweepLine(void)                 // destructor
                  { Tree.clear();}

    SLseg*   add( Event* );
    SLseg*   find( Event* );
    int      intersect( SLseg*, SLseg*  );
    void     remove( SLseg* );
};

SLseg* SweepLine::add( Event* E )
{
    // fill in SLseg element data
    SLseg* s = new SLseg;
    s->edge  = E->edge;

    // if it is being added, then it must be a LEFT edge event
    // but need to determine which endpoint is the left one
    QPointF* v1 = &((*Pn)[s->edge]);
    QPointF* v2 = &((*Pn)[s->edge+1]);
    if (xyorder( v1, v2) < 0) { // determine which is leftmost
        s->lP = *v1;
        s->rP = *v2;
    }
    else {
        s->rP = *v1;
        s->lP = *v2;
    }
    s->above = (SLseg*)0;
    s->below = (SLseg*)0;

    // add a node to the balanced binary tree
    auto nd = Tree.insert(s);
    auto nx = nd;
    nx++;
    auto np = nd;
    if (np == Tree.begin())
    {
        np = Tree.end();
    }
    else {
    --np;
    }
    if (nx != Tree.end()) {
        s->above = (SLseg*)(*nx);
        s->above->below = s;
    }
    if (np != Tree.end()) {
        s->below = (SLseg*)*np;
        s->below->above = s;
    }
    return s;
}
#include<algorithm>

SLseg* SweepLine::find( Event* E )
{
    // need a segment to find it in the tree
    SLseg* s = new SLseg;
    s->edge  = E->edge;
    s->above = (SLseg*)0;
    s->below = (SLseg*)0;

    auto nd = std::find_if(Tree.begin(), Tree.end(), [&ed=s->edge](SLseg* seg){return ed==seg->edge;});
    delete s;
    if (nd == Tree.end())
        return nullptr;

    return *nd;
}

void SweepLine::remove( SLseg* s )
{
    // remove the node from the balanced binary tree
    auto nd = Tree.find(s);
    if (nd == Tree.end())
        return;       // not there

    // get the above and below segments pointing to each other
    auto nx = nd;
    ++nx;
    if (nx != Tree.end()) {
        SLseg* sx = (SLseg*)(*nx);
        sx->below = s->below;
    }
    auto np = nd;
    --np;
    if (np != --Tree.begin()) {
        SLseg* sp = (SLseg*)(*np);
        sp->above = s->above;
    }
    Tree.erase(nd);       // now  can safely remove it
    delete s;
}

// test intersect of 2 segments and return: 0=none, 1=intersect
int SweepLine::intersect( SLseg* s1, SLseg* s2)
{
    if (s1 == (SLseg*)0 || s2 == (SLseg*)0)
        return FALSE;       // no intersect if either segment doesn't exist

    // check for consecutive edges in polygon
    int e1 = s1->edge;
    int e2 = s2->edge;
    if (((e1+1)%nv == e2) || (e1 == (e2+1)%nv))
        return FALSE;       // no non-simple intersect since consecutive
    if(e1 == Pn->size() - 2)
    {
        return FALSE;
    }

    // test for existence of an intersect point
    float lsign, rsign;
    lsign = isLeft(s1->lP, s1->rP, s2->lP);    //  s2 left point sign
    rsign = isLeft(s1->lP, s1->rP, s2->rP);    //  s2 right point sign
    if (lsign * rsign > 0) // s2 endpoints have same sign  relative to s1
        return FALSE;       // => on same side => no intersect is possible
    lsign = isLeft(s2->lP, s2->rP, s1->lP);    //  s1 left point sign
    rsign = isLeft(s2->lP, s2->rP, s1->rP);    //  s1 right point sign
    if (lsign * rsign > 0) // s1 endpoints have same sign  relative to s2
        return FALSE;       // => on same side => no intersect is possible
    // the segments s1 and s2 straddle each other
    return TRUE;            // => an intersect exists
}
//=============================================
// simple_Polygon(): test if a Polygon is simple or not
//     Input:  Pn = a polygon with n vertices V[]
//     Return: FALSE(0) = is NOT simple
//             TRUE(1)  = IS simple
int
simple_Polygon( QPolygonF Pn)
{
    EventQueue  Eq(Pn);
    SweepLine   SL(Pn);
    Event*      e;                  // the current event
    SLseg*      s;                  // the current SL segment

    // This loop processes all events in the sorted queue
    // Events are only left or right vertices since
    // No new events will be added (an intersect => Done)
    e = Eq.next();
    while ( e) {         // while there are events
        if (e->type == LEFT) {      // process a left vertex
            s = SL.add(e);          // add it to the sweep line
            if (SL.intersect(  s, s->above)){
                 return FALSE;      // Pn is NOT simple
            }
            if (SL.intersect(  s, s->below))
                 return FALSE;      // Pn is NOT simple
        }
        else {                      // processs a right vertex
            s = SL.find(e);
            if (SL.intersect(  s->above, s->below))
                 return FALSE;      // Pn is NOT simple
            SL.remove(s);           // remove it from the sweep line
        }
        e = Eq.next();
    }
    return TRUE;      // Pn IS simple
}
int main()
{
    QPolygonF polygon;
    std::list<QPointF*> list1;
    for(int i = 0; i < 10; i++)
    {
        QPointF* pp = new QPointF(5855552+i*10, 25542343+i*20);
        list1.insert(list1.end(), pp);
        polygon << *pp;
    }
    QPointF* pp = new QPointF(5855552+100, 25542343);
    list1.insert(list1.end(), pp);
    polygon << *pp;
    QPointF* pp2 = new QPointF(5855552, 25542343);
    list1.insert(list1.end(), pp2);
    polygon << *pp2;

    QPolygonF noSim;
    noSim << QPointF(5, 5) << QPointF(0, 7.5) << QPointF(5, 10) << QPointF(10,7.5) << QPointF(20, 7.5) << QPointF(0,5) ;

    qDebug() <<"Polygon is" <<  simple_Polygon(noSim);
    qDebug() << "STOP";
    return 0;
}

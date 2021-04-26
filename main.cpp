//Shamos - Hoey Algorithm
/*Aлгоритм поиска пересекающихся отрезков в полигоне с помощью заметающей прямой
 (SweepLine), считается, что соседние и первый и последний отрезки не могут
пересекаться, а имеют общую точку.
структура TREE - содержит упорядоченые по xy указатели на отрезки, при таком
порядке пересекаться могут только соседние в TREE отрезки, что и проверяется
в intersect(). Проверку достаточно производить в точках начала и конца абсцис
отрезков, эти точки реализованы структурой Event и храняться в упорядочено по xy в
EventQueue. В функции simple_Polygon(const QPolygonF&) в основном цикле по EventQueue
: если точка - точка начала отрезка, то она добавляется в TREE и проверяется
со своими соседями на пересечение. Если эта точка - точка конца отрезка, то
соответствующая ей точка начала отрезка удаляется из TREE.

Algorithm of search of the intersect segments in the polygons by means of the
sweeping line (SweepLine), is considered that the next both, and first and last
segments cannot to be intersected, and have the common point.
the structure of TREE - contains the pointers ordered on xy on segments, at it
order only the next segments in TREE can be intersected, as it is checked
in intersect (). Check is enough to be made in points of the beginning and end of
abscissae segments, these points are realized by structure of Event and are stored in is ordered on xy in
EventQueue. As simple_Polygon(const QPolygonF&) in a main loope on EventQueue
: if a point - is a the beginning of a segment, then it is added to
TREE and checked with the neighbors on crossing. If this point is end of a segment,
then the point of the beginning of a segment corresponding to it is removed from TREE.
*/

#include <iostream>
#include <algorithm>
#include <set>
#include <QPolygonF>

#define TREE std::multiset<SLseg*, SLseg_compare>

#define LEFT_sxflib    0
#define RIGHT_sxflib   1


const double EPS_sxflib = 1E-9;

// xyorder(): determines the xy lexicographical order of two points
//      returns: (+1) if p1 > p2; (-1) if p1 < p2; and  0 if equal
int xyorder( const QPointF* p1, const QPointF* p2 )
{
    if (p1->x() - p2->x() > EPS_sxflib) return 1;
    if (p2->x() - p1->x() > EPS_sxflib) return (-1);
    if (p1->y() - p2->y() > EPS_sxflib ) return 1;
    if (p2->y() - p1->y() > EPS_sxflib ) return (-1);
    return 0;
}

// isLeft(): tests if point P2 is Left|On|Right of the line P0 to P1.
//      returns: >0 for left, 0 for on, and <0 for  right of the line.
inline double isLeft( QPointF P0, QPointF P1, QPointF P2 )
{
    return (P1.x() - P0.x())*(P2.y() - P0.y()) -
            (P2.x() - P0.x())*(P1.y() -  P0.y());
}

typedef struct _event Event;
struct _event {
    int      edge;
    int      type;          // event type: LEFT_sxflib or RIGHT_sxflib vertex
    const QPointF*   eV;
};

int E_compare( const void* v1, const void* v2 ) // qsort compare two events
{
    Event**    pe1 = (Event**)v1;
    Event**    pe2 = (Event**)v2;
    return xyorder( (*pe1)->eV, (*pe2)->eV );
}

class EventQueue {
    int      ne;                // total number of events in array
    int      ix;                // index of next event on queue
    Event*  Edata;              // array of all events
    Event**  Eq;                // sorted list of event pointers
    QPolygonF* P;

public:
    EventQueue(const QPolygonF& P)
    {
        ix = 0;
        ne = (P.size() - 1) * 2;           // vertex events for each point
        Edata = new Event[ne];
        Eq =  new Event*[ne];
        for (int i=0; i < ne; i++)           // init Eq array pointers
        {
            Eq[i] = &Edata[i];
        }

        // Initialize event queue with edge segment endpoints
        for (int j = 0, i = 0; i < P.size() - 1; )
        {
            Eq[j]->eV = &(P[i]);
            i++;
            j++;
            Eq[j]->eV = &(P[i]);
            Eq[j-1]->edge = i - 1;
            Eq[j]->edge   = i - 1;
            if (xyorder( Eq[j-1]->eV, Eq[j]->eV) < 0)  {
                Eq[j-1]->type = LEFT_sxflib;
                Eq[j]->type   = RIGHT_sxflib;
            }
            else {
                Eq[j - 1]->type = RIGHT_sxflib;
                Eq[j]->type    = LEFT_sxflib;
            }
            j++;
        }
        qsort(Eq, ne, sizeof(Event**), E_compare);
    }

    ~EventQueue(void)           { delete[] Eq; delete[] Edata;}
    const Event*   next()
    {
        if (ix >= ne)
        {
            return NULL;
        }
        else {
            return Eq[ix++];
        }
    }
};

typedef struct _SL_segment SLseg;
struct _SL_segment {
    int      edge;          // number of edge
    QPointF  lP;            // leftmost vertex pointi
    QPointF  rP;            // rightmost vertex point
};


//comparator for TREE
struct SLseg_compare{
    bool operator()(SLseg* up, SLseg* down)
    {
        double minUY = std::min(up->lP.y(), up->rP.y());
        double minDU = std::min(down->lP.y(), down->rP.y());
        return minUY > minDU;
    }
};

//comparator for find_if algoritm
struct SLseg_equal_by_edge{
    int val;
    bool operator()(SLseg* seg)
    {
        return val==seg->edge;
    }
};

class SweepLine {
    int      nv;            // number of vertices in polygon
    const QPolygonF* Pn;
    TREE Tree;
public:
    SweepLine(const QPolygonF& P)   { nv = (P.size() - 1); Pn = &P; }
    ~SweepLine()
    {
        for (TREE::iterator it = Tree.begin(); it != Tree.end(); ++it)
        {
            delete *it;
        }
        Tree.clear();
    }

    TREE::iterator   add(const Event* E)
    {
        SLseg* s = new SLseg;
        s->edge  = E->edge;
        const QPointF* v1 = &((*Pn)[s->edge]);
        const QPointF* v2 = &((*Pn)[s->edge+1]);
        if (xyorder( v1, v2) < 0)
        {
            s->lP = *v1;
            s->rP = *v2;
        }
        else {
            s->rP = *v1;
            s->lP = *v2;
        }
        return Tree.insert(s);
    }

    TREE::iterator   find(const Event* E)
    {
        SLseg s;
        s.edge = E->edge;
        SLseg_equal_by_edge eqSl;
        eqSl.val = s.edge;
        TREE::iterator nd = std::find_if(Tree.begin(), Tree.end(), eqSl);
        return nd;
    }

    bool intersect( TREE::iterator it1, TREE::iterator it2)
    {
        if(it1 == Tree.end() || it2 == Tree.end())
        {
            return false;
        }
        if(it1 == --Tree.begin() || it2 == --Tree.begin())
        {
            return false;
        }
        SLseg* s1 = *it1;
        SLseg* s2 = *it2;
        if (s1 == NULL || s2 == NULL)
        {
            return false;
        }

        int e1 = s1->edge;
        int e2 = s2->edge;
        if (((e1+1)%nv == e2) || (e1 == (e2+1)%nv))
        {
            return false;       // no non-simple intersect since consecutive
        }

        float lsign, rsign;
        lsign = isLeft(s1->lP, s1->rP, s2->lP);    //  s2 left point sign
        rsign = isLeft(s1->lP, s1->rP, s2->rP);    //  s2 right point sign
        if (lsign * rsign > 0) // s2 endpoints have same sign  relative to s1
            return false;       // => on same side => no intersect is possible
        lsign = isLeft(s2->lP, s2->rP, s1->lP);    //  s1 left point sign
        rsign = isLeft(s2->lP, s2->rP, s1->rP);    //  s1 right point sign
        if (lsign * rsign > 0) // s1 endpoints have same sign  relative to s2
            return false;       // => on same side => no intersect is possible
        // the segments s1 and s2 straddle each other
        return true;
    }

    void     remove( TREE::iterator nd)
    {
        delete *nd;
        Tree.erase(nd);
    }

};




// test intersect of 2 segments and return: 0=none, 1=intersect
//=============================================
// simple_Polygon(): test if a Polygon is simple or not (simple polygon have not
//                   self-intersect and other polygon inside)
//     Input:  Pn = a polygon with n vertices V[]
//     Return: FALSE(0) = is NOT simple
//             TRUE(1)  = IS simple
bool simple_Polygon(const QPolygonF& Pn)
{
    EventQueue  Eq(Pn);
    SweepLine   SL(Pn);
    const Event*      e;                  // the current event
    TREE::iterator  s;                  // the current SL segment

    // This loop processes all events in the sorted queue
    // Events are only left or right vertices since
    // No new events will be added (an intersect => Done)
    e = Eq.next();
    while ( e) {
        if (e->type == LEFT_sxflib) {
            s = SL.add(e);
            TREE::iterator above = s;
            above--;
            if (SL.intersect(  s, above)){
                 return false;
            }
            TREE::iterator below = s;
            below++;
            if (SL.intersect(  s, below))
                 return false;
        }
        else {
            s = SL.find(e);
            TREE::iterator above = s;
            above--;
            TREE::iterator below = s;
            below++;
            if (SL.intersect(  above, below))
                 return false;
            SL.remove(s);
        }
        e = Eq.next();
    }
    return true;
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
    noSim << QPointF(5, 5) << QPointF(0, 7.5) << QPointF(5, 10)
          << QPointF(10,7.5) << QPointF(10, 5)  << QPointF(5,5);

    QPolygonF p2;
    p2 << QPointF(5256196.73, 22615493.58) << QPointF(5256208.02, 22615535.35)
       << QPointF(5256231.33, 22615529.05) << QPointF(5256164.15, 22615375.26)
       << QPointF(5256220.04, 22615535.28);// << QPointF(5256196.73, 22615493.58);

    std::cout << "Polygon is ";
    if (simple_Polygon(p2))
    {
        std::cout << "simple!" << std::endl;
    }
    else {
        std::cout << "NOT simple!" << std::endl;
    }

    return 0;
}


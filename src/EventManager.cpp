#include "EventManager.h"
#include <cstdio>
#define EMPTY -1

struct EventItem{
    EventRef previous_;
    EventRef next_;
    size_t   qIndex_;
};

EventManager::EventManager(void):
    currentIndex_(0), baseIndex_(0),
    nCBTEvents_(0)
{}

EventManager::~EventManager(void){
}

void EventManager::resize(size_t nPart){
    events_.resize(nPart);
    eventItems_.resize(nPart);
    //Check these numbers
    leafs_.resize(nPart + 1);
    nodes_.resize(2 * nPart);
}

void EventManager::init(void){
    //Instrument queue
    {
        double tmax = 0.0;
        double tmin[2] = {100000.0};
        for(size_t i = 0; i < events_.size(); ++i){
            if(events_[i].empty()) continue;

            double time = events_[i].top().time_;
            double temp = tmin[0];
            tmin[0] = std::min(tmin[0], time);
            if(tmin[0] != temp) tmin[1] = temp;
            tmax = std::max(tmax, time);
        }
        double dtmin = tmin[1] - tmin[0];
        double dtmax = tmax - tmin[0];

        scaleFactor_ = 0.1 / dtmin;
        llSize_      = scaleFactor_ * dtmax;
    }

    llQueue_.resize(llSize_ + 1, EMPTY);

    for(size_t i = 0; i < events_.size(); ++i){
        if(!events_[i].empty()) insertInEventQ(i);
    }
}

//NOTE: Need to properly clear
void EventManager::clear(void){
    for(size_t i = 0; i < events_.size(); ++i) clear(i);
}

void EventManager::clear(size_t pID){
    events_[pID].clear();
}

bool EventManager::empty(size_t pID)const{
    return events_[pID].empty();
}

void EventManager::push(size_t pID, const ParticleEvent& event){
    events_[pID].push(event);
}

void EventManager::update(size_t pID){
    deleteFromEventQ(pID);
    insertInEventQ(pID);
}

void EventManager::insertInEventQ(size_t eRef){
    const ParticleEvent& event = events_[eRef].top();
    size_t index = (size_t)(scaleFactor_ * event.time_ - baseIndex_); 

    if(index > llSize_ - 1){
        index -= llSize_;
        if(index + 1 >= currentIndex_) index = llSize_;
    }

    EventItem& eItem = eventItems_[eRef];

    eItem.qIndex_ = index;

    if(index == currentIndex_) cbtInsert(eRef); //Insert in binary tree
    else{                                       //Insert in linked-list
        EventRef oldFirst = llQueue_[index];

        eItem.previous_ = EMPTY;
        eItem.next_     = oldFirst;
        llQueue_[index] = eRef;

        if(oldFirst != EMPTY) eventItems_[oldFirst].previous_ = eRef;
    }
}

void EventManager::processOverflowList(void){
    size_t index    = llSize_;
    EventRef eRef   = llQueue_[index];
    llQueue_[index] = EMPTY;

    while(eRef != EMPTY){
        EventRef next = eventItems_[eRef].next_;
        insertInEventQ(eRef);
        eRef = next;
    }
}

void EventManager::deleteFromEventQ(size_t eRef){
    const EventItem& eItem = eventItems_[eRef];
    size_t index = eItem.qIndex_;
    if(index == currentIndex_) cbtDelete(eRef);
    else{
        EventRef previous = eItem.previous_;
        EventRef next     = eItem.next_;
        if(previous == EMPTY) llQueue_[index] = next;
        else eventItems_[previous].next_ = next;

        if(next != EMPTY) eventItems_[next].previous_ = previous;
    }
}

ParticleEvent EventManager::getNextEvent(void){
    while(nCBTEvents_ == 0){
        ++currentIndex_;
        if(currentIndex_ == llSize_){
            currentIndex_  = 0;
            baseIndex_    += llSize_;
            processOverflowList();
        }

        //populate binary tree
        for(EventRef eRef = llQueue_[currentIndex_]; eRef != EMPTY; eRef = eventItems_[eRef].next_){
            cbtInsert(eRef);
        }

        llQueue_[currentIndex_] = EMPTY;
    }

    EventRef eRef = nodes_[1]; //The root contains the earliest event

    return events_[eRef].pop();
}

void EventManager::cbtUpdate(EventRef eRef){
    size_t father;
    for(father = leafs_[eRef] / 2; (nodes_[father] == eRef) && (father > 0); father /= 2){
        EventRef leftNode  = nodes_[2 * father];
        EventRef rightNode = nodes_[2 * father + 1];

        nodes_[father] = (events_[leftNode].top().time_ < events_[rightNode].top().time_)? leftNode : rightNode;
    }

    for( ; father > 0; father /= 2){
        EventRef oldWinner = nodes_[father];

        EventRef leftNode  = nodes_[2 * father];
        EventRef rightNode = nodes_[2 * father + 1];

        nodes_[father] = (events_[leftNode].top().time_ < events_[rightNode].top().time_)? leftNode : rightNode;

        if(nodes_[father] == oldWinner) return;
    }
}

void EventManager::cbtInsert(EventRef eRef){
    if(nCBTEvents_ == 0){
        nodes_[1] = eRef;
        ++nCBTEvents_;
        return;
    }

    EventRef lastNode           = nodes_[nCBTEvents_];
    nodes_[2 * nCBTEvents_]     = lastNode;
    nodes_[2 * nCBTEvents_ + 1] = eRef;

    leafs_[lastNode] = 2 * nCBTEvents_;
    leafs_[eRef]     = 2 * nCBTEvents_ + 1;

    ++nCBTEvents_;
    cbtUpdate(lastNode);
}

void EventManager::cbtDelete(EventRef eRef){
    if(nCBTEvents_ < 2){
        nodes_[1] = 0;
        leafs_[0] = 1;
        --nCBTEvents_;
        return;
    }

    size_t lastLeaf = 2 * nCBTEvents_ - 1;
    if(nodes_[lastLeaf - 1] == eRef){
        leafs_[nodes_[lastLeaf]] = lastLeaf / 2;
        nodes_[lastLeaf / 2] = nodes_[lastLeaf];
        cbtUpdate(nodes_[lastLeaf]);
        --nCBTEvents_;
        return;
    }

    leafs_[nodes_[lastLeaf - 1]] = lastLeaf / 2;
    nodes_[lastLeaf / 2] = nodes_[lastLeaf - 1];
    cbtUpdate(nodes_[lastLeaf - 1]);

    if(nodes_[lastLeaf] != eRef){
        nodes_[leafs_[eRef]] = nodes_[lastLeaf];
        leafs_[nodes_[lastLeaf]] = leafs_[eRef];
        cbtUpdate(nodes_[lastLeaf]);
    }

    --nCBTEvents_;
}


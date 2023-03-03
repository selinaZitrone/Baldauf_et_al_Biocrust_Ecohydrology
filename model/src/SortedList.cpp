/***************************************************************************************//**
 * SortedList.cpp
 * written by Britta Tietjen
 *
 * This is just because I didn't know the official list in c++ when I
 * programmed this
 *******************************************************************************************/

#include <cstddef>
#include <iostream>

#include "SortedList.h"

using namespace std;

/*******************************************************************************************
 * is the list empty
 *******************************************************************************************/
bool CellList::empty(){

	return (first->next==NULL);
}

/*******************************************************************************************
 * is the actual element the last one of the list
 *******************************************************************************************/
bool CellList::endpos(){

	return (position->next==NULL);
}

/*******************************************************************************************
 * start at the first list element
 *******************************************************************************************/
void CellList::reset(){

	position = first;
}

/*******************************************************************************************
 * one step forward in the list
 *******************************************************************************************/
void CellList::advance(){

	if (!endpos())
		position = position->next;
}

/*******************************************************************************************
 * destructor
 *******************************************************************************************/
void CellList::delAll(){

	reset();
	while(!empty())
		del();
}

/*******************************************************************************************
 * insert new element at the end of the list
 *******************************************************************************************/
void CellList::insertEnd(Element* newEl){
	while(!endpos())
		advance();
	Element* help = new Element();
	help->next = position->next;
	position->next = newEl;
	position->next->previous = position;
	newEl->next = help->next;
	delete help;
}

/******************************************************************************************
 * insert new element at the right place
 *******************************************************************************************/
void CellList::insertSorted(Element* newEl){

	reset();
	//advance until as well the number of contributing flows as the number of contributing flows of
	//the cell which provides runon is bigger than the numbers of the previous list element and
	//smaller than the following element
	while((!endpos()) && (position->next->noOfCells < newEl->noOfCells)){
		advance();
	}
	while((!endpos()) && (!(position->next->noOfCells > newEl->noOfCells)) && (position->next->noOfCellsPrev < newEl->noOfCellsPrev)){
		advance();
	}
	Element* help = new Element();
	help->next = position->next;
	position->next = newEl;
	position->next->previous = position;
	newEl->next = help->next;
	delete help;
}

/*******************************************************************************************
 * delete the actual element
 *******************************************************************************************/
void CellList::del(){

	if (!endpos()){
		Element* help = new Element();
		help->next = position->next;
		position->next = position->next->next;
		delete help;
	}
}

/*******************************************************************************************
 * constructor with an empty list
 *******************************************************************************************/
CellList::CellList(){

	first = new Element();
	first->xval = -1;
	first->yval = -1;
	first->next = NULL;
	first->previous = NULL;
	position = first;
}

/*******************************************************************************************
 * destructor
 *******************************************************************************************/
CellList::~CellList(){

	reset();
	while(!empty())
		del();
}


/*******************************************************************************************
 * constructor of the element
 *******************************************************************************************/
Element::Element(int x, int y, int cells, int cellsPrev){

	xval = x;
	yval = y;
	noOfCells = cells;
	noOfCellsPrev = cellsPrev;
}

/*******************************************************************************************
 * default constructor of an element
 *******************************************************************************************/
Element::Element(){

}

/*******************************************************************************************
 * destructor
 *******************************************************************************************/
Element::~Element(){

}


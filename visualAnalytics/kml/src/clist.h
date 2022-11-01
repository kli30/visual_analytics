/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#ifndef __CLIST__H
#define __CLIST__H
#include <iterator>
#include <new>
#include <fstream>
#include <cassert>


using namespace std;

namespace KML
{

	template<typename T>
	struct clist_node
	{
		typedef clist_node<T>* node_pointer;
		node_pointer prev;
		node_pointer next;
		T data;
	};

	template<typename T >
	struct clist_iterator : public iterator<bidirectional_iterator_tag,T>
	{
		typedef clist_node<T>* link_type;
		typedef size_t size_type;
		typedef clist_iterator<T>  self;


		link_type node;  // link to node;

		//constructor
		clist_iterator (link_type x) : node(x){}
		clist_iterator () : node(NULL){};
		clist_iterator (const clist_iterator& x) : node(x.node) {}

		bool operator == (const self& x) const { return node == x.node ;}
		bool operator != (const self& x ) const { return node != x.node ;}


		T& operator *() const {return (*node).data; }
		T* operator ->() const { return &(operator*());}
		self& operator ++()
		{
			node = (link_type) ((*node).next);
			return *this;
		}
		self operator ++ (int)
		{
			self tmp= *this;
			++*this;
			return tmp;
		}
		self& operator --()
		{
			node = (link_type)((*node).prev);
			return *this;
		}
		self operator --(int)
		{
			self tmp = *this;
			--*this;
			return tmp ;
		}
		self& operator += (int n)
		{
			for (int i=0; i<n; ++i)
			{
				operator++();
			}
			return *this;
		}
		self& operator -= (int n)
		{
			for (int i=0; i<n; ++i)
			{
				operator--();
			}
			return *this;
		}
	};

	template<typename T >
	class clist
	{
	protected:
		typedef clist_node<T> node_type;
	public:
		typedef node_type* link_type;

	protected:

		link_type rootNodeLink;
		link_type currentNodeLink;

		link_type create_node(const T& x)
		{
			link_type p = new node_type;
			p->data= (T)(x);
			return p;
		}
		void destroy_node(link_type p)
		{
			delete p;
		}

		size_t clist_size;

	public:
		typedef clist_iterator<T> iterator;
		typedef const clist_iterator<T> const_iterator;
		typedef size_t size_type;
		//typedef iterator::size_type size_type;
		//typedef iterator::distance_type distance_type;
		//typedef iterator::difference_type difference_type ;


	public:

		clist()
		{
			rootNodeLink= NULL;
			currentNodeLink=NULL;
			clist_size=0;
		}

		clist(  const clist<T>& src)
		{
			clist_size=0;
			rootNodeLink=NULL;
			currentNodeLink=NULL;

 			size_t n= src.size();
			iterator itSrc=src.begin();
			for (size_t index=0; index< n; ++index)
			{
				this->push_back(*itSrc);
				++itSrc;
			} // end of for loop
		}

		~ clist()
		{
			if (clist_size>0)
			{
				iterator tmpIt=begin();
				++tmpIt;
				while (tmpIt!= begin())
				{
					erase(tmpIt++);
				}

				destroy_node(begin().node);
				clist_size=0;
			}

		}

		iterator begin() {return (link_type)(rootNodeLink); }
		iterator begin() const {return (link_type)(rootNodeLink); }
		iterator end() {return rootNodeLink->prev; }
		iterator end() const {return rootNodeLink->prev; }
		bool empty() const { return clist_size==0; }
		size_type size()  const
		{
			return clist_size;
		}

		inline void push_back(const T& x)
		{
			++clist_size;

			if (1==clist_size) // first one;
			{
				rootNodeLink= create_node(x);
				rootNodeLink->next=rootNodeLink;
				rootNodeLink->prev=rootNodeLink;
				currentNodeLink=rootNodeLink;
			}
			else
			{

				link_type tmp_node= create_node(x);
				tmp_node->next=rootNodeLink;

				tmp_node->prev=rootNodeLink->prev;
				rootNodeLink->prev->next=tmp_node;
				rootNodeLink->prev=tmp_node;
			}
		}

		 void copy_data( clist<T>& src)
		{
			// assume size==
			assert(clist_size==src.size());
			iterator itSrc=src.begin();
			iterator itDes=this->begin();
			for (size_t index=0;index<clist_size; ++index)
			{
				*itDes=*itSrc;
				++itDes;
				++itSrc;
			} // end of for loop
		}

		void clear()
		{
			if ( clist_size > 1)
			{
				iterator it= this->begin();
				++it;
				while(it!= this->end())
				{
					iterator current=it;
					++it;
					erase(current);
				} // end of for loop

				// destroy the begin and end;
				destroy_node(this->end().node);
				destroy_node(this->begin().node);
			}
			else if (clist_size==1)
			{
				destroy_node(this->begin().node);
			}
			else
			{
				//doing nothing;
			}

			clist_size=0;
		}

		iterator erase(iterator position)
		{
			if (clist_size>1)
			{
				--clist_size;
				link_type next= link_type( position.node->next );
				link_type prex= link_type( position.node->prev );
				prex->next= next;
				next->prev=prex;
				destroy_node(position.node);
				return iterator(next);
			}
			else if (clist_size==1)  // only the head node;
			{
				--clist_size;
				destroy_node(position.node);
				return iterator();
			}
			else
			{
				return iterator();
			}

		}

		//void pop_front() { erase(begin());}

		T& current_value(void)
		{
			return currentNodeLink->data;
		}
		T& prev_value(void)
		{
			return currentNodeLink->prev->data;
		}
		T& next_value(void)
		{
			return currentNodeLink->next->data;
		}
		void move_on(void)
		{
			currentNodeLink=currentNodeLink->next;
		}
		void set_current_pos(iterator newPos)
		{
			this->currentNodeLink= newPos.node;
		}
		void reset_current_pos(void)
		{
			this->set_current_pos(rootNodeLink);
		}

		iterator get_current_pos(void)
		{
			return iterator(this->currentNodeLink);
		}
		void set_current_value(T newValue)
		{
			currentNodeLink->data= newValue;
		}
		void erase_current()
		{
			erase(iterator(currentNodeLink));
		}
		void move_back(void)
		{
			currentNodeLink=currentNodeLink->prev;
		}

	};


	template< typename T>
	fstream& operator<< (fstream& os,  clist<T>& src)
	{
		src.set_current_pos(src.begin());
		for(int index=0; index <= src.size(); ++index)
		{
			os<< src.current_value()<<"\t";
			src.move_on();
		}
		os<<endl;
		return os;
	}

/*	template< typename T>
	ofstream& operator<< (ofstream& os, const clist<T>& src)
	{
		for (clist<T>::iterator itSrc= src.begin(); itSrc != src.end(); ++itSrc)
		{
			os<<*itSrc<<"\t";
		}

		os<<*(src.end())<<endl;
		return os;
	}

	template<typename T>
	ostream& operator<< (ostream& os, clist<T>& src)
	{
		for (clist<T>::iterator itSrc= src.begin(); itSrc != src.end(); ++itSrc)
		{
			os<<*itSrc<<"\t";
		}

		os<<*(src.end())<<endl;
		return os;

	}
*/

} // end of KML;



#endif

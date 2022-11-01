/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#pragma  once
namespace KML
{

	// a class contain edge information like two points of the edge, the nbr cells of this edge;
	class CEdgeInfoType
	{
	private:
		int pointId1;
		int pointId2;
		int cellId1;
		int cellId2;
	public:

		friend  class CTriSurface;
		CEdgeInfoType(): pointId2(0), pointId1(0){};
		CEdgeInfoType(int x,int y):pointId2(y),pointId1(x){};
		~CEdgeInfoType(){};

		bool operator==(const CEdgeInfoType& oth) const
		{
			return (this->pointId1==oth.pointId1 && this->pointId2==oth.pointId2 );
		}
		bool operator !=(const CEdgeInfoType& oth) const
		{
			return ! (this->operator ==(oth));
		}
		bool operator > (const CEdgeInfoType& oth) const
		{
			return (this->pointId2 > oth.pointId2 || ( this->pointId2 == oth.pointId2 && this->pointId1 > oth.pointId1) );
		}
		bool operator <(const CEdgeInfoType& oth) const
		{
			return oth.operator >(*this);
		}
		bool operator >=(const CEdgeInfoType& oth) const
		{
			return (this->operator >(oth) || this->operator ==(oth));
		}
		bool operator <=(const CEdgeInfoType& oth) const
		{
			return (this->operator <(oth) || this->operator ==(oth));
		}
		void Sort(void)
		{
			if( this->pointId1> this->pointId2)
			{
				int tmp=this->pointId1;
				this->pointId1=this->pointId2;
				this->pointId2=tmp;
			}
		}
		int& GetPoint1(void) {return this->pointId1;};
		int& GetPoint2(void) {return this->pointId2;};
		int& GetCellId1(void) { return this->cellId1; };
		int& GetCellId2(void) { return this->cellId2;};

	};
}

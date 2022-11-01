#ifndef ___indexer___
#define  ___indexer___
#include "Vector3D.h"

namespace KML{

class CIndexerPrivacy;
class CIndexer{
private:
	CIndexerPrivacy* data;
public:
	CIndexer(void);
	CIndexer(const char* mhdFile);
        CIndexer(string mhdFile);
        CIndexer(const CIndexer& oth);
	~CIndexer();
	void SetOffset(const Vector3D<float>& offset);
	void SetDims( const Vector3D<float>&  dims);
	void SetSize(const Vector3D<size_t>& size);
        Vector3D<size_t> GetIndex(const Vector3D<float>& pos)const;
        Vector3D<size_t> GetIndex(float pos[3])const;
        Vector3D<size_t> GetIndex(double pos[3])const;
        
        Vector3D<float> GetWorldCood(const Vector3D<int>& grid);
        Vector3D<float> GetWorldCood(const Vector3D<size_t>& grid);
        
	const Vector3D<size_t>& GetSize(void) const;
	const Vector3D<float>& GetDims(void) const;
	const Vector3D<float>& GetOffset(void) const;
	void ResetDims( const Vector3D<float>& newDims);
        void ResetDims(float newRes);
	bool operator == (const CIndexer& oth)const;
        bool operator != (const CIndexer& oth)const;
	CIndexer& operator = (const KML::CIndexer& oth);
        size_t GetHistIndex(Vector3D<size_t>& idx)const; 
        size_t GetHistIndex(const Vector3D<float>& pos)const; 
	size_t GetNumGrids(void)const;
        bool IsInitialized(void)const; 
        void OutputIndexerAppend(string fileName)const;
};

        

}
#endif

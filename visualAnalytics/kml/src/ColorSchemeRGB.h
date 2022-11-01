#pragma once
// this class defines several colors that are distinctive and distinguishable; 
// the purpos of this class is to build a visualization scheme for surfac rendering or other work; 
// each color is defined to have RGB values; 

#include "Vector3D.h"
#include <list>
#include <vector>
#include "newran.h"


using namespace std;
using namespace NEWRAN;
using namespace KML;

namespace KML  // included in namespace LKM just like many other libs;
{

	typedef Vector3D<float> RGBTYPE;

	class CColorSchemeRGB
	{
	public:
		CColorSchemeRGB(void);
		~CColorSchemeRGB(void);	
		

	// interface definition 

		void SetColorNums(int numOfColors); // set the number of colors that are needed; 
		unsigned int GetMaxNumOfColors(void) const ; // get the capacity of this color scheme; 
		const RGBTYPE& GetRandomColor(void)  ; // get a random color from all those colors
		const RGBTYPE& GetRandomColorWithExclusion(const RGBTYPE& exclusionColor)  ; // get a random color that is not the exclusion
		const RGBTYPE& GetRandomColorWithExclusionS(const list<RGBTYPE>& exclusions)  ; // get a random color that is not any of the exclusions;
		unsigned int GetRandomIndex(void);
		const RGBTYPE& GetColorByIndex(int index); // get the preset color;
	private:
		void InitializeData(void); 



		//data;
		vector<RGBTYPE> m_colors ;
		unsigned int m_maxNumOfColors ;
		unsigned int m_numOfColorsNeeded;
		Uniform m_uniVariable;


	}; // end of CColorScheme; 

} //end of LKM

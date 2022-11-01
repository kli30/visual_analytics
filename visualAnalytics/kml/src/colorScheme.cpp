/*Copyright 2011 kaiming li ( likaiming@gmail.com) under GPL V3.*/

#include "ColorSchemeRGB.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <ctime> 


using namespace KML;
using namespace std; 

namespace KML
{
	CColorSchemeRGB::CColorSchemeRGB(void)
	{

		// give a random number as seed for this variable using rand() function from stdlib; 
		m_colors.resize(255);
		srand( (unsigned)time( NULL ) );
		double seed= rand() /(1.0+RAND_MAX);
		m_uniVariable.Set(seed);
		InitializeData();
	}
	CColorSchemeRGB::~CColorSchemeRGB()
	{

	}
	void CColorSchemeRGB::InitializeData()
	{

//  		m_colors[7]=RGBTYPE(154,200,237)/255;
// 		m_colors[8]=RGBTYPE(255,128,64)/255;
// 		m_colors[9]=RGBTYPE(99,90,51)/255;
// 		m_colors[10]=RGBTYPE(106,118,249)/255;
// 		m_colors[11]=RGBTYPE(54,244,188)/255;
// 		m_colors[12]=RGBTYPE(198,193,243)/255;
// 		m_colors[13]=RGBTYPE(81,193,61)/255;
//		m_colors[14]=RGBTYPE(0,243,0)/255;
//		m_colors[15]=RGBTYPE(112,233,145)/255;
//		m_colors[16]=RGBTYPE(145,64,126)/255;
//		m_colors[17]=RGBTYPE(121,103,222)/255;
//		m_colors[18]=RGBTYPE(200,255,255)/255;
//		m_colors[19]=RGBTYPE(133,177,24)/255;
// 				m_colors[14]=RGBTYPE(255,255,128)/255;
// 				m_colors[13]=RGBTYPE(251,175,93)/255;
// 				m_colors[15]=RGBTYPE(158,107,82)/255;
// 				m_colors[16]=RGBTYPE(140,98,57)/255;
// 				m_colors[17]=RGBTYPE(184,40,50)/255;
// 				m_colors[18]=RGBTYPE(216,81,113)/255;
// 				m_colors[19]=RGBTYPE(86,63,127)/255;

// 		m_colors[0]=RGBTYPE(131,44,78)/255;
// 		m_colors[1]=RGBTYPE(241,43,17)/255;
// 		m_colors[2]=RGBTYPE(126,31,21)/255;
// 		m_colors[3]=RGBTYPE(70,93,250)/255;
// 		m_colors[4]=RGBTYPE(19,164,195)/255;
// 		m_colors[5]=RGBTYPE(209,38,159)/255;
// 		m_colors[6]=RGBTYPE(88,233,132)/255;
// 		m_colors[7]=RGBTYPE(154,200,237)/255;
// 		m_colors[8]=RGBTYPE(22,172,19)/255;

		m_colors[0]=RGBTYPE(1,0,0);
		m_colors[1]=RGBTYPE(0,1,0);
		m_colors[2]=RGBTYPE(0,0,1);
		m_colors[3]=RGBTYPE(0,1,1);
		m_colors[4]=RGBTYPE(1,1,0);
		m_colors[5]=RGBTYPE(1,0.5,0);
		m_colors[6]=RGBTYPE(1,0,1); 
		m_colors[7]=RGBTYPE(0.5,0,1); 
		m_colors[8]=RGBTYPE(0.5,0,0); 


		m_colors[9]=RGBTYPE(99,90,51)/255;
		m_colors[10]=RGBTYPE(106,118,249)/255;
		m_colors[11]=RGBTYPE(54,244,188)/255;
		m_colors[12]=RGBTYPE(198,193,243)/255;
		m_colors[13]=RGBTYPE(81,193,61)/255;
		m_colors[14]=RGBTYPE(0,243,0)/255;
		m_colors[15]=RGBTYPE(112,233,145)/255;
		m_colors[16]=RGBTYPE(145,64,126)/255;
		m_colors[17]=RGBTYPE(121,103,222)/255;
		m_colors[18]=RGBTYPE(200,255,255)/255;
		m_colors[19]=RGBTYPE(133,177,24)/255;
        	m_colors[20]=RGBTYPE(197,62,87)/255;
		m_colors[21]=RGBTYPE(75,77,226)/255;
		m_colors[22]=RGBTYPE(166,101,172)/255;
		m_colors[23]=RGBTYPE(239,59,213)/255;
		m_colors[24]=RGBTYPE(198,110,171)/255;
		m_colors[25]=RGBTYPE(40,71,34)/255;
		m_colors[26]=RGBTYPE(191,153,35)/255;
		m_colors[27]=RGBTYPE(204,55,143)/255;
		m_colors[28]=RGBTYPE(50,252,63)/255;
		m_colors[29]=RGBTYPE(192,219,228)/255;
		m_colors[30]=RGBTYPE(100,110,32)/255;
		m_colors[31]=RGBTYPE(60,251,166)/255;
		m_colors[32]=RGBTYPE(61,115,201)/255;
		m_colors[33]=RGBTYPE(121,38,62)/255;
		m_colors[34]=RGBTYPE(156,251,121)/255;
		m_colors[35]=RGBTYPE(189,97,122)/255;
		m_colors[36]=RGBTYPE(25,151,88)/255;
		m_colors[37]=RGBTYPE(198,181,113)/255;
		m_colors[38]=RGBTYPE(24,245,140)/255;
		m_colors[39]=RGBTYPE(147,162,199)/255;
		m_colors[40]=RGBTYPE(180,0,0)/255;
		m_colors[41]=RGBTYPE(144,106,78)/255;
		m_colors[42]=RGBTYPE(144,124,154)/255;
		m_colors[43]=RGBTYPE(33,165,9)/255;
		m_colors[44]=RGBTYPE(29,96,164)/255;
		m_colors[45]=RGBTYPE(141,91,144)/255;
		m_colors[46]=RGBTYPE(41,156,43)/255;
		m_colors[47]=RGBTYPE(74,222,212)/255;
		m_colors[48]=RGBTYPE(228,151,137)/255;
		m_colors[49]=RGBTYPE(167,176,67)/255;
		m_colors[50]=RGBTYPE(207,48,107)/255;
		m_colors[51]=RGBTYPE(214,35,166)/255;
		m_colors[52]=RGBTYPE(122,96,128)/255;
		m_colors[53]=RGBTYPE(0,0,90)/255;
		m_colors[54]=RGBTYPE(154,186,142)/255;
		m_colors[55]=RGBTYPE(204,150,68)/255;
		m_colors[56]=RGBTYPE(140,201,226)/255;
		m_colors[57]=RGBTYPE(17,204,231)/255;
		m_colors[58]=RGBTYPE(42,76,42)/255;
		m_colors[59]=RGBTYPE(214,136,9)/255;
		m_colors[60]=RGBTYPE(51,91,158)/255;
		m_colors[61]=RGBTYPE(139,39,209)/255;
		m_colors[62]=RGBTYPE(6,96,157)/255;
		m_colors[63]=RGBTYPE(109,133,95)/255;
		m_colors[64]=RGBTYPE(100,250,152)/255;
		m_colors[65]=RGBTYPE(56,203,222)/255;
		m_colors[66]=RGBTYPE(3,188,106)/255;
		m_colors[67]=RGBTYPE(51,46,199)/255;
		m_colors[68]=RGBTYPE(40,201,10)/255;
		m_colors[69]=RGBTYPE(152,102,58)/255;
		m_colors[70]=RGBTYPE(156,84,154)/255;
		m_colors[71]=RGBTYPE(96,47,76)/255;
		m_colors[72]=RGBTYPE(36,2,225)/255;
		m_colors[73]=RGBTYPE(0,255,255)/255;
		m_colors[74]=RGBTYPE(89,37,27)/255;
		m_colors[75]=RGBTYPE(21,6,241)/255;
		m_colors[76]=RGBTYPE(231,220,37)/255;
		m_colors[77]=RGBTYPE(17,166,187)/255;
		m_colors[78]=RGBTYPE(40,23,31)/255;
		m_colors[79]=RGBTYPE(194,244,7)/255;
		m_colors[80]=RGBTYPE(27,109,78)/255;
		m_colors[81]=RGBTYPE(225,193,129)/255;
		m_colors[82]=RGBTYPE(194,224,127)/255;
		m_colors[83]=RGBTYPE(187,59,13)/255;
		m_colors[84]=RGBTYPE(223,128,172)/255;
		m_colors[85]=RGBTYPE(154,126,150)/255;
		m_colors[86]=RGBTYPE(11,225,27)/255;
		m_colors[87]=RGBTYPE(147,2,98)/255;
		m_colors[88]=RGBTYPE(49,129,197)/255;
		m_colors[89]=RGBTYPE(177,232,171)/255;
		m_colors[90]=RGBTYPE(108,5,54)/255;
		m_colors[91]=RGBTYPE(128,49,164)/255;
		m_colors[92]=RGBTYPE(67,85,179)/255;
		m_colors[93]=RGBTYPE(243,192,222)/255;
		m_colors[94]=RGBTYPE(62,162,216)/255;
		m_colors[95]=RGBTYPE(183,86,12)/255;
		m_colors[96]=RGBTYPE(228,220,0)/255;
		m_colors[97]=RGBTYPE(218,88,31)/255;
		m_colors[98]=RGBTYPE(29,92,52)/255;
		m_colors[99]=RGBTYPE(211,129,208)/255;
		m_colors[100]=RGBTYPE(222,181,164)/255;
		m_colors[101]=RGBTYPE(243,117,132)/255;
		m_colors[102]=RGBTYPE(83,2,249)/255;
		m_colors[103]=RGBTYPE(44,40,136)/255;
		m_colors[104]=RGBTYPE(80,87,40)/255;
		m_colors[105]=RGBTYPE(61,51,183)/255;
		m_colors[106]=RGBTYPE(49,107,197)/255;
		m_colors[107]=RGBTYPE(62,1,131)/255;
		m_colors[108]=RGBTYPE(21,179,66)/255;
		m_colors[109]=RGBTYPE(36,221,77)/255;
		m_colors[110]=RGBTYPE(115,96,118)/255;
		m_colors[111]=RGBTYPE(6,230,50)/255;
		m_colors[112]=RGBTYPE(32,60,153)/255;
		m_colors[113]=RGBTYPE(255,151,47)/255;
		m_colors[114]=RGBTYPE(159,110,145)/255;
		m_colors[115]=RGBTYPE(97,184,130)/255;
		m_colors[116]=RGBTYPE(225,25,168)/255;
		m_colors[117]=RGBTYPE(185,115,214)/255;
		m_colors[118]=RGBTYPE(98,212,159)/255;
		m_colors[119]=RGBTYPE(112,57,62)/255;
		m_colors[120]=RGBTYPE(158,18,33)/255;
		m_colors[121]=RGBTYPE(183,216,248)/255;
		m_colors[122]=RGBTYPE(244,89,147)/255;
		m_colors[123]=RGBTYPE(244,236,110)/255;
		m_colors[124]=RGBTYPE(216,104,74)/255;
		m_colors[125]=RGBTYPE(165,229,151)/255;
		m_colors[126]=RGBTYPE(119,48,32)/255;
		m_colors[127]=RGBTYPE(176,253,185)/255;
		m_colors[128]=RGBTYPE(170,183,96)/255;
		m_colors[129]=RGBTYPE(120,6,240)/255;
		m_colors[130]=RGBTYPE(87,215,33)/255;
		m_colors[131]=RGBTYPE(201,24,192)/255;
		m_colors[132]=RGBTYPE(39,217,243)/255;
		m_colors[133]=RGBTYPE(131,30,221)/255;
		m_colors[134]=RGBTYPE(166,249,246)/255;
		m_colors[135]=RGBTYPE(36,75,134)/255;
		m_colors[136]=RGBTYPE(187,194,76)/255;
		m_colors[137]=RGBTYPE(22,99,31)/255;
		m_colors[138]=RGBTYPE(131,77,231)/255;
		m_colors[139]=RGBTYPE(200,236,200)/255;
		m_colors[140]=RGBTYPE(142,97,198)/255;
		m_colors[141]=RGBTYPE(160,63,207)/255;
		m_colors[142]=RGBTYPE(119,114,126)/255;
		m_colors[143]=RGBTYPE(63,113,99)/255;
		m_colors[144]=RGBTYPE(124,128,153)/255;
		m_colors[145]=RGBTYPE(48,45,191)/255;
		m_colors[146]=RGBTYPE(209,213,236)/255;
		m_colors[147]=RGBTYPE(242,86,175)/255;
		m_colors[148]=RGBTYPE(108,166,63)/255;
		m_colors[149]=RGBTYPE(100,159,107)/255;
		m_colors[150]=RGBTYPE(13,144,104)/255;
		m_colors[151]=RGBTYPE(142,247,226)/255;
		m_colors[152]=RGBTYPE(194,246,193)/255;
		m_colors[153]=RGBTYPE(205,216,91)/255;
		m_colors[154]=RGBTYPE(45,70,141)/255;
		m_colors[155]=RGBTYPE(96,119,115)/255;
		m_colors[156]=RGBTYPE(178,14,47)/255;
		m_colors[157]=RGBTYPE(78,208,135)/255;
		m_colors[158]=RGBTYPE(32,6,189)/255;
		m_colors[159]=RGBTYPE(182,155,114)/255;
		m_colors[160]=RGBTYPE(161,228,240)/255;
		m_colors[161]=RGBTYPE(120,19,179)/255;
		m_colors[162]=RGBTYPE(28,42,126)/255;
		m_colors[163]=RGBTYPE(65,171,238)/255;
		m_colors[164]=RGBTYPE(237,176,195)/255;
		m_colors[165]=RGBTYPE(243,157,209)/255;
		m_colors[166]=RGBTYPE(3,57,221)/255;
		m_colors[167]=RGBTYPE(115,84,85)/255;
		m_colors[168]=RGBTYPE(83,192,27)/255;
		m_colors[169]=RGBTYPE(212,192,97)/255;
		m_colors[170]=RGBTYPE(189,196,205)/255;
		m_colors[171]=RGBTYPE(223,133,207)/255;
		m_colors[172]=RGBTYPE(66,95,127)/255;
		m_colors[173]=RGBTYPE(4,40,207)/255;
		m_colors[174]=RGBTYPE(230,0,0)/255;
		m_colors[175]=RGBTYPE(155,81,167)/255;
		m_colors[176]=RGBTYPE(72,182,237)/255;
		m_colors[177]=RGBTYPE(155,10,215)/255;
		m_colors[178]=RGBTYPE(197,18,222)/255;
		m_colors[179]=RGBTYPE(156,207,160)/255;
		m_colors[180]=RGBTYPE(107,165,60)/255;
		m_colors[181]=RGBTYPE(70,63,224)/255;
		m_colors[182]=RGBTYPE(118,17,192)/255;
		m_colors[183]=RGBTYPE(26,249,93)/255;
		m_colors[184]=RGBTYPE(137,172,114)/255;
		m_colors[185]=RGBTYPE(159,36,96)/255;
		m_colors[186]=RGBTYPE(130,199,113)/255;
		m_colors[187]=RGBTYPE(238,141,191)/255;
		m_colors[188]=RGBTYPE(187,196,46)/255;
		m_colors[189]=RGBTYPE(151,191,58)/255;
		m_colors[190]=RGBTYPE(162,31,252)/255;
		m_colors[191]=RGBTYPE(243,235,67)/255;
		m_colors[192]=RGBTYPE(180,66,6)/255;
		m_colors[193]=RGBTYPE(204,186,31)/255;
		m_colors[194]=RGBTYPE(217,79,152)/255;
		m_colors[195]=RGBTYPE(104,249,101)/255;
		m_colors[196]=RGBTYPE(32,182,224)/255;
		m_colors[197]=RGBTYPE(240,46,234)/255;
		m_colors[198]=RGBTYPE(212,186,250)/255;
		m_colors[199]=RGBTYPE(71,254,253)/255;
		m_colors[200]=RGBTYPE(3,204,96)/255;
		m_colors[201]=RGBTYPE(38,136,12)/255;
		m_colors[202]=RGBTYPE(73,205,193)/255;
		m_colors[203]=RGBTYPE(176,135,163)/255;
		m_colors[204]=RGBTYPE(28,45,55)/255;
		m_colors[205]=RGBTYPE(95,107,216)/255;
		m_colors[206]=RGBTYPE(87,227,209)/255;
		m_colors[207]=RGBTYPE(70,66,99)/255;
		m_colors[208]=RGBTYPE(67,89,144)/255;
		m_colors[209]=RGBTYPE(195,138,249)/255;
		m_colors[210]=RGBTYPE(230,32,179)/255;
		m_colors[211]=RGBTYPE(138,225,96)/255;
		m_colors[212]=RGBTYPE(202,145,108)/255;
		m_colors[213]=RGBTYPE(10,143,193)/255;
		m_colors[214]=RGBTYPE(12,83,244)/255;
		m_colors[215]=RGBTYPE(63,251,184)/255;
		m_colors[216]=RGBTYPE(220,175,254)/255;
		m_colors[217]=RGBTYPE(218,245,147)/255;
		m_colors[218]=RGBTYPE(96,216,137)/255;
		m_colors[219]=RGBTYPE(204,234,93)/255;
		m_colors[220]=RGBTYPE(223,7,19)/255;
		m_colors[221]=RGBTYPE(48,237,61)/255;
		m_colors[222]=RGBTYPE(107,120,44)/255;
		m_colors[223]=RGBTYPE(151,171,57)/255;
		m_colors[224]=RGBTYPE(208,75,105)/255;
		m_colors[225]=RGBTYPE(180,56,90)/255;
		m_colors[226]=RGBTYPE(250,68,186)/255;
		m_colors[227]=RGBTYPE(104,189,50)/255;
		m_colors[228]=RGBTYPE(168,103,23)/255;
		m_colors[229]=RGBTYPE(39,13,183)/255;
		m_colors[230]=RGBTYPE(136,122,152)/255;
		m_colors[231]=RGBTYPE(201,7,142)/255;
		m_colors[232]=RGBTYPE(182,99,231)/255;
		m_colors[233]=RGBTYPE(180,244,169)/255;
		m_colors[234]=RGBTYPE(148,102,18)/255;
		m_colors[235]=RGBTYPE(102,174,79)/255;
		m_colors[236]=RGBTYPE(60,68,187)/255;
		m_colors[237]=RGBTYPE(94,174,119)/255;
		m_colors[238]=RGBTYPE(12,59,180)/255;
		m_colors[239]=RGBTYPE(108,254,148)/255;
		m_colors[240]=RGBTYPE(188,230,232)/255;
		m_colors[241]=RGBTYPE(86,10,178)/255;
		m_colors[242]=RGBTYPE(162,250,181)/255;
		m_colors[243]=RGBTYPE(231,71,251)/255;
		m_colors[244]=RGBTYPE(118,230,155)/255;
		m_colors[245]=RGBTYPE(93,203,155)/255;
		m_colors[246]=RGBTYPE(56,81,180)/255;
		m_colors[247]=RGBTYPE(43,185,1)/255;
		m_colors[248]=RGBTYPE(171,122,52)/255;
		m_colors[249]=RGBTYPE(42,140,128)/255;
		m_colors[250]=RGBTYPE(51,118,36)/255;
		m_colors[251]=RGBTYPE(37,75,90)/255;
		m_colors[252]=RGBTYPE(12,82,170)/255;
		m_colors[253]=RGBTYPE(193,149,27)/255;
		m_colors[254]=RGBTYPE(230,13,149)/255;

//		m_colors[0]=RGBTYPE(220,220,220)/255;
//		m_colors[1]=RGBTYPE(120,120,120)/255;
//		m_colors[2]=RGBTYPE(10,59,118)/255;
//		m_colors[3]=RGBTYPE(67,149,209)/255;
//		m_colors[4]=RGBTYPE(153,217,234)/255;
//		m_colors[5]=RGBTYPE(0,118,163)/255;
//		m_colors[6]=RGBTYPE(13,104,107)/255;
//		m_colors[7]=RGBTYPE(0,169,157)/255;
//		m_colors[8]=RGBTYPE(116,164,2)/255;
//		m_colors[9]=RGBTYPE(132,135,28)/255;
//		m_colors[10]=RGBTYPE(217,213,111)/255;
//		m_colors[11]=RGBTYPE(255,244,104)/255;
//		m_colors[12]=RGBTYPE(255,194,14)/255;
//		m_colors[14]=RGBTYPE(235,97,25)/255;
//		m_colors[13]=RGBTYPE(251,175,93)/255;
//		m_colors[15]=RGBTYPE(158,107,82)/255;
//		m_colors[16]=RGBTYPE(140,98,57)/255;
//		m_colors[17]=RGBTYPE(184,40,50)/255;
//		m_colors[18]=RGBTYPE(216,81,113)/255;
//		m_colors[19]=RGBTYPE(86,63,127)/255;

//		m_colors.push_back(RGBTYPE(0.196078,0,0));
//		m_colors.push_back(RGBTYPE(0,0.392157,0));
//		m_colors.push_back(RGBTYPE(0.588235,0.588235,0.392157));
//		m_colors.push_back(RGBTYPE(0.784314,0.588235,0.392157));
//		m_colors.push_back(RGBTYPE(0,0.784314,0.392157));
//		m_colors.push_back(RGBTYPE(0.196078,0.784314,0.392157));
//		m_colors.push_back(RGBTYPE(0.392157,0.784314,0.392157));
//		m_colors.push_back(RGBTYPE(0.196078,0.392157,0));
//		m_colors.push_back(RGBTYPE(0.588235,0.784314,0.392157));
//		m_colors.push_back(RGBTYPE(0.784314,0.784314,0.392157));
//		m_colors.push_back(RGBTYPE(0,0,0.588235));
//		m_colors.push_back(RGBTYPE(0.196078,0,0.588235));
//		m_colors.push_back(RGBTYPE(0.392157,0,0.588235));
//		m_colors.push_back(RGBTYPE(0.392157,0.392157,0));
//		m_colors.push_back(RGBTYPE(0.588235,0,0.588235));
//		m_colors.push_back(RGBTYPE(0.784314,0,0.588235));
//		m_colors.push_back(RGBTYPE(0,0.196078,0.588235));
//		m_colors.push_back(RGBTYPE(0.196078,0.196078,0.588235));
//		m_colors.push_back(RGBTYPE(0.392157,0.196078,0.588235));
//		m_colors.push_back(RGBTYPE(0.588235,0.196078,0.588235));
//		m_colors.push_back(RGBTYPE(0.588235,0.392157,0));
//		m_colors.push_back(RGBTYPE(0.784314,0.196078,0.588235));
//		m_colors.push_back(RGBTYPE(0,0.392157,0.588235));
//		m_colors.push_back(RGBTYPE(0.784314,0.392157,0));
//		m_colors.push_back(RGBTYPE(0.196078,0.392157,0.588235));
//		m_colors.push_back(RGBTYPE(0.392157,0.392157,0.588235));
//		m_colors.push_back(RGBTYPE(0,0.588235,0));
//		m_colors.push_back(RGBTYPE(0.588235,0.39157,0.588235));
//		m_colors.push_back(RGBTYPE(0.196078,0.588235,0));
//		m_colors.push_back(RGBTYPE(0.784314,0.392157,0.588235));
//		m_colors.push_back(RGBTYPE(0.392157,0.588235,0));
//		m_colors.push_back(RGBTYPE(0.588235,0.588235,0));
//		m_colors.push_back(RGBTYPE(0,0.588235,0.588235));
//		m_colors.push_back(RGBTYPE(0.392157,0,0));
//		m_colors.push_back(RGBTYPE(0.196078,0.588235,0.588235));
//		m_colors.push_back(RGBTYPE(0.784314,0.588235,0));
//		m_colors.push_back(RGBTYPE(0,0.784314,0));
//		m_colors.push_back(RGBTYPE(0.196078,0.784314,0));
//		m_colors.push_back(RGBTYPE(0.392157,0.588235,0.588235));
//		m_colors.push_back(RGBTYPE(0.588235,0.588235,0.588235));
//		m_colors.push_back(RGBTYPE(0.392157,0.784314,0));
//		m_colors.push_back(RGBTYPE(0.588235,0.784314,0));
//		m_colors.push_back(RGBTYPE(0.784314,0.784314,0));
//		m_colors.push_back(RGBTYPE(0.588235,0,0));
//		m_colors.push_back(RGBTYPE(0,0,0.196078));
//		m_colors.push_back(RGBTYPE(0.196078,0,0.196078));
//		m_colors.push_back(RGBTYPE(0.392157,0,0.196078));
//		m_colors.push_back(RGBTYPE(0.588235,0,0.196078));
//		m_colors.push_back(RGBTYPE(0.784314,0,0.196078));
//		m_colors.push_back(RGBTYPE(0,0.196078,0.196078));
//		m_colors.push_back(RGBTYPE(0.196078,0.196078,0.196078));
//		m_colors.push_back(RGBTYPE(0.392157,0.196078,0.196078));
//		m_colors.push_back(RGBTYPE(0.588235,0.196078,0.196078));
//		m_colors.push_back(RGBTYPE(0.784314,0,0));
//		m_colors.push_back(RGBTYPE(0.784314,0.196078,0.196078));
//		m_colors.push_back(RGBTYPE(0,0.392157,0.196078));
//		m_colors.push_back(RGBTYPE(0.196078,0.392157,0.196078));
//		m_colors.push_back(RGBTYPE(0,0.196078,0));
//		m_colors.push_back(RGBTYPE(0.392157,0.392157,0.196078));
//		m_colors.push_back(RGBTYPE(0.588235,0.392157,0.196078));
//		m_colors.push_back(RGBTYPE(0.784314,0.392157,0.196078));
//		m_colors.push_back(RGBTYPE(0,0.588235,0.196078));
//		m_colors.push_back(RGBTYPE(0.196078,0.588235,0.196078));
//		m_colors.push_back(RGBTYPE(0.392157,0.588235,0.196078));
//		m_colors.push_back(RGBTYPE(0.588235,0.588235,0.196078));
//		m_colors.push_back(RGBTYPE(0.196078,0.196078,0));
//		m_colors.push_back(RGBTYPE(0.784314,0.588235,0.196078));
//		m_colors.push_back(RGBTYPE(0,0.784314,0.196078));
//		m_colors.push_back(RGBTYPE(0.196078,0.784314,0.196078));
//		m_colors.push_back(RGBTYPE(0.392157,0.784314,0.196078));
//		m_colors.push_back(RGBTYPE(0.588235,0.784314,0.196078));
//		m_colors.push_back(RGBTYPE(0.784314,0.784314,0.196078));
//		m_colors.push_back(RGBTYPE(0.392157,0.196078,0));
//		m_colors.push_back(RGBTYPE(0,0,0.392157));
//		m_colors.push_back(RGBTYPE(0.196078,0,0.392157));
//		m_colors.push_back(RGBTYPE(0.392157,0,0.392157));
//		m_colors.push_back(RGBTYPE(0.588235,0,0.392157));
//		m_colors.push_back(RGBTYPE(0.784314,0,0.392157));
//		m_colors.push_back(RGBTYPE(0.588235,0.196078,0));
//		m_colors.push_back(RGBTYPE(0,0.196078,0.392157));
//		m_colors.push_back(RGBTYPE(0.196078,0.196078,0.392157));
//		m_colors.push_back(RGBTYPE(0.392157,0.196078,0.392157));
//		m_colors.push_back(RGBTYPE(0.588235,0.196078,0.392157));
//		m_colors.push_back(RGBTYPE(0.784314,0.196078,0.392157));
//		m_colors.push_back(RGBTYPE(0.784314,0.196078,0));
//		m_colors.push_back(RGBTYPE(0,0.392157,0.392157));
//		m_colors.push_back(RGBTYPE(0.196078,0.392157,0.392157));
//		m_colors.push_back(RGBTYPE(0.392157,0.392157,0.392157));
//		m_colors.push_back(RGBTYPE(0.588235,0.392157,0.392157));
//		m_colors.push_back(RGBTYPE(0.784314,0.392157,0.392157));
//		m_colors.push_back(RGBTYPE(0,0.588235,0.392157));
//		m_colors.push_back(RGBTYPE(0.196078,0.588235,0.392157));
//		m_colors.push_back(RGBTYPE(0.392157,0.588235,0.392157));


		m_maxNumOfColors=m_colors.size();
		m_numOfColorsNeeded=m_maxNumOfColors; //default;


	}

	void CColorSchemeRGB::SetColorNums(int numOfColors)
	{
		if (numOfColors>0 && numOfColors<=m_maxNumOfColors)
		{
			m_numOfColorsNeeded=numOfColors;
		}
		else
		{
			cerr<<"ERROR! Number of colors should be between the range: "<<"1~"<<m_maxNumOfColors<<endl;
			exit(EXIT_FAILURE);			
		}
	}

	unsigned int CColorSchemeRGB::GetMaxNumOfColors()const
	{
		return m_maxNumOfColors;
	}

    const RGBTYPE & CColorSchemeRGB::GetColorByIndex(int index)
    {
    	if(index>=m_maxNumOfColors || index <0 )
    	{
    		cerr<<"Invalid index: "<<index<<endl;
    		cerr<<" should be 0~"<<m_maxNumOfColors<<endl;
    		exit(EXIT_FAILURE);
    	}
    	return m_colors[index];
    }
	unsigned int CColorSchemeRGB::GetRandomIndex() 
	{
		unsigned int index=0;

		double randomValue=0; 		
		randomValue=m_uniVariable.Next();
		randomValue*=m_numOfColorsNeeded; // generate range 0~m_numOfColorsNeeded

		index=floor(randomValue);  //this will generate 0,1,2,...,m_numOfColorsNeeded-1;
		return index;

	}
	const RGBTYPE& CColorSchemeRGB::GetRandomColor() 
	{
		unsigned int index=GetRandomIndex();
		return m_colors[index];
	}


	const RGBTYPE& CColorSchemeRGB::GetRandomColorWithExclusion(const KML::RGBTYPE &exclusionColor) 
	{

	
		unsigned int index=GetRandomIndex();
		while (exclusionColor== m_colors[index])
		{
			index=GetRandomIndex();
		}
		
		return m_colors[index];
	}
	
	const RGBTYPE& CColorSchemeRGB::GetRandomColorWithExclusionS(const std::list<RGBTYPE> &exclusions) 
	{
		unsigned int index=GetRandomIndex();
		std::list<RGBTYPE>::const_iterator itResult;

		do 
		{
			index=GetRandomIndex();
			itResult=find(exclusions.begin(),exclusions.end(),m_colors[index]);			
		}
		while(itResult!=exclusions.end());

		return m_colors[index];
	}



}

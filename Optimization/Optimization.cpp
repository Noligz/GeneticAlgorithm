#include <stdio.h>
#include <math.h>

#include <ga/GASimpleGA.h>
#include <ga/GABin2DecGenome.h>
#include <ga/std_stream.h>

#include "Parameters.h"

//#define cout STD_COUT
#include <iostream>
#include <string>
using namespace std;

QiKuai qiKuaiList;
Gang gangList;
GuangKongLv guangKongLvList;

const double PI = 3.1415926;

const double factor = pow(10, 15);

void InitParameterList()
{
	qiKuaiList.push_back(QiKuaiComponent("MU10", 2.5, 4, "CB20", 9.6, 5));
	qiKuaiList.push_back(QiKuaiComponent("MU15", 2.79, 5, "CB30", 10.2, 6));
	qiKuaiList.push_back(QiKuaiComponent("MU20", 3.1, 6, "CB40", 14.1, 7));

	gangList.push_back(GangComponent("HPB235", 210, 2, 0.6));
	gangList.push_back(GangComponent("HRB335", 300, 3, 0.53));

	guangKongLvList.push_back(1.0 / 3.0);
	guangKongLvList.push_back(0.5);
	guangKongLvList.push_back(1.0);
}


float Objective(GAGenome &);

void main()
{
	InitParameterList();

	unsigned int seed = 0;
	GARandomSeed(seed);

	//砌块，墙体长度，灌孔率，钢筋类型，配筋率，纵向边缘钢筋直径，水平钢筋直径，水平钢筋间距
	float paraMin [] = { 0.000000, 1200, 0.000000, 0.000000, 0.007, 8, 8, 100 };
	float paraMax [] = { 2.999999, 5000, 2.999999, 1.999999, 0.06, 25, 25, 1000 };
	unsigned int parameterCount = sizeof(paraMin) / sizeof(float);

	GABin2DecPhenotype map;
	for (int i = 0; i < parameterCount; i++)
		map.add(8, paraMin[i], paraMax[i]);

	GABin2DecGenome genome(map, Objective);

	int popsize = 150;
	int ngen = 200;
	float pmut = 0.01;
	float pcross = 0.6;

	GASimpleGA ga(genome);
	ga.populationSize(popsize);
	ga.nGenerations(ngen);
	ga.pMutation(pmut);
	ga.pCrossover(pcross);
	//ga.scoreFilename("bog.dat");
	ga.flushFrequency(50);	// dump scores to disk every 50th generation
	ga.evolve(seed);

	genome.initialize();

	genome = ga.statistics().bestIndividual();
	cout << "运算结果:\n";

	int qiKuaiIdx = genome.phenotype(0);
	QiKuaiComponent qiKuai = qiKuaiList[qiKuaiIdx];
	cout << "砌块: " << qiKuai.qiKuai.name << endl;

	cout << "墙体长度: " << genome.phenotype(1) << endl;

	int guangKongLvIdx = genome.phenotype(2);
	cout << "灌孔率: " << guangKongLvList[guangKongLvIdx] << endl;

	int gangIdx = genome.phenotype(3);
	GangComponent gang = gangList[gangIdx];
	cout << "钢筋类型: " << gang.gang.name << endl;

	cout << "配筋率: " << genome.phenotype(4) << endl;

	cout << "纵向边缘钢筋直径: " << genome.phenotype(5) << endl;

	cout << "水平钢筋直径: " << genome.phenotype(6) << endl;

	cout << "水平钢筋间距: " << genome.phenotype(7) << endl;

	cout << "总造价: " << factor / genome.score() << endl;

	cout << "\n\n"; cout.flush();

	system("pause");
}

bool Restrict(GABin2DecGenome & genome)
{
	//0砌块，1墙体长度，2灌孔率，3钢筋类型，4配筋率，5纵向边缘钢筋直径，6水平钢筋直径，7水平钢筋间距
	int qiKuaiIdx = genome.phenotype(0);
	QiKuaiComponent qiKuai = qiKuaiList[qiKuaiIdx];
	float f = qiKuai.qiKuai.strength;
	float fc = qiKuai.hunNingTu.strength;

	float h = genome.phenotype(1);

	int guangKongLvIdx = genome.phenotype(2);
	float gkl = guangKongLvList[guangKongLvIdx];

	int gangIdx = genome.phenotype(3);
	GangComponent gang = gangList[gangIdx];
	float fy = gang.gang.strength;
	float ipsironB = gang.ipsiron;

	float pw = genome.phenotype(4);

	float d = genome.phenotype(5);

	float ds = genome.phenotype(6);

	float l = genome.phenotype(7);

	float fg = f + 0.6*0.46*gkl*fc;
	float fvg = 0.2*pow(fg, 0.55);
	float en = 17 / 16 * 1000 + 4500 * 4500 / h / 2200 * (1 - 0.022 * 4500 / h);
	float x = (1600 * 1000 + fy*pw * 190 * (h - 300)) / ((fg + 1.5*fy*pw) * 190);

	if (x < 600)
		return false;
	float tmp = fg * 190 * x*(h - 300 - x / 2) + 3 / 4 * PI*d*d*fy*(h - 600);
	if (x < ipsironB*(h - 300))
	{
		if (1600 * 1000 * en > tmp - fy*pw*pow((h - 300 - x / 2), 2) / 2)
			return false;
	}
	else
	{
		if (1600 * 1000 * en > tmp)
			return false;
	}

	float lumda = 17 / 4 * 1000 / (h - 300);
	if (lumda < 1.5)
		lumda = 1.5;
	if (lumda > 2.2)
		lumda = 2.2;
	float ph = 2 * PI / 4 * ds*ds / 190 / l;
	if (ph < 0.0007)
		return false;
	else
	{
		if (0.25*fg * 190 * h < 1600)
		{
			if ((0.6*fvg * 190 * (h - 300) + 0.12*0.25*fg * 190 * h * 1000) / (lumda - 0.5) + 0.9*fy * 2 * PI / 4 * ds*ds / l*(h - 300) < 400)
				return false;
		}
		else
		{
			if ((0.6*fvg * 190 * (h - 300) + 0.12 * 1600 * 1000) / (lumda - 0.5) + 0.9*fy * 2 * PI / 4 * ds*ds / l*(h - 300) < 400)
				return false;
		}
	}

	if (0.25*fg * 190 * h < 400)
		return false;

	return true;
}

float Objective(GAGenome& g)
{
	GABin2DecGenome & genome = (GABin2DecGenome &) g;
	if (!Restrict(genome))
		return 0;

	//0砌块，1墙体长度，2灌孔率，3钢筋类型，4配筋率，5纵向边缘钢筋直径，6水平钢筋直径，7水平钢筋间距
	int qiKuaiIdx = genome.phenotype(0);
	QiKuaiComponent qiKuai = qiKuaiList[qiKuaiIdx];
	float k1 = qiKuai.qiKuai.price;
	float k2 = qiKuai.hunNingTu.price;

	float h = genome.phenotype(1);

	int guangKongLvIdx = genome.phenotype(2);
	float gkl = guangKongLvList[guangKongLvIdx];

	int gangIdx = genome.phenotype(3);
	GangComponent gang = gangList[gangIdx];
	float k3 = gang.gang.price;

	float pw = genome.phenotype(4);

	float ds = genome.phenotype(6);

	float l = genome.phenotype(7);

	float p1 = k1*4.5*0.19*h / 1000;
	float p2 = k2*4.5*0.19*h*gkl*0.46 / 1000;
	float p3 = k3* pw * 190 * 200 / gkl*h*gkl * 200 * 4500;
	int n = 4500 / l;
	float p4 = k3 * 2 * h*PI / 4 * ds*ds*n;
	float K = p1 + p2 + p3 + p4;

	return factor / K;
}

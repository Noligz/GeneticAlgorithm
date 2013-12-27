#include <vector>

class Component
{
public:
	std::string name;
	float strength;
	float price;
public:
	Component(std::string n, float s, float p)
		:name(n), strength(s), price(p){}
};

class GangComponent
{
public:
	Component gang;
	float ipsiron;
public:
	GangComponent(std::string n, float s, float p, float i)
		:gang(n, s, p), ipsiron(i){}
};
typedef std::vector<GangComponent> Gang;

class QiKuaiComponent
{
public:
	Component qiKuai;
	Component hunNingTu;
public:
	QiKuaiComponent(std::string nq, float sq, float pq, std::string nh, float sh, float ph)
		:qiKuai(nq, sq, pq), hunNingTu(nh, sh, ph){}
};
typedef std::vector<QiKuaiComponent> QiKuai;

typedef std::vector<float> GuangKongLv;
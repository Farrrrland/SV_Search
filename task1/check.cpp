#include<bits/stdc++.h>
using namespace std;
struct node
{
	string s, n1, n2;
	int pos11, pos12, pos21, pos22;
}sv[55], tp;
bool ans[55];
inline bool ck(node a, node b)
{
	if(a.s != b.s)
		return 0;
	if(a.s != "TRA")
	{
		if(a.n1 != b.n1)
			return 0;
		return 1.0 * max(0, min(a.pos12, b.pos12) - max(a.pos11, b.pos11)) / max(a.pos12 - a.pos11, b.pos12 - b.pos11) > 0.6;
	}
	if(a.n1 != b.n1)
		swap(b.n1, b.n2), swap(b.pos11, b.pos21), swap(b.pos12, b.pos22);
	if(a.n1 != b.n1 || a.n2 != b.n2)
		return 0;
	if(1.0 * max(0, min(a.pos12, b.pos12) - max(a.pos11, b.pos11)) / max(a.pos12 - a.pos11, b.pos12 - b.pos11) > 0.6)
		return 1.0 * max(0, min(a.pos22, b.pos22) - max(a.pos21, b.pos21)) / max(a.pos22 - a.pos21, b.pos22 - b.pos21) > 0.6;
	return 0;
	
}
int main()
{
	ifstream in;
	in.open("task1_svs.bed"); //open standard answer file
	int cnt = 50;
	for(int i = 0; i < cnt; i++)
	{
		in >> sv[i].s;
		in >> sv[i].n1 >> sv[i].pos11 >> sv[i].pos12;
		if(sv[i].s == "TRA")
			in >> sv[i].n2 >> sv[i].pos21 >> sv[i].pos22;
	}
	in.close();
	in.open("sv.bed"); //open your answer file
	int num = 0;
	while(num++ < 120)
	{
		in >> tp.s;
		in >> tp.n1 >> tp.pos11 >> tp.pos12;
		if(tp.s == "TRA")
			in >> tp.n2 >> tp.pos21 >> tp.pos22;
		if(in.eof())
			break;
		for(int i = 0; i < cnt; i++)
			if(ck(sv[i], tp))
				ans[i] = 1;	
	}
	int r = 0;
	for(int i = 0; i < cnt; i++)
		r += ans[i];
	cout << "the score is " << r;
} 

#include "HighPrecisionTimer.hpp"
#include <iostream>

using namespace Panda;
using namespace std;

int64_t time1, time2;
HighPrecisionTimer timer;

void TestTime1()
{
	int32_t a(0), b(1), c(2), d(3), e(4), f(5), g(6), h(7), i(8), j(9), k(10), l(11), m(12), n(13);
	int32_t o(14), p(15), q(16), r(17), s(18), t(19), u(20), v(21), w(22), x(23), y(24), z(25);
	int32_t aa, bb, cc, dd, ee, ff, gg, hh, ii, jj, kk, ll, mm, nn, oo, pp, qq, rr, ss, tt, uu, vv, ww, xx, yy, zz;

	timer.Start();
	aa = a + b;
	bb = aa + c;
	cc = bb + d;
	dd = cc + e;
	ee = dd + f;
	ff = ee + g;
	gg = ff + h;
	hh = gg + g;
	ii = hh + h;
	jj = ii + i;
	kk = jj + j;
	ll = kk + k;
	mm = ll + l;
	nn = mm + m;
	oo = nn + n;
	pp = oo + o;
	qq = pp + p;
	rr = qq + q;
	ss = rr + r;
	tt = ss + s;
	uu = tt + t;
	vv = uu + u;
	ww = vv + v;
	xx = ww + w;
	yy = xx + x;
	zz = yy + y;
	timer.Stop();
	time1 = timer.TotalClockCount();
}

void TestTime2()
{
	int32_t a(0), b(1), c(2), d(3), e(4), f(5), g(6), h(7), i(8), j(9), k(10), l(11), m(12), n(13);
	int32_t o(14), p(15), q(16), r(17), s(18), t(19), u(20), v(21), w(22), x(23), y(24), z(25);
	int32_t aa, bb, cc, dd, ee, ff, gg, hh, ii, jj, kk, ll, mm, nn, oo, pp, qq, rr, ss, tt, uu, vv, ww, xx, yy, zz;

	timer.Start();
	aa = a + b;
	bb = c + d;
	cc = e + f;
	dd = g + h;
	ee = i + j;
	ff = k + l;
	gg = m + n;
	hh = o + p;
	ii = o + s;
	jj = o + z;
	kk = a + z;
	ll = b + m;
	mm = w + y;
	nn = i + y;
	oo = q + o;
	pp = w + x;
	qq = o + h;
	rr = l + a;
	ss = z + s;
	tt = a + t;
	uu = u + u;
	vv = v + w;
	ww = o + i;
	xx = w + e;
	yy = o + h;
	zz = h + g;
	timer.Stop();
	time2 = timer.TotalClockCount();
}

int main()
{
	TestTime1();
	TestTime2();

    cout << "timer1 = " << time1 << endl;
    cout << "timer2 = " << time2 << endl;
}
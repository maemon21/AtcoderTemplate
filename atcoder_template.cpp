#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <string>
#include <sstream>
#include <complex>
#include <vector>
#include <list>
#include <queue>
#include <deque>
#include <stack>
#include <map>
#include <set>
#include <functional>
#include <iomanip>
#include <limits>
//#include <bits/stdc++.h>

using namespace std;
typedef long long ll; //10^18
typedef unsigned long ul;

using Graph = vector<vector<ll>>;

typedef map<ll, ll> m;
typedef multimap<ll, ll> mm;
typedef set<ll> s;
typedef multiset<ll> ms;
typedef priority_queue<ll> pq;
typedef queue<ll> q;
typedef deque<ll> dq;
typedef list<ll> lst;
typedef pair<ll, ll> p;

#define EPS (1e-7)
#define INF (1e9)
#define PI (acos(-1))
#define MOD 1000000007LL
#define WALL '#'

//#define and &&a
//#define or ||
//#define not !
//#define neq !=
//#define eq ==

#define REP(i, n) for (ll i = 0; i < n; i++)	   // from 0 to n
#define REPR(i, n) for (ll i = n; i >= 0; i--)	   // from n to 0
#define FOR(i, m, n) for (ll i = m; i < n; i++)   // from m to n
#define FORR(i, m, n) for (ll i = m; i >= n; i--) // from m to n
#define DBG(a) cout << #a << " : " << a << "\n";
#define MSG(a) cout << a << "\n";
#define ALL(v) v.begin(), v.end()
#define SZ(x) ((int)(x).size())
#define PNT(a) printf("%lld", (a))

#define pb push_back //配列などの最後に要素を追加
#define mp make_pair
#define lb lower_bound
#define ub upper_bound
#define FST first
#define SND second

template <class T>
bool chmax(T& a, const T& b)
{
	if (a < b)
	{
		a = b;
		return 1;
	}
	return 0;
}
template <class T>
bool chmin(T& a, const T& b)
{
	if (b < a)
	{
		a = b;
		return 1;
	}
	return 0;
}

vector<ll> dy = { 0, 0, -1, 1, 0 };
vector<ll> dx = { -1, 1, 0, 0, 0 };

//swap(a, b);
//sort(arr, arr + n);	//昇順
//sort(arr, arr+n, greater<int>());	//降順
//max(a, b);
//min(a, b);

//upper_bound(a, a+n, k)	//配列aの中で、kより大きい値が初めて現れる位置へのポインタ
//upper_bound(ALL(v), k)	//STLvの中で、kより大きい値が初めて現れる位置へのポインタ
//lower_bound(a, a+n, k)
//lower_bound(ALL(v), k)	//STLvの中で、kの以上値が初めて現れる位置へのポインタ
//lower_bound(ALL(v),k) - upper_bound(ALL(v),k)	//二分探索を用いて、ある列aに含まれる数kの個数を求める

// n個のデータをvectorで取得
vector<ll> INV(ll n)
{
	vector<ll> v(n);
	REP(i, n)
		cin >> v[i];
	return v;
}

// n個のデータをvectorで取得
vector<vector<ll>> INV2(ll n, ll m)
{
	vector<vector<ll>> v(n, vector<ll>(m));
	REP(i, n)
	{
		REP(j, m)
		{
			cin >> v[i][j];
		}
	}
	return v;
}

// index が条件を満たすかどうか
bool isOK(vector<ll>& v, int index, int key)
{
	if (v[index] >= key)
		return true;
	else
		return false;
}

// 汎用的な二分探索
// sort
ll bs(vector<ll>& v, ll key)
{
	int ng = -1;	//「index = 0」が条件を満たすこともあるので、初期値は -1
	int ok = SZ(v); // 「index = a.size()-1」が条件を満たさないこともあるので、初期値は a.size()

	/* ok と ng のどちらが大きいかわからないことを考慮 */
	while (abs(ok - ng) > 1)
	{
		int mid = (ok + ng) / 2;

		if (isOK(v, mid, key))
			ok = mid;
		else
			ng = mid;
	}
	return ok;
}

// 最大公約数
ll gcd(ll a, ll b)
{
	if (a < b)
		swap(a, b);
	ll r = a % b;
	while (r != 0)
	{
		a = b;
		b = r;
		r = a % b;
	}
	return b;
}

// 最小公倍数
void lcm(ll a, ll b) {

}

// 素数判定
bool is_prime(ll n)
{
	bool flg = true;
	if (n <= 1)
		flg = false;
	else if (n == 2)
		flg = true;
	else if (n % 2 == 0)
		flg = false;
	else
	{
		for (ll i = 3; i * i <= n; i += 2)
		{
			if (n % i == 0)
				flg = false;
		}
	}
	return flg;
}

// 素因数分解
// iで割った回数をcnt_pf[i - 1]に格納している
// cnt_pf[0]に入力が素数だった場合にその素数が入る
vector<ll> prime_factorization(ll n)
{
	vector<ll> cnt_pf(sqrt(n), 0);
	FOR(i, 1, SZ(cnt_pf))
	{
		while (n % (i + 1) == 0)
		{
			cnt_pf[i]++;
			n /= (i + 1);
		}
		if (n == 1)
			break;
	}
	if (n != 1)
	{
		cnt_pf[0] = n;
	}
	return cnt_pf;
}

// 迷路のマップ情報をベクトル化する
// 通れるところを0に、壁を-1にする
// スタート地点からの距離を格納するときなどに使う
vector<vector<ll>> map_vec(vector<string>& str)
{
	// SZ(str[0] = SZ(str[distance(str.begin(), max_element(ALL(str)))])
	vector<vector<ll>> v(SZ(str), vector<ll>(SZ(str[distance(str.begin(), max_element(ALL(str)))]), (int)INF));
	REP(i, SZ(str))
	{
		REP(j, SZ(str[i]))
		{
			if (str[i][j] == WALL)
				v[i][j] = -1;
			//else	v[i][j] = INF;	// if (str[i][j] == '.')
		}
	}
	return v;
}
// str中のWALL='#'の数を数える
ll cnt_wall(vector<string> str)
{
	ll cnt = 0;
	REP(i, SZ(str))
	{
		REP(j, SZ(str[i]))
		{
			if (str[i][j] == WALL)
				cnt++;
		}
	}
	return cnt;
}

// 迷路用幅優先探索
// フィールドの広さと壁の位置を受け取り、ゴールへの最短距離を返す
ll bfs_maze(vector<string>& str, ll s_y, ll s_x, ll g_y, ll g_x)
{
	struct Corr
	{
		ll y;
		ll x;
		ll depth;
	};
	queue<Corr> que;
	// SZ(str[0] = SZ(str[distance(str.begin(), max_element(ALL(str)))])
	vector<vector<ll>> v(SZ(str), vector<ll>(SZ(str[distance(str.begin(), max_element(ALL(str)))])));
	v = map_vec(str);

	que.push({ s_y, s_x, 1 });
	while (!que.empty())
	{
		Corr now = que.front();
		que.pop();
		if (now.y == g_y && now.x == g_x)
			break;

		REP(i, 4)
		{
			Corr next = { now.y + dy[i], now.x + dx[i], now.depth + 1 };
			// SZ(v[0] = SZ(v[distance(v.begin(), max_element(ALL(v)))])
			if (0 <= (int)next.y && (int)next.y < SZ(v) && 0 <= (int)next.x && (int)next.x < SZ(v[distance(v.begin(), max_element(ALL(v)))]) && v[(int)next.y][(int)next.x] == INF)
			{
				v[(int)next.y][(int)next.x] = next.depth;
				que.push(next);
			}
		}
	}

	return v[(int)g_y][(int)g_x];
}

// 累積和
vector<ll> cumulative_sum(vector<ll> a)
{
	vector<ll> v(SZ(a) + 1);
	v[0] = 0;
	REP(i, SZ(a))
	{
		v[i + 1] = v[i] + a[i];
	}
	return v;
}

// 繰り返し二乗法 a^n
ll iterativePower(ll a, ll n) {
	ll res = 1;
	while (n > 0) {
		if (n & 1)	res = (res * a) % MOD;
		a = a * a % MOD;
		n >>= 1;
	}
	return res;
}

// 階乗
ll factorial(ll n) {
	ll res = 1;

	REP(i, n) {
		res = res * (i + 1) % MOD;
	}
	return res;
}

// MODの逆元
vector<ll> MODInverse(ll n, ll fact) {
	vector<ll> res(n + 1);
	res[n] = iterativePower(fact, MOD - 2);
	REPR(i, n - 1) {
		res[i] = res[i + 1] * (i + 1) % MOD;
	}
	return res;
}

// 二項係数nCr
ll comb(ll n, ll r) {
	ll fact = factorial(n);
	vector<ll> fact_inv;
	fact_inv = MODInverse(n, fact);
	return (fact * fact_inv[r]) % MOD * fact_inv[n - r] % MOD;
}


// 文字列を連続した文字ごとに分解
vector<pair<char, ll>> decompose_str(string s)
{
	vector<pair<char, ll>> moji_cnt;
	moji_cnt.pb(mp(s[0], 0));
	REP(i, SZ(s))
	{
		if (moji_cnt.back().first == s[i])
		{
			moji_cnt.back().second++;
		}
		else
		{
			moji_cnt.pb(mp(s[i], 1));
		}
	}
	return moji_cnt;
}

// 解答のベクトル出力(空白区切り)
void ans_vec(vector<ll> ans)
{
	REP(i, SZ(ans))
	{
		cout << ans[i] << " ";
	}
	cout << endl;
}

//
void dinamic_programming(void)
{
}

// 総和の公式：Σk
ll totalSumFirst(ll x, ll y)
{
	x = x - 1;
	return (y * (y + 1) - x * (x + 1)) / 2;
}

// 総和の公式：Σk^2
ll totalSumSecond(ll x, ll y)
{
	x = x - 1;
	return (y * (y + 1) * (2 * y + 1) - x * (x + 1) * (2 * x + 1)) / 6;
}

// 総和の公式：Σk^3
ll totalSumThird(ll x, ll y)
{
	return pow(totalSumFirst(x, y), 2);
}

// 約数
vector<ll> makeDivisors(ll n)
{
	vector<ll> divisors;
	for (ll i = 1; i * i <= n; i++)
	{
		if (n % i == 0)
		{
			divisors.pb(i);
			if (i != n / i)
			{
				divisors.pb(n / i);
			}
		}
	}
	sort(ALL(divisors));

	return divisors;
}

// 尺取り法
ll shakutori(vector<ll>& v, ll x) {
	ll res = 0;
	ll n = SZ(v);

	ll sum = 0;
	ll right = 0;

	REP(left, n) {

		while (right < n && sum + v[right] <= x) {
			sum += v[right];
			right++;
		}
		res += (right - left);

		if (right == left)	right++;
		else sum -= v[left];
	}

	return res;
}



// 深さ優先探索
void dfs(const Graph& g, ll x) {

}

// 幅優先探索
void bfs() {

}

//
// main関数
//

signed main()
{
	//cin.tie(0);
	//ios::sync_with_stdio(false);

	// 変数（scala）取得
	ll w, h, n;
	cin >> w >> h >> n;

	// 変数（vector）取得
	 //vector<ll> a(n);
	 //a = INV(n);
	////m=2;
	//vector<vector<ll>> vec(n, vector<ll>(m));
	//vec = INV2(n, m);

	// vector<vector<ll>> a(3);
	// REP(i, a)
	// {
	// 	ll tmp;
	// 	cin >> tmp;
	// 	a[0].pb(tmp);
	// }
	// REP(i, b)
	// {
	// 	ll tmp;
	// 	cin >> tmp;
	// 	a[1].pb(tmp);
	// }
	// REP(i, c)
	// {
	// 	ll tmp;
	// 	cin >> tmp;
	// 	a[2].pb(tmp);
	// }

	// 文字列取得
	// string s;
	// cin >> s;

	// 文字列（vector）取得
	// vector<string> str(n);
	// REP(i, n)
	// {
	// 	cin >> str[i];
	// }

	// グラフ取得
	//Graph g(n);
	//REP(i, k) {
	//	ll from, weight;
	//	cin >> from >> weight;
	//	g[from].pb(to);
	//	//g[to].pb(from);
	//}


	//
	// 実装部分
	//

	ll ans = 0;


	ll w_min, w_max, h_min, h_max;
	w_min = h_min = 0;
	w_max = w;
	h_max = h;

	REP(i, n) {
		ll x, y, a;
		cin >> x >> y >> a;

		if (a == 1) {
			w_min = max(w_min, x);
		}
		else if (a == 2) {
			w_max = min(w_max, x);
		}
		else if (a == 3) {
			h_min = max(h_min, y);
		}
		else if (a == 4) {
			h_max = min(h_max, y);
		}
	}

	

	ll tmp = (h_max - h_min) * (w_max - w_min);
	if ((h_max - h_min) >= 0 && (w_max - w_min) >= 0) ans = tmp;



	//
	// 実装部分おわり
	//

	// 解答出力
	//cout << fixed << setprecision(10);
	MSG(ans);
	//ans_vec(ans);

	return 0;
}

//
// memo
//

//for(ll i=0; i<n;i++)　// ループ
//cout << << endl;	// 出力
//sort(ALL(a), greater<ll>()); // 降順
// abs(k)	// 絶対値
// sqrt(n)	// ルート
// ceil(n)	// 切り上げ
// floor(n)	// 切り捨て
// round(n)	// 四捨五入

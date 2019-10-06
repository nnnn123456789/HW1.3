#include <memory>
#include <functional>
#include <cstring>
//#include <tuple>
#include <cmath>

#pragma once
#pragma pack(push)
#pragma pack(4)
#pragma warning(disable:4634)
enum class exception
{

};
enum class matrix_exception
{
	matrix_strange
};

template<typename _Data, int n> class vector;
template<typename _Data, int n> class matrix;



template<typename _Data, int n>
class vector
{

protected:
	using _MyType = vector<_Data, n>;
	constexpr static int length = n;
	std::unique_ptr<_Data> head;
private:
	_MyType& me = *this;
protected:
	vector(std::nullptr_t para) : head(nullptr)
	{

	}
public:
	vector() : head(new _Data[length])
	{

	}
	vector(std::initializer_list<_Data> l) : vector()
	{
		for (long int i = 0; i < l.size() && i < length; i++)
		{
			head.get()[i] = l.begin()[i];
		}
	}
	vector(const _MyType& r) : head(new _Data[length])
	{
		std::memcpy(this->head.get(), r.head.get(), length * sizeof(_Data));
	}
	vector(_MyType&& r) : head(nullptr)
	{
		this->head.swap(r.head);
	}
	_MyType& operator=(const _MyType& r)
	{
		std::memcpy(this->head.get(), r.head.get(), length * sizeof(_Data));
	}
	_MyType& operator=(_MyType&& r)
	{
		this->head.swap(r.head);
	}
	virtual _Data& operator[](const int n)
	{
		return head.get()[n];
	}
	virtual _Data& locate(const int n)
	{
		return head.get()[n];
	}
	virtual const _Data& operator[](const int n) const
	{
		return head.get()[n];
	}
	virtual const _Data& locate(const int n) const
	{
		return head.get()[n];
	}
	_MyType operator+(const _MyType& r)const
	{
		const _MyType& l = *this;
		_MyType ret;
		for (int i = 0; i < length; i++)
		{
			ret[i] = l[i] + r[i];
		}
		return ret;
	}
	_MyType operator-(const _MyType& r)const
	{
		const _MyType& l = *this;
		_MyType ret;
		for (int i = 0; i < length; i++)
		{
			ret[i] = l[i] - r[i];
		}
		return ret;
	}
	virtual void init(std::function<_Data(int)> fun)
	{
		for (int i = 0; i < length; i++)
		{
			me[i] = fun(i);
		}
	}
	friend class matrix<_Data, n>;


};



template<typename _Data, int n>
class matrix : public vector<_Data, n* n>
{
protected:
	using _MyBase = vector<_Data, n* n>;
	using _MyType = matrix< _Data, n>;
	using _VecType = vector<_Data, n>;
	constexpr static int length = n;
	constexpr static int size = length * length;
	std::unique_ptr<_Data> head;
private:
	_MyType& me = *this;

public:
	matrix() : _MyBase() //head(new _Data[size])
	{

	}
	matrix(std::initializer_list<_Data> l) : matrix()
	{
		for (long int i = 0; i < l.size() && i < size; i++)
		{
			head.get()[i] = l.begin()[i];
		}
	}
	matrix(const _MyType& r) : matrix()
	{
		std::memcpy(this->head.get(), r.head.get(), size * sizeof(_Data));
	}
	matrix(_MyType&& r) : _MyBase(nullptr)
	{
		this->head.swap(r.head);
	}
	_MyType& operator=(const _MyType& r)
	{
		std::memcpy(me.head.get(), r.head.get(), size * sizeof(_Data));
		return me;
	}
	_MyType& operator=(_MyType&& r)
	{
		me.head.swap(r.head);
		return me;
	}
protected:
	_Data& operator[](const int n)
	{
		return head.get()[n];
	}
	const _Data& operator[](const int n) const
	{
		return head.get()[n];
	}
public:
	_Data& locate(const int x, const int y)
	{
		return head.get()[x * length + y];
	}
	const _Data& locate(const int x, const int y) const
	{
		return head.get()[x * length + y];
	}
	_MyType operator+(const _MyType& r)const
	{
		const _MyType& l = *this;
		_MyType ret;
		for (int i = 0; i < size; i++)
		{
			ret[i] = l[i] + r[i];
		}
		return ret;
	}
	_MyType operator-(const _MyType& r)const
	{
		const _MyType& l = *this;
		_MyType ret;
		for (int i = 0; i < size; i++)
		{
			ret[i] = l[i] - r[i];
		}
		return ret;
	}
	///<summary>将矩阵按照所给函数初始化</summary>
	///<param name="fun">初始化矩阵使用的函数</param>
	void init(std::function<_Data(int, int)> fun)
	{
		for (int i = 0; i < length; i++)
			for (int j = 0; j < length; j++)
			{
				me.locate(i, j) = fun(i, j);
			}
	}
	_VecType operator*(const _VecType& v) const
	{
		_VecType ret;
		for (int i = 0; i < length; i++)
		{
			ret[i] = 0;
			for (int j = 0; j < length; j++)
			{
				ret[i] += v[j] * me.locate(i, j);
			}
		}
		return ret;
	}
	_MyType operator*(const _MyType& r) const
	{
		_MyType ret;
		for (int i = 0; i < length; i++)
			for (int k = 0; k < length; k++)
			{
				ret.locate(i, k) = 0;
				for (int j = 0; j < length; j++)
				{
					ret.locate(i, k) += me.locate(i, j) * r.locate(j, k);
				}
			}
		return ret;
	}
	//原地转置
	void dotrans()
	{
		for (int i = 0; i < length; i++)
			for (int j = 0; j < i; j++)
			{
				auto temp = me.locate(i, j);
				me.locate(i, j) = me.locate(j, i);
				me.locate(j, i) = temp;
			}
	}

	///<summary>计算矩阵的转置但不改变原值</summary>
	///<return>返回转置</return>
	_MyType trans() const
	{
		_MyType ret;
		for (int i = 0; i < length; i++)
			for (int j = 0; j < length; j++)
				ret.locate(i, j) = me.locate(j, i);
		return ret;
	}

	///<seealso cref="https://www.cnblogs.com/xiaoxi666/p/6421228.html"></seealso>
	///<summary>求逆, 并返回</summary>
	_MyType LUP_solve_inverse() const
	{
		auto LUP_Solve = [](const _MyType& L, const _MyType& U, const vector<int, length> P, const _VecType& b)->_VecType
		{
			_VecType x, y;

			//正向替换
			for (int i = 0; i < length; i++)
			{
				y[i] = b[P[i]];
				for (int j = 0; j < i; j++)
				{
					y[i] = y[i] - L[i * length + j] * y[j];
				}
			}
			//反向替换
			for (int i = length - 1; i >= 0; i--)
			{
				x[i] = y[i];
				for (int j = length - 1; j > i; j--)
				{
					x[i] = x[i] - U[i * length + j] * x[j];
				}
				x[i] /= U[i * length + i];
			}
			return x;
		};
		auto LUP_Descomposition = [](_MyType& A, _MyType& L, _MyType& U, vector<int, length>& P)->void
		{
			int row = 0;
			P.init([](int i) { return i; });
			for (int i = 0; i < length - 1; i++)
			{
				double p = 0.0;
				for (int j = i; j < length; j++)
				{
					if (std::fabs(A[j * length + i]) > p)
					{
						p = std::fabs(A[j * length + i]);
						row = j;
					}
				}
				if (0 == p)
				{
					throw(matrix_exception::matrix_strange);
				}
				//交换P[i]和P[row]
				int tmp = P[i];
				P[i] = P[row];
				P[row] = tmp;
				double tmp2 = 0.0;
				for (int j = 0; j < length; j++)
				{
					//交换A[i][j]和 A[row][j]
					tmp2 = A[i * length + j];
					A[i * length + j] = A[row * length + j];
					A[row * length + j] = tmp2;
				}

				//以下同LU分解
				double u = A[i * length + i], l = 0.0;
				for (int j = i + 1; j < length; j++)
				{
					l = A[j * length + i] / u;
					A[j * length + i] = l;
					for (int k = i + 1; k < length; k++)
					{
						A[j * length + k] = A[j * length + k] - A[i * length + k] * l;
					}
				}
			}

			//构造L和U
			for (int i = 0; i < length; i++)
			{
				for (int j = 0; j <= i; j++)
				{
					if (i != j)
					{
						L[i * length + j] = A[i * length + j];
					}
					else
					{
						L[i * length + j] = 1;
					}
				}
				for (int k = i; k < length; k++)
				{
					U[i * length + k] = A[i * length + k];
				}
			}

		};

		//创建矩阵A的副本，注意不能直接用A计算，因为LUP分解算法已将其改变
		_MyType me_copy;
		_MyType ret; //最终的逆矩阵（还需要转置）
		_VecType inv_A_each;//矩阵逆的各列
		//_MyType B;
		_VecType b;//b阵为B阵的列矩阵分量
		_MyType L;
		_MyType U;
		vector<int, length> P;
		me_copy = me;
		LUP_Descomposition(me_copy, L, U, P);
		for (int i = 0; i < length; i++)
		{
			//构造单位阵的每一列
			b.init([]([[maybe_unused]] int x) { return 0; });
			b[i] = 1;

			inv_A_each = LUP_Solve(L, U, P, b);
			std::memcpy(ret.head.get() + i * length, inv_A_each.head.get(), length * sizeof(_Data));//将各列拼接起来

		}
		ret.dotrans();//由于现在根据每列b算出的x按行存储，因此需转置

		return ret;
	}

	_MyType L() const
	{
		_MyType ret = *this;
		for (int i = 0; i < this->length; i++)
		{
			for (int j = i; j < this->length; j++)
			{
				ret.locate(i, j) = 0;
			}
		}
		return ret;
	}
	_MyType D() const
	{
		_MyType ret;
		for (int i = 1; i < this->length; i++)
		{

			ret.locate(i, i) = this->locate(i, i);

		}
		return ret;
	}
	_MyType R() const
	{
		_MyType ret = *this;
		for (int i = 0; i < this->length; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				ret.locate(i, j) = 0;
			}
		}
		return ret;
	}
};


#pragma pack(pop)
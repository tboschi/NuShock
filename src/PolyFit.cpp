#include "tools/PolyFit.h"

PolyFit::PolyFit(std::vector<double> ax, std::vector<double> &dt, int n) :
	_order(n),
	_size(ax.size()),
	axis(std::vector<double>(ax.size())),
	data(std::vector<double>(dt.size())),
	beta(Eigen::VectorXd(_order+1))
{
	_weighted = (_size == dt.size() / 2);

	std::cout << "Loaded " << _size << " points" << std::endl;
	if (_weighted)
		std::cout << "Fitting with error weighting" << std::endl;
	else
		std::cout << "Fitting with fixed weigh" << std::endl;

	//sorting Y in X
	std::vector<int> index(_size);
	for (size_t i = 0; i < _size; ++i)
		index[i] = i;

	std::sort(index.begin(), index.end(), [&](int a, int b) { return ax[a] < ax[b]; });

	for (size_t i = 0; i < _size; ++i)
	{
		int j = index[i];
		axis[i] = ax[j];
		data[i] = dt[j];
		if (_weighted)
			data[i+_size] = dt[j+_size];
	}
}

std::vector<double> PolyFit::GetAxis()
{
	return axis;
}

std::vector<double> PolyFit::GetData()
{
	return data;
}

double PolyFit::ff(double x)
{
	double ret = beta(_order);
	for (size_t n = 0; n < _order; ++n)
		ret = x*ret + beta(_order-1 - n);

	return ret;
}

double PolyFit::fe(double x)
{
	Eigen::MatrixXd v = VarBeta();
	Eigen::VectorXd jac(_order+1);

	for (size_t n = 0; n < _order+1; ++n)
		jac(n) = pow(x, n);

	return (jac.transpose() * v * jac).norm();
}

double PolyFit::fd(double x)
{
	double ret = _order * beta(_order);
	for (size_t n = 1; n < _order; ++n)
		ret = x*ret + (_order-n) * beta(_order-1 - n);

	return ret;
}

Eigen::VectorXd PolyFit::Value()	//vector of residual
{
	Eigen::VectorXd val(_size);
	for (size_t i = 0; i < _size; ++i)
		val(i) = data[i];

	return val;
}

Eigen::VectorXd PolyFit::Residual()	//vector of residual
{
	Eigen::VectorXd res(_size);
	for (size_t i = 0, e = _size; i < _size; ++i, ++e)
		res(i) = (data[i] - ff(axis[i]));

	return res;
}

Eigen::MatrixXd PolyFit::Jacobian()	//Jacobian matrix
{
	Eigen::MatrixXd jac(_size, _order+1);
	for (size_t i = 0, e = _size; i < _size; ++i, ++e)
		for (size_t n = 0; n < _order+1; ++n)
			jac(i, n) = pow(axis[i], n);

	return jac;
}

Eigen::MatrixXd PolyFit::Weight()	//Jacobian matrix
{
	if (_weighted)
	{
		Eigen::MatrixXd weight = Eigen::MatrixXd::Zero(_size, _size);
		for (size_t i = 0, e = _size; i < _size; ++i, ++e)
			weight(i, i) =  1.0 / data[e];

		return weight;
	}
	else
		return Eigen::MatrixXd::Identity(_size, _size);
}

//Least square solver
Eigen::VectorXd PolyFit::Delta()	//this returns increment to new beta params
{
	Eigen::MatrixXd j = Jacobian();
	Eigen::MatrixXd r = Residual();
	Eigen::MatrixXd w = Weight();

	return (j.transpose() * w * j).ldlt().solve(j.transpose() * w * r);	//using Norm eq.
}

Eigen::VectorXd PolyFit::Beta0()	//this returns increment to new beta params
{
	Eigen::MatrixXd j = Jacobian();
	Eigen::MatrixXd y = Value();
	Eigen::MatrixXd w = Weight();

	return (j.transpose() * w * j).ldlt().solve(j.transpose() * w * y);	//using Norm eq.
}

Eigen::MatrixXd PolyFit::VarBeta()	//this returns increment to new beta params
{
	Eigen::MatrixXd j = Jacobian();
	Eigen::MatrixXd r = Residual();
	Eigen::MatrixXd w = Weight();

	double res = (r.transpose() * w * r).norm() / (_size - _order - 1);
	return res * (j.transpose() * w * j).inverse();
}

std::vector<double> PolyFit::Solve(bool verbose)
{
	std::vector<double> res = LeastSquare(verbose);

	if (verbose)
		Print(res);

	return res;
}

std::vector<double> PolyFit::LeastSquare(bool verbose)
{
	//find starting point for a, b, c, and d
	int c = 0;
	double err = 1e-9;
	Eigen::VectorXd delta = Eigen::VectorXd::Constant(_order+1, 1);

	beta = Beta0();

	while (delta.norm() > err && c < 100)
	{
		++c;
		delta = Delta();
		beta += delta;

		if (verbose)
		{
			std::cout << "Iteration " << c << std::endl;
			std::cout << "\tdelta :  " << delta.transpose() << std::endl;
			std::cout << "\tbeta  :  " << beta.transpose() << std::endl;
			std::cout << "\terror :  " << delta.norm() << std::endl;
		}
	}

	Eigen::MatrixXd var;
	if (c >= 100)
		std::cout << "can't converge" << std::endl;
	else
		var = VarBeta();

	std::vector<double> vb(beta.data(), beta.data() + beta.size());
	vb.insert(vb.end(), var.data(), var.data() + var.size());

	return vb;
}

void PolyFit::Print(std::vector<double> &r)
{
	//size of r is order+1 + (order+1)²
	//p0, p1, p2... e00, e01, e02...e10, e11, e12...e20, e21, e22...
	int s = _order+1;

	std::cout << "Fit result" << std::endl;
	for (int i = 0; i < s; ++i)
		std::cout << "p[" << i << "]\t" << r[i] << "\t+/- "
			  << sqrt(r[i + s * (1 + i)]) << "\t("
			  << 100.*sqrt(r[i + s * (1 + i)])/std::abs(r[i]) << ")\n";

	//size
	//e00	e01	e02
	//e10	e11	e12
	//e20	e21	e22
	std::cout << "Correlation matrix" << std::endl;
	for (int i = 0; i < s; ++i)
	{
		for (int j = 0; j < i+1; ++j)
			std::cout << r[s + i + s * j] / sqrt(r[s + i + s*i] * r[s + j + s*j]) << "\t";
		std::cout << std::endl;
	}
}

void PolyFit::Reset()
{
	beta = Eigen::VectorXd::Constant(_order+1, 1.0);
}

#pragma once
#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>

#define PAIR_PAIR std::pair<std::pair<double, double>, double>
#define VEC_PAIR std::vector<PAIR_PAIR>

class Chebyshev
{
public:
  Chebyshev(int n, int m, int K) : N(n + 1), M(m + 1)
  {
    h = (b - a) / n;
    k = (d - c) / m;

    double lambdaMin = 4.0 / (h * h) * pow(sin(M_PI / (2.0 * n)), 2) + 4.0 / (k * k) * pow(sin(M_PI / (2.0 * m)), 2);
    double lambdaMax = 4.0 / (h * h) * pow(sin(M_PI * (n - 1) / (2.0 * n)), 2) + 4.0 / (k * k) * pow(sin(M_PI * (m - 1) / (2.0 * m)), 2);

    //tau = std::vector<double>(K);
    for (size_t i = 0; i < K; ++i)
    {
      tau.push_back(2.0 / (lambdaMin + lambdaMax + (lambdaMax - lambdaMin) * cos(M_PI * (1.0 + 2.0 * i) / (2 * K))));
    }
  }

  virtual ~Chebyshev() = default;

  void discrepancy()
  {
    double coeffA = -2.0 / (h * h);

    for (size_t i = 1; i < N - 1; ++i)
    {
      for (size_t j = 1; j < M - 1; ++j)
      {
        r[i][j].second = f[i][j].second + (numericalSol[i - 1][j].second - 2.0 * numericalSol[i][j].second + numericalSol[i + 1][j].second) / (h * h)
                      + (numericalSol[i][j - 1].second - 2.0 * numericalSol[i][j].second + numericalSol[i][j + 1].second) / (k * k);
      }
    }
  }


  double inf_discrepancy()
  {
    double tmp = 0;
    for (size_t i = 0; i < N; ++i)
    {
      for (size_t j = 0; j < M; ++j)
      {
        if (tmp < abs(r[i][j].second))
          tmp = r[i][j].second;
      }
    }
    return tmp;
  }

  void run()
  {
    for (size_t i = 0; i < 10000; ++i)
    {
      discrepancy();
      for (size_t k = 0; k < N - 1; ++k)
      {
        for (size_t j = 0; j < M - 1; ++j)
        {
          numericalSol[k][j].second += tau[i % tau.size()] * r[k][j].second;
        }
      }
    }
    //printf("\n\t |||||||||||||||||||||||||||||| \n\t\t Numerical Solution \n");
    //printV(numericalSol);
  }

protected:
  double virtual foo(const std::pair<double, double>&) = 0;
  double virtual foo(const double&, const double&) = 0;
  double virtual mu1(const double&) = 0;
  double virtual mu2(const double&) = 0;
  double virtual mu3(const double&) = 0;
  double virtual mu4(const double&) = 0;

  void printV(const std::vector<VEC_PAIR>& matr) const
  {
    for (int i = M - 1; i >= 0; --i)
    {
      for (int j = 0; j < N; ++j)
      {
        std::cout << matr[j][i].second << " ";
      }
      printf("\n");
    }
    printf("\n");
    //for (int i = M - 1; i >= 0; --i)
    //{
    //  for (int j = 0; j < N; ++j)
    //  {
    //    std::cout << '(' << matr[j][i].first.first << ", " << matr[j][i].first.second << ") ";
    //  }
    //  printf("\n");
    //}
    printf("\n");
  }
  void initCondition()
  {
    double posX = a;
    double posY = c;
    for (size_t i = 0; i < N; ++i, posX += h)
    {
      r.push_back(VEC_PAIR{});
      f.push_back(VEC_PAIR{});
      numericalSol.push_back(VEC_PAIR{});
      for (size_t j = 0; j < M; ++j, posY += k)
      {
        numericalSol[i].push_back(PAIR_PAIR{ std::pair<double, double>{posX, posY}, 0.0 });
        f[i].push_back(PAIR_PAIR{ std::pair<double, double>{posX, posY}, 0.0 });
        f[i][j].second = foo(f[i][j].first);
        r[i].push_back(PAIR_PAIR{ std::pair<double, double>{posX, posY}, 0.0 });
      }
      posY = c;
    }

    for (size_t i = 0; i < N; ++i)
    {
      numericalSol[i][0].second = mu3(numericalSol[i][0].first.first);
      numericalSol[i][M - 1].second = mu4(numericalSol[i][M - 1].first.first);
    }
    for (size_t i = 0; i < M; ++i)
    {
      numericalSol[0][i].second = mu1(numericalSol[0][i].first.second);
      numericalSol[N - 1][i].second = mu2(numericalSol[N - 1][i].first.second);
    }
  }

protected:
  double a = -1;
  double b = 1;
  double h; // шаг по x \in [a, b] = (b-a)/n, N - разбиение сетки
  double c = -1;
  double d = 1;
  double k; // шаг по y \in [c, d] = (d-c)/m, M - разбиение сетки
  int N;
  int M;
  std::vector<double> tau;
  std::vector<VEC_PAIR> numericalSol; // ~V
  std::vector<VEC_PAIR> f; // 
  std::vector<VEC_PAIR> r;
};

class Test_function : public Chebyshev
{
public:
  Test_function(int n, int m, int K) : Chebyshev(n, m, K)
  {
    initCondition();
    trueSol = numericalSol;
    for (size_t i = 0; i < M; ++i)
    {
      for (size_t j = 0; j < N; ++j)
      {
        trueSol[j][i].second = true_sol(numericalSol[j][i].first);
      }
    }

    //discrepancy();
    //printV(numericalSol);
    //printf("\n\t |||||||||||||||||||||||||||||| \n\t\t True Solution \n");
    //printV(trueSol);
    //printf("\n\t |||||||||||||||||||||||||||||| \n\t\t Discrepancy (r) \n");
    //printV(r);
  }

  double true_sol(const double& x, const double& y)
  {
    return exp(1.0 - x * x - y * y);
  }
  double true_sol(const std::pair<double, double>& xy)
  {
    return exp(1.0 - xy.first * xy.first - xy.second * xy.second);
  }

protected:
  double foo(const double& x, const double& y) override
  {
    return 4.0 * exp(1.0 - x * x - y * y) * (1.0 - x * x - y * y);
  }
  double foo(const std::pair<double, double>& xy) override
  {
    return 4.0 * exp(1.0 - xy.first * xy.first - xy.second * xy.second) * (1.0 - xy.first * xy.first - xy.second * xy.second);
  }
  double mu1(const double& y) override
  {
    return exp(-y * y);
  }
  double mu2(const double& y) override
  {
    return exp(-y * y);
  }
  double mu3(const double& x) override
  {
    return exp(-x * x);
  }
  double mu4(const double& x) override
  {
    return exp(-x * x);
  }

protected:
  std::vector<VEC_PAIR> trueSol;
};

#pragma once

#include <complex>
#include <vector>

using Complex = std::complex<double>;

class ExponentialTable {
public:
  ExponentialTable() {}

  Complex exp(Complex z, int n) {
    if (z != base) {
      buildTable(z, 0, n);
      base = z;
    } else if (n >= table.size()) {
      buildTable(z, table.size(), n);
    }
    return table[n];
  }

private:
  Complex base = 0.0;
  std::vector<Complex> table;

  void buildTable(Complex base, int start, int end) {
    if (start == 0) {
      table.assign(1, 1.0);
      start = 1;
    }
    table.resize(end + 1);
    for (int i = start; i <= end; ++i) {
      table[i] = table[i - 1] * base;
    }
  }
};

class InverseExponentialTable {
public:
  InverseExponentialTable() {}

  Complex inv_exp(Complex z, int n) {
    if (z != base) {
      buildTable(z, 0, n);
      base = z;
    } else if (n >= table.size()) {
      buildTable(z, table.size(), n);
    }
    return table[n];
  }

private:
  Complex base = 0.0;
  std::vector<Complex> table;

  void buildTable(Complex base, int start, int end) {
    if (start == 0) {
      table.assign(1, 1.0);
      start = 1;
    }
    table.resize(end + 1);
    for (int i = start; i <= end; ++i) {
      table[i] = table[i - 1] / base;
    }
  }
};

class BinomialTable {
public:
  BinomialTable() {}

  long long binomial(int n, int k) {
    if (k < 0 || k > n) {
      return 0;
    }
    if (n != this->n) {
      buildTable(n);
    }

    return table[k];
  }

private:
  int n = 0;
  std::vector<long long> table = {1};

  void buildTable(int n) {
    table.assign(n + 1, 0);
    table[0] = 1;
    for (int i = 1; i <= n; i++) {
      for (int j = i; j >= 1; j--) {
        table[j] += table[j - 1];
      }
    }
    this->n = n;
  }
};

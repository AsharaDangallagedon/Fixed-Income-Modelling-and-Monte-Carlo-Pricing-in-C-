#### Fixed Income Modelling and Monte Carlo Pricing in C++

This repository contains a C++17 implementation of an end-to-end fixed income workflow combining:

- Discount curve bootstrapping from market swap rates
- Vasicek short-rate term-structure analytics
- Exact short-rate simulation
- Monte Carlo pricing of bonds and interest-rate options
- Parameter sensitivity analysis

The project links deterministic curve construction with stochastic short-rate modelling in a single, reproducible workflow.



## Project Overview

The objective was to build a clear and internally consistent fixed income framework starting from observable swap market inputs and moving into model-based pricing under a short-rate process.

The workflow is:

1. Bootstrap a discount curve from quoted par swap rates  
2. Recover zero rates and implied swap pricing quantities  
3. Implement analytical Vasicek bond pricing formulas  
4. Simulate short-rate paths using exact Gaussian transitions  
5. Price bonds and option-style payoffs via Monte Carlo  
6. Analyse sensitivity to mean reversion and volatility  

The focus is clarity, validation, and clean implementation rather than model complexity.



## 1. Discount Curve Construction

The discount curve is bootstrapped from semi-annual par swap rates using:

- Sequential solving of discount factors
- Linear interpolation on discount factors
- Semi-annual fixed-leg cashflows

Market inputs span maturities from 1 to 10 years.

The curve can be queried for:

- Zero-coupon bond prices
- Continuously-compounded zero rates
- Swap valuation
- Implied par rates at non-quoted maturities

The resulting discount factor curve is monotone decreasing and produces an upward-sloping zero rate curve under the chosen inputs.



## 2. Vasicek Term Structure Model

The short rate follows the Vasicek process:

dr_t = a(b − r_t)dt + σ dW_t

The model provides:

- Closed-form zero-coupon bond pricing
- Analytical spot rate curve
- Analytical forward rate curve

These analytical outputs are used to validate the Monte Carlo engine.



## 3. Monte Carlo Engine

Short-rate paths are generated using the exact conditional Gaussian transition of the Vasicek model rather than an Euler discretisation scheme.

Pathwise discount factors are computed using trapezoidal integration.

The engine prices:

- 5-year zero-coupon bond
- Short-rate call option
- Bond call option

Monte Carlo bond prices are compared against analytical benchmarks to verify correctness.



## 4. Sensitivity Analysis

A small parameter grid explores:

- Mean reversion speed (a)
- Volatility (σ)

Observed behaviour:

- Higher volatility increases option values
- Faster mean reversion reduces option prices
- Stronger mean reversion lowers long-maturity bond prices in this parameter setup

Results are consistent with financial intuition.



## Code Structure

The implementation is organised into small, modular classes:

- DiscountCurve — stores discount factors and performs interpolation
- Swap — generates cashflows and solves for par rates
- Bootstrapper — constructs the discount curve
- Vasicek — analytical bond pricing and simulation logic
- ProjectConfig — centralised configuration of market inputs, model parameters, and Monte Carlo settings

Language standard: C++17



## Limitations

- Linear interpolation on discount factors (simple but basic)
- Vasicek model allows negative rates
- No calibration to live market data
- Option estimates require larger simulation budgets for tighter convergence

This project is designed as a clear and inspectable fixed income modelling framework rather than a production-grade rates library.

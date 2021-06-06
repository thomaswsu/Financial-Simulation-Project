# Financial-Simulation-Project

This project implements a delta hedging strategy from a self-funding portfolio. As the delta of our portfolio changes, the portfolio is rebalanced such that the delta of the new portfolio is zero. Common random numbers (CRN) and the finite difference method is used to approximate delta. The repository pulls data from QuantMod to gather historical data and generates a market simulation using Heston's Stochastic Volatility Model with parameters generated from the gathered data. Our code works with European and Barrier options.

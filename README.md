# **Coverage Probability of Confidence Intervals for Correlation Coefficients**

## **Overview**
This project performs a **Monte Carlo simulation** to evaluate the **coverage probability** of confidence intervals for the **sample correlation coefficient** (\( r_n \)). The study involves:

- Generating random samples from a **bivariate normal distribution** with different true correlation coefficients (\( \rho \)).
- Computing the **sample correlation coefficient** (\( r_n \)) for each sample.
- Applying **Fisher’s transformation** (\( g(r_n) = \text{atanh}(r_n) \)) to stabilize variance.
- Constructing **confidence intervals** for \( \rho \) using Fisher’s transformed values.
- **Estimating coverage probability** by checking how often the true \( \rho \) falls within the interval.

---

## **Files**
- `correlation_simulation.R`: The main R script that runs the simulation and generates confidence intervals.
- `README.md`: This file, explaining the project structure and purpose.

---

## **Methodology**
### **1. Simulation Setup**
- The study considers sample sizes: **\( n \in \{15, 25\} \)**.
- The true correlation coefficient (\( \rho \)) takes values: **\( \{0.0, 0.2, 0.4, 0.6, 0.8\} \)**.
- Each scenario is simulated **1,000 times**.

### **2. Statistical Procedures**
- Generate bivariate normal data using `MASS::mvrnorm()`.
- Compute **sample correlation \( r_n \)**.
- Apply **Fisher’s transformation** for variance stabilization:
- Compute **coverage probability**, i.e., the proportion of times the interval contains the true \( \rho \).

---

## **Usage**
### **1. Clone the Repository**
```bash
git clone https://github.com/yourusername/correlation-coverage.git
cd correlation-coverage

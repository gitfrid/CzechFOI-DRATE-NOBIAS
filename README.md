# CzechFOI-DRATE-NOBIAS

**CzechFOI-DRATE-NOBIAS** is a data analysis project focused on identifying and correcting **non-random boundary condition bias** in observational vaccine and mortality data from the Czech Republic.

## Solving the Non-Random Boundary Condition Bias

The **non-random boundary condition bias** occurs when a homogeneous group is **artificially split into subgroups based on a non-random time point**, such as the time of vaccination with the **constraint death day > last dose day**. These subgroups are then compared in terms of outcomes (e.g., death rates), leading to **illusory differences** that **do not reflect any true causal effect**.

This bias is particularly dangerous because it can produce **misleading conclusions** even when no real treatment effect exists.

### Methode
A population with a 
constant mortality rate is simulated, and the real vaccine doses 1–7 are randomly assigned to individuals, repeating the process until the death day of the selected individual is after the last vaccination. 
This creates a dataset with a constant homogeneous mortality rate and a realistic vaccination schedule.
<br><br>The simulated (HR=1) and real dataset is then used with scientific methods such as the time-varying Cox model to prove whether the methode eliminate this bias in the data, or produce distorted results.

### Scientific Foundations

A few serious scientists, genuinely interested in knowledge, have addressed this bias rigorously, including **Miguel Hernán** and **James Robins**. 
<br>Their textbook:

> **"Causal Inference: What If"**  
> by Miguel Hernán & James Robins  
> 📘 https://miguelhernan.org/whatifbook

thankfully is available **open-access** and provides theoretical and practical guidance on this topic — especially in **Chapters 19 to 21**, where they explain how to correct such biases using **causal models and time-aligned person-day analysis**.

This project draws heavily on those principles, applying them to large-scale real-world data to avoid common pitfalls like **immortal time bias**, **lag bias**, and **non-random group splitting**.

[**See also CzechFOI-DRATE_EXAM project for Investigation of the Bias**](https://github.com/gitfrid/CzechFOI-DRATE_EXAM/tree/main)

_________________________________________
## Comparison vaccinated vs. unvaccinated with eliminated bias - using different Methodes  

I have not yet found a scientific Cox/Poisson method that completely eliminates the bias **caused by the restriction “date of death > date of last dose”**, which is inevitably present in the real observational data, in order to fairly compare vaccinated and unvaccinated individuals and calculate efficacy.

**Assuming that the bias in both data sets is comparable and both produce similar biased results, it follows that the true unbiased effect of vaccinated vs unvaccinated individuals is negligible.**

<br>Methods such as the G-formula or target trial emulation with clone censoring weighting (CCW) should be able to correct this bias, but require in-depth knowledge to apply.
Most scientific studies on this topic do not publish their data or the code used for evaluation, and probably use standard methods that do not correct for this bias.
<br>Since most scientists use R code, I will shortly be publishing R-code scripts.
_________________________________________

## The Solution of Non-Random Boundary Condition Bias

This is what **Hernán & Robins** say:

> **“Don’t compare vaccinated people with unvaccinated people. Compare person-time under vaccination vs. person-time without vaccination.”**

---

### Why is this not a problem for comparison?

Because you no longer say:

> **“Group A (vx) vs. Group B (uvx)”**

But rather:

> **“What is the hazard rate during vaccinated periods compared to unvaccinated periods”**

This is **fair** because:

- All individuals contribute to the risk time of **both groups** (depending on their vaccination status over time).
- No one is excluded **"from the outset"** due to the timing of vaccination.

---

### The **wrong** model (which is often used):

You put all vaccinated people in one group and all unvaccinated people in another – and then compare the two groups.

But:

> ❌ **That’s unfair!**

Because:

- The **bias is already baked into the observational real-world data** due to non-random splitting into the two groups.

That means:

- It doesn’t help to randomly assign doses only to the vaccinated group, as the bias is already introduced by real world constraint death day > last dose day.
- If you randomly assigne the doses to the entire population, it would eliminate the bias — but this would destroy the real relation of vx/uvx grouping you actually want to compare.
- Therefore, you need to find another way to eliminate the bias from the raw data — and **Hernán & Robins** showed how to do this.
- 
_________________________________________

### FW) cox time-varying Methode with Kaplan–Meier (KM) survival curve plot
Phyton script [FW) cox time-varying.py](https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Py%20Scripts/FW%29%20cox%20time-varying.py) [FW) cox time-varying.R](https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/R%20Scripts/FG%29%20simulate%20deaths%20doses%20curves.R)**-> hasn't been checked for errors now**
<br>
<br>**Simulated data (expected HR~1 / no effect-placebo)** [Results TXT](https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FW%29%20cox%20time-varying/FW-FG%29%20case3_sim_deaths_sim_real_doses_with_constraint%20AG70%20cox%20time-varying.TXT)
<br>There is still a slight distortion, as vx and uvx must theoretically overlap. 
<br>
<img src=https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FW)%20cox%20time-varying/FW-FG)%20case3_sim_deaths_sim_real_doses_with_constraint%20AG70%20cox%20time-varying.png width="1280" height="auto">
<br>

<br>**Vs. real Czech data** [Results TXT](https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FW%29%20cox%20time-varying/FW%29%20Vesely_106_202403141131_AG70%20cox%20time-varying.TXT)
<br>Real data showing only a very small effect, which corresponds approximately to the magnitude of the residual distortion of the simulated data.
<br>
<img src=https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FW%29%20cox%20time-varying/FW%29%20Vesely_106_202403141131_AG70%20cox%20time-varying.png width="1280" height="auto">
<br>
_________________________________________

### FW) CoX Time-Varying Survival Analysis Results:

- **Study Period sart date 2020-01-01 until END_MEASURE 1095 days**
- **Follow-up window for life years saved calculation:** up to day 734

**Notes**
- The model uses **time-varying covariates** to avoid immortal time bias.
- `vaccinated_time` allows modeling of **waning protection** over time.
---

### Simulated data (expected HR~1 / no effect-placebo): 

Model Convergence Converged after 5 iterations

### Model Coefficients and Hazard Ratios

| Covariate        | Coef     | exp(Coef) (HR) | 95% CI (HR)        | p-value        | Interpretation                         |
|------------------|----------|----------------|---------------------|----------------|----------------------------------------|
| `vaccinated`     | -0.147   | **0.863**      | 0.843 – 0.884       | 3.1e-34        | **13.7% lower death hazard**           |
| `t`              | -0.00097 | 0.999          | 0.999 – 0.999       | ≈ 0            | Slightly decreasing baseline hazard    |
| `vaccinated_time`| -0.00048 | 0.9995         | 0.9995 – 0.9996     | 7.3e-130       | Very mild waning over time             |


### Key Metrics

- **Vaccine Effectiveness (VE):** `1 - HR = 1 - 0.863 = 13.7%` **still small bias present as expected HR is ~1**
- **Statistical significance:** all covariates are **highly significant** (p < 10⁻³⁰)

--

### vs. Real Czech data: 

Model Convergence Converged after 5 iterations

### Model Coefficients and Hazard Ratios

| Covariate        | Coef     | exp(Coef) (HR) | 95% CI (HR)        | p-value         | Interpretation                        |
|------------------|----------|----------------|---------------------|------------------|---------------------------------------|
| `vaccinated`     | -0.077   | **0.926**      | 0.904 – 0.948       | 1.56e-10         | **7.4% lower death hazard**           |
| `t`              | -0.00094 | 0.999          | 0.999 – 0.999       | ≈ 0              | Slight baseline decline over time     |
| `vaccinated_time`| -0.00046 | 0.9995         | 0.9995 – 0.9996     | 1.44e-118        | Very mild waning effect               |


### Key Metrics

- **Vaccine Effectiveness (VE):** `1 - HR = 1 - 0.926 = 7.4%` **probably because still small bias present as the simulated data show**
- **All covariates statistically significant** at extremely low p-values (p < 1e-10)

_________________________________________
### FY) cox time-varying Methode per Dose with Kaplan–Meier (KM) survival curve plot
Phyton script [FY) cox time-varying per dose.py](https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Py%20Scripts/FY%29%20cox%20time-varying%20per%20dose.py)
<br>
<br>**Simulated data (expected HR~1 / no effect-placebo)** [Results TXT](https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FY%29%20cox%20time-varying%20per%20Dose/FY-FG%29%20case3_sim_deaths_sim_real_doses_with_constraint%20AG70%20cox%20time-varying%20per%20dose.TXT)
<br>There is still a slight distortion, as vx and uvx should theoretically overlap horizontally. 
<br>
<img src=https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FY)%20cox%20time-varying%20per%20Dose/FY-FG)%20case3_sim_deaths_sim_real_doses_with_constraint%20AG70%20cox%20time-varying%20per%20dose.png width="1280" height="auto">
<br>

<br>**Vs. real Czech data per Dose** [Results TXT](https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FY%29%20cox%20time-varying%20per%20Dose/FY%29%20real%20data%20Vesely_106_202403141131_AG70%20cox%20time-varying%20per%20dose.TXT)
<br>Real data showing only a small effect, approximately a little bigger then the magnitude of the residual distortion of the simulated data.
<br>
<img src=https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FY%29%20cox%20time-varying%20per%20Dose/FY%29%20real%20data%20Vesely_106_202403141131_AG70%20cox%20time-varying%20per%20dose.png width="1280" height="auto">
<br>
_________________________________________
### FZ) Poisson Methode with Kaplan–Meier (KM) survival curve plot
Phyton script [FZ) poisson.py](https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Py%20Scripts/FZ%29%20poisson.py)
<br>
<br>**Simulated data (expected HR~1 / no effect-placebo)** [Results TXT](https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FZ%29%20poisson/FZ-FG%29%20case3_sim_deaths_sim_real_doses_with_constraint%20AG70%20poisson.TXT)
<br>There is still a slight distortion, as vx and uvx must theoretically overlap horizontally. 
<br>
<img src=https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FZ%29%20poisson/FZ-FG%29%20case3_sim_deaths_sim_real_doses_with_constraint%20AG70%20poisson_KM_survival.png width="1280" height="auto">
<br>

<br>**Vs. real Czech data** [Results TXT](https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FZ%29%20poisson/FZ%29%20real%20data%20Vesely_106_202403141131_AG70%20poisson.TXT)
<br>Real data showing only a small effect, which corresponds approximately to the magnitude of the residual distortion of the simulated data.
<br>
<img src=https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FZ%29%20poisson/FZ%29%20real%20data%20Vesely_106_202403141131_AG70%20poisson_KM_survival.png width="1280" height="auto">
<br>
_________________________________________

## FZ) Poisson Regression

**Study Period sart date 2020-01-01 until END_MEASURE 1095 days**

**Notes**
- This model uses **person-day data**, expanding individual records by time.
- A **Poisson model with log link** is typically appropriate for **rare event count data** like deaths per day.
- The null result here may indicate **insufficient effect** in the simulated data (HR ~ 1), or structural bias due to **model or data design**.
---

### Simulated Deaths with random real dose schedule AG70 (expected HR~1 / no effect-placebo) Poisson Methode

**PerfectSeparationWarning**:  
  `Perfect separation or prediction detected, parameter may not be identified.`  
  This indicates that the model perfectly predicts the outcome for some observations — a sign of either overfitting or a non-informative predictor.

### Coefficients and IRRs

| Covariate   | Coef        | IRR (exp(coef)) | 95% CI        | p-value | Interpretation |
|-------------|-------------|------------------|---------------|---------|----------------|
| `const`     | ≈ 0         | **1.000**        | [1.000, 1.000]| 1.000   | No baseline death risk change |
| `vaccinated`| ≈ 0         | **1.000**        | [1.000, 1.000]| 1.000   | **No effect detected**        |
| `age_c`     | 0           | –                | [0, 0]        | NaN     | Possibly constant or missing  |

### Interpretation

- **No difference** in death incidence between vaccinated and unvaccinated individuals was detected in this simulation.
- The coefficient for `vaccinated` is extremely close to zero, with an **Incidence Rate Ratio (IRR) of 1.000**, meaning **no relative change in death rate**.
- **Perfect separation** warning suggests either:
  - There is no variation in the outcome with respect to predictors.
  - The data or simulation design causes deterministic outcomes (e.g., deaths assigned in a fixed or rule-based way).
- `age_c` was dropped due to data only for AG70

--

### vs. Real Czech data AG70 Poisson Methode

### Coefficients and Incidence Rate Ratios (IRRs)

| Covariate     | Coef     | IRR (exp(coef)) | 95% CI IRR         | p-value | Interpretation |
|---------------|----------|------------------|---------------------|---------|----------------|
| `const`       | -0.0003  | ~1.000           | [0.999, 1.000]      | 0.009   | Slightly lower baseline risk |
| `vaccinated`  | -0.0013  | **0.999**        | [0.998, 0.999]      | <0.001  | **Statistically significant reduction** in death rate for vaccinated individuals |
| `age_c`       | 0        | –                | [0, 0]              | NaN     | Possibly dropped due to no variation or collinearity |


### Interpretation

- **Vaccinated individuals** show a **statistically significant** but **very small reduction** in daily death rates compared to unvaccinated.
  - **IRR = 0.999** → about **0.1% lower** death rate.
- The **pseudo R² is low (0.0335)**, indicating that vaccination explains only a small portion of the variability in death counts.
- The model converged after 5 iterations with valid diagnostics.
  
____________________________________________________

## Essentials to Eliminate the Bias

### Use person-time analysis instead of person-group comparison:
- Don’t divide individuals into **"vaccinated"** vs. **"unvaccinated"** groups.
- Instead, **split each individual's time** into **vaccinated** and **unvaccinated periods**.

### Ensure that follow-up time is correctly aligned:
- The **start of follow-up** must be defined equally for all individuals (e.g., a shared baseline date or event).
- Avoid **immortal time bias** — no one should be counted as "at risk" before they are actually at risk.
- For example use **Target Trial Emulation (TTE)** veryone in the emulated trial starts at the same time — e.g., the date when all were eligible for vaccination. This eliminates biases caused by unequal follow-up periods or survival prerequisites.

### Allow individuals to contribute to both exposure states:
- A person can contribute **unvaccinated person-time** before receiving a dose and **vaccinated person-time** afterward.
- This ensures comparisons are **within the same population**, not between fundamentally different groups.

### Don’t condition on future information:
- You must **not use future events** (e.g., dose received later) to define exposure status at earlier times.

_________________________________________
## Software Requirements:

The aggregation is handled directly by Python scripts, which can generate aggregated CSV files very quickly.
For coding questions or help, visit https://chatgpt.com.

- [Python 3.12.5](https://www.python.org/downloads/) to run the scripts.
- [Visual Studio Code 1.92.2](https://code.visualstudio.com/download) to edit and run scripts.


## Disclaimer:
**The results have not been checked for errors. Neither methodological nor technical checks or data cleansing have been performed.**

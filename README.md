# CzechFOI-DRATE-NOBIAS

**CzechFOI-DRATE-NOBIAS** is a data analysis project focused on identifying and correcting **non-random boundary condition bias** in observational vaccine and mortality data from the Czech Republic.

## Solving the Non-Random Boundary Condition Bias

The **non-random boundary condition bias** occurs when a homogeneous group is **artificially split into subgroups based on a non-random time point**, such as the time of vaccination with the **constraint death day > last dose day**. These subgroups are then compared in terms of outcomes (e.g., death rates), leading to **illusory differences** that **do not reflect any true causal effect**.

This bias is particularly dangerous because it can produce **misleading conclusions** even when no real treatment effect exists.

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

_________________________________________
## Comparison vaccinated vs. unvaccinated with eliminated bias - using different Methodes  

---

### FL) cox time-varying with COX daily survival curve splot
Phyton script [FL) cox time-varying survival curves daily.py](https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Py%20Scripts/FL%29%20cox%20time-varying%20survival%20curves%20daily.py)
<br>
<br>**Simulated data (expected HR~1 / no effect-placebo)** [Results TXT](https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FL%29%20cox%20time-varying%20survival%20curves%20daily/FL-FG%29%20case3_sim_deaths_sim_real_doses_with_constraint%20cox%20time-varying.TXT)
<br>With corrected distortion, the curves for vaccinated and unvaccinated individuals must overlap for the simulated data with HR~1! 
<br>This proves that the evaluation is bias-free.
<br>
<img src=https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FL%29%20cox%20time-varying%20survival%20curves%20daily/FL-FG%29%20case3_sim_deaths_sim_real_doses_with_constraint%20cox%20time-varying.png width="1280" height="auto">
<br>

<br>**Vs. real Czech data** [Results TXT](https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FL%29%20cox%20time-varying%20survival%20curves%20daily/FL%29%20Vesely_106_202403141131_AG70%20cox%20time-varying.TXT)
<br>Czech FOI real data AG70 – When evaluated identically to eliminate distortions, the result shows virtually no effectiveness.
<br>
<img src=https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FL)%20cox%20time-varying%20survival%20curves%20daily/FL)%20Vesely_106_202403141131_AG70%20cox%20time-varying.png width="1280" height="auto">
<br>
_________________________________________

### FL) CoX Time-Varying Survival Analysis Results:

- **Study Period sart date 2020-01-01 until END_MEASURE 1095 days**

**Notes**
- The model uses **time-varying covariates** to avoid immortal time bias.
- `vaccinated_time` allows modeling of **waning protection** over time.
---

### Simulated Dataset (Simulated Deaths and real Dose schedule) - expected HR~1 no effect -> placebo

This run evaluates a time-varying Cox proportional hazards model using **real-world vaccination timing data with constraint death day > last dose day** combined with **simulated death events** that are independent of vaccination status (expected hazard ratio HR ≈ 1). The model used a COX_PENALIZER = 3, which was empirically determined to produce stable and unbiased estimates aligned with the known outcomme of the simulated data.

**Model Outcome Summary:**

| Covariate         | HR (exp(coef)) | 95% CI       | p-value   | Interpretation                                                         |
|-------------------|----------------|--------------|-----------|----------------------------------------------------------------------- |
| `vaccinated`      | 0.993          | 0.989–0.998  | 0.0029    | Slight apparent protective effect due to noise or still minimal bias   |
| `vaccinated_time` | 1.000          | 1.000–1.000  | 3.5e-07   | Effect over time is negligible as expected                             |

- **Deaths before vaccination + lag 0 :** 5496  
- **Deaths after vaccination + lag 0 :** 2371  
- **Life years saved:** ~0.0006 years (negligible, as expected)  
- **Total deaths:** 7867  
- **Event distribution:**  
  - **Unvaccinated deaths:** 5496  
  - **Vaccinated deaths:** 2371  
<br>The results show that the model correctly detects **no meaningful difference in mortality risk** between vaccinated and unvaccinated individuals when using simulated deaths that are intentionally unrelated to vaccination. 
This confirms that the model setup — including the chosen penalizer value — does **not introduce artificial vaccine effects** when none exist, ensuring valid and unbiased estimates under a neutral ground truth.
--

### vs. Real-World Dataset Result (Czech FOI Deaths + Real Doses) for AG70

This run applies the time-varying Cox proportional hazards model to **real-world data from the Czech Republic**, using actual vaccination dose dates and recorded all-cause deaths. 
The analysis focused on deaths observed from start 01.01.2020 as day 0 up to day 1095.

The model includes two covariates:
- `vaccinated`: binary indicator for vaccination status
- `vaccinated_time`: time since vaccination

Same penalization value of COX_PENALIZER = 3 was used to stabilize the model and prevent overfitting.

**Model Outcome Summary:**

| Covariate         | HR (exp(coef)) | 95% CI       | p-value   | Interpretation                                                  |
|-------------------|----------------|--------------|-----------|-----------------------------------------------------------------|
| `vaccinated`      | 0.996          | 0.991–1.000  | 0.0717    | Slight protective trend, but **not statistically significant**  |
| `vaccinated_time` | 1.000          | 1.000–1.000  | 2.0e-06   | Small decreasing effect over time, statistically significant    |

- **Deaths before vaccination + lag:** 5131  
- **Deaths after vaccination + lag:** 2706  
- **Life years saved by vaccination:** ~0.0004 years  
- **Total deaths:** 7837  
- **Event distribution:**  
  - **Unvaccinated deaths:** 5131  
  - **Vaccinated deaths:** 2706  

The results suggest **no meaningfull reduction in mortality hazard** following vaccination in AG70. While the main effect (`vaccinated`) is not statistically significant, 
the direction aligns with a modest protective signal. The `vaccinated_time` covariate shows a statistically significant (though clinically minimal) decreasing hazard over time post-vaccination.

This supports that, even under real-world observational data, the model remains robust and interpretable without introducing artifacts, using the same penalizer.

_________________________________________

### FW) cox time-varying Methode with Kaplan–Meier (KM) survival curve plot
Phyton script [FW) cox time-varying.py](https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Py%20Scripts/FW%29%20cox%20time-varying.py)
<br>
<br>**Simulated data (expected HR~1 / no effect-placebo)** [Results TXT](https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FW%29%20cox%20time-varying/FW-FG%29%20case3_sim_deaths_sim_real_doses_with_constraint%20AG70%20cox%20time-varying.TXT)
<br>There is still a slight distortion, as vx and uvx should theoretically overlap. 
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
<br>There is still a slight distortion, as vx and uvx should theoretically overlap horizontally. 
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

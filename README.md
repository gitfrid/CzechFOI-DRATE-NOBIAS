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
## Comparison vaccinated vs. unvaccinated with eliminated bias - using different Methodes  

---
### FW) cox time-varying Methode 
Phyton script [FW) cox time-varying.py](https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Py%20Scripts/FW%29%20cox%20time-varying.py)
<br>
<br>**Simulated data (expected HR~1 / no effect-placebo)** [Results TXT](https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FW%29%20cox%20time-varying/FW-FG%29%20case3_sim_deaths_sim_real_doses_with_constraint%20AG70%20cox%20time-varying.TXT)
<br>There is still a slight distortion, as vx and uvx should theoretically overlap horizontally. 
<br>
<img src=https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FW)%20cox%20time-varying/FW-FG)%20case3_sim_deaths_sim_real_doses_with_constraint%20AG70%20cox%20time-varying.png width="1280" height="auto">
<br>

<br>**Vs. real Czech data** [Results TXT](https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FW%29%20cox%20time-varying/FW%29%20Vesely_106_202403141131_AG70%20cox%20time-varying.TXT)
<br>Real data showing only a very small effect, which corresponds approximately to the magnitude of the residual distortion of the simulated data.
<br>
<img src=https://github.com/gitfrid/CzechFOI-DRATE-NOBIAS/blob/main/Plot%20Results/FW%29%20cox%20time-varying/FW%29%20Vesely_106_202403141131_AG70%20cox%20time-varying.png width="1280" height="auto">
<br>
_________________________________________
### FY) cox time-varying Methode per Dose
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
### FZ) Poisson Methode
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

## The Solution of Non-Random Boundary Condition Bias

This is what **Hernán & Robins** say:

> **“Don’t compare vaccinated people with unvaccinated people. Compare person-time under vaccination vs. person-time without vaccination.”**

---

### Why is this not a problem for comparison?

Because you no longer say:

> **“Group A (vx) vs. Group B (uvx)”**

But rather:

> **“What is the hazard rate during vaccinated periods compared to unvaccinated periods?”**

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

- It doesn’t help to randomly assign doses only to the vaccinated group — the bias is already introduced by the way vaccinated and unvaccinated groups were selected.
- If you randomly assigned doses to the entire population, it would eliminate the bias — but this would destroy the real grouping you actually want to compare.
- Therefore, you need to find another way to eliminate the bias from the raw data — and **Hernán & Robins** showed how to do this.

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


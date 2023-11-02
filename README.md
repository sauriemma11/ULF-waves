# ULFW - Ultra Low Frequency Waves

[//]: # (Magnetospheric Plasma Wave-Particle Interactions)

---


<div align="center">

![ulfwave_img.png](docs%2Fulfwave_img.png)

The ULFW project focuses on the analysis of Ultra-Low-Frequency (ULF) waves in
space
using data from the Geostationary Operational Environmental Satellite (GOES)
magnetic field sensors. It seeks to uncover insights into the relationship
between observed plasma wave parameters and their impact on electron
diffusion within Earth's Radiation belt.

This repo was created by students in the Fall 2023 semester of Software
Engineering for Scientists ([SWE4S](https://github.com/swe4s)) course at
University of Colorado, Boulder with the intent of helping answer this question
in future studies.

The figure below shows the important calculations and outputs. 
1) The first panel is a timeseries of the parallel component of the magnetic field filtered for the frequency band of interest. This shows when increases in amplitude occur. 
2) The second panel is a spectrogram which has time on the x axis, frequency on the y and the color is the power (blue is low, red is high). The red bursts of high power in the lower frequencies are the plasma waves we are interested in. You can also see the bursts correspond with amplitude increase in the above panel. 
3) The third panel is the average power spectral density (PSD) in the frequency band over time. This calculation is done in 3 hour intervals, hence why it is a line plot with steps. You can see that when the amplitude and power in the above plots increase, the average PSD also increases. 
4) The final panel is tau, or the electron radial diffusion time scale. Again, this correlated nicely as times with higher power have shorter timescales (more power -> diffuse electrons quicker). 
![example_output.png](docs%2Fexample_output.png)

</div>

---

<!-- TOC -->
### TOC
* [ULFW - Ultra Low Frequency Waves](#ulfw---ultra-low-frequency-waves)
    * [Project Goals](#project-goals)
    * [Key Features](#key-features)
    * [Usage](#usage)
    * [Updates](#updates)
<!-- TOC -->


---

### Project Goals

* Calculate observed wave properties and how long it would take an electron to
  diffuse one Earth radius towards Earth due to wave-particle interactions (also
  known as <i>diffusion time</i>)
* Use these parameters to better inform radiation belt models that currently
  rely on proxy indices

### Key Features

![codebase_ulf.png](docs%2Fcodebase_ulf.png)

* Data Analysis : ULFW uses GOES L2 high-res magnetic field data to perform analysis of
  Ultra-Low-Frequency (ULF) waves in space.
* Diffusion Time Calculation : Calculate the time it takes for an electron to
  diffuse one Earth radii towards Earth due to wave-particle interaction.
* Data Visualization : Generate time series plots
* User-Friendly : This project was designed to be user-friendly and with best
code practices.


---

## Usage

- Clone this repo and ensure all [environment](https://github.com/sauriemma11/ULF-waves/blob/main/env.yml) requirements are fulfilled
- Run from main project directory, execute with `$ ./run.sh`. **Expected run time for full day: 1-2 minutes.**

#### Examples how to run:
Mutable user inputs are:
- `filename` (.nc data file -- required)
- `fband` (frequency band in Hz -- default is [0.001, 0.01])
- `timespan` (time per window in hours, must go evenly into 24 -- default is 1).

From the main repository directory, run:
```shell
python main.py --filename foo.nc
```

## Updates

<details>
<summary>Expand for version release updates</summary>

### V 1.0
First full draft before the code review. Split one file that runs everything into different modules, created initial unit and functional tests, a main file to call all the functions, and a run.sh file.

</details>

## LICENSE

MIT License

Copyright (c) 2023 sauriemma11

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

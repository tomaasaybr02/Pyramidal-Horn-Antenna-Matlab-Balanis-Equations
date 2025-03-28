# Pyramidal Horn Antenna Project

## Overview
This repository contains all the code, figures, and documentation associated with the design, simulation, and experimental analysis of **Ka-band and C-band pyramidal horn antennas**. This project was developed for the *Advanced Antenna Techniques and Measurements* course at DTU, and explores electromagnetic theory, numerical simulations using HFSS, and real-world measurements conducted in anechoic chambers.

## Project Focus
The main objective is to compare theoretical, simulated, and experimental antenna performance across the **Ka-band (26.5–40 GHz)** and **C-band (3.95–5.85 GHz)**. Key parameters such as directivity, beamwidth, side lobe levels, and cross-polarization were analyzed.

## Repository Contents

- `balanis_equations.m` — MATLAB script implementing theoretical Balanis equations  
- `30430_TomasNiyaz_reportProject.pdf` — Final academic report with full analysis  
- `Figures/` — Folder containing figures, radiation patterns, and measurement setups  

> 🔍 *Note:* The `Figures/` folder includes all key radiation patterns (theoretical, simulated, experimental), antenna CAD views, and measurement environment photos for both frequency bands.

## Getting Started

To reproduce the theoretical results:
1. Open `balanis_equations.m` in MATLAB.
2. Run the script to compute the far-field patterns of the pyramidal horn using Balanis' method.
3. Results will include 2D radiation patterns, directivity plots, and comparative performance data.


## Figures and Visual Documentation

All major plots and measurement setup photos can be found in the `Figures/` folder, including:
- HFSS 3D designs
- Theoretical vs simulated vs experimental radiation patterns
- Anechoic chamber setups for Ka-band and C-band
- HPBW and FSLL extraction visuals

## License

This project is academic and open for educational purposes. Please cite the original report if you build upon this work.

## Authors

**Tomás Ruiz Aybar** – Theoretical modeling, MATLAB code  

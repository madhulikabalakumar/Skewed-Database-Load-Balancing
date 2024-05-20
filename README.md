# Skewed-Database-Load-Balancing

## Overview

This project explores and analyzes large flight and taxi datasets using advanced data processing and visualization techniques. The primary focus is on implementing effective bucket allocation strategies and creating heatmaps to visualize joint frequency distributions of attributes.

## Methodology

The project involves:
- Data loading and preprocessing
- Frequency-based bucket allocation 
- Heatmap generation for visualizing joint frequency distributions

## Results

The results include heatmaps and visualizations of the joint frequency distributions of various attribute pairs, providing insights into transportation trends and patterns.

## Getting Started

### Prerequisites
- Python 3.x
- Jupyter Notebook
- Libraries: `pandas`, `numpy`, `matplotlib`, `seaborn`

### Installation

Clone the repository:
```bash
git clone https://github.com/madhulikabalakumar/Skewed-Database-Load-Balancing
cd Skewed-Database-Load-Balancing
```

## Dataset Sources 

Flight Dataset: https://www.kaggle.com/datasets/polartech/flight-data-with-1-million-or-more-records
Taxi Dataset: https://www.kaggle.com/datasets/microize/newyork-yellow-taxi-trip-data-2020-2019?resource=download

## Usage 

Download the datasets first.  
Open and run the Jupyter notebooks:
```bash
jupyter notebook flight_data_analysis.ipynb
jupyter notebook taxi_data_analysis.ipynb
```
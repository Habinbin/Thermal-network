# Thermal-network
This files attribute the function and class for 1 dimensional thermal network analysis

## Installation
```
pip install git+https://github.com/Habinbin/Thermal-network.git
```

## Usage

Here's a basic example of how to use the library:

```python
import thermal_simulation as ts

# Create layers
layer1 = ts.SetLayer(L=0.1, dx=0.01, k=0.5, c=1000, rho=2000)
layer2 = ts.SetLayer(L=0.2, dx=0.02, k=1.0, c=800, rho=2500)

# Create construction
wall = ts.SetConstruction(layer1, layer2)

# Use functions
temp_celsius = 25
temp_kelvin = ts.D2K(temp_celsius)

# More examples...
```

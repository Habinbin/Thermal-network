# Thermal-network
This files attribute the function and class for 1 dimensional thermal network analysis

## Installation
```
pip install git+https://github.com/Habinbin/Thermal-network.git
```

## Usage

Here's a basic example of how to use the library:

```python
import thermal_simulation as tm

# Create layers
layer1 = tm.SetLayer(L=0.1, dx=0.01, k=0.5, c=1000, rho=2000)
layer2 = tm.SetLayer(L=0.2, dx=0.02, k=1.0, c=800, rho=2500)

# Create construction
wall = tm.SetConstruction(layer1, layer2)

# Use functions
temp_celsius = 25
temp_kelvin = tm.D2K(temp_celsius)

# More examples...
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.
```

5. LICENSE:
```
MIT License

Copyright (c) [year] [fullname]

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

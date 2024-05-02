# id15_NXstress

The `id15_NXstress` plugin is designed for converting raw data from the ID15 beamline into the NXstress format, suitable for use in the SOFT-AIS application. This tool analyzes raw data, and requires calibration files for both energy and angle to ensure accurate results.

## Installation

Before running the script, ensure you have Python installed along with the necessary libraries. It's recommended to use a virtual environment.

## Setup

To use the script, you need to create a `config.yml` file in the same directory as the script. This configuration file should specify paths to your data and calibration files, along with other necessary parameters.

### Configuration File Structure

Here is the required structure of the `config.yml` file:

```yaml
file_path: 'C:\path_to_file\Ni_Ring_0001.h5'
det_calib_file_angle: 'C:\path_to_file\angleCalib.h5'
det_calib_file_energy: 'C:\path_to_file\energyCalib.h5'
with_cradle: false
lattice: 'cubic'
phase_name: 'Nickel'
scanNbForRotation: 80
experimental_identifier: 'ihme19'
collection_identifier: 'write something here'
```

## Example Usage

The `id15_NXstress` package is designed to be flexible and easy to integrate into your projects. Below is a practical example of how to utilize the `initialize_system` function provided by the package to configure and execute the main functionality using a `config.yml` file located in the same directory as your script.

### Setting Up Your Environment

Before running the example, ensure that you have a `config.yml` file in the same directory as your script. This configuration file should contain all the necessary settings that `NXstressFromRaw` requires to operate correctly.

### Running the Example

Here's a simple script to demonstrate how to initialize and run the package:

```python
from easistrain.id15_NXstress import initialize_system
import os

# Determine the path to the configuration file
config_path = os.path.join(os.path.dirname(__file__), 'config.yml')

# Initialize the system with the configuration file
nx_stress = initialize_system(config_path)

# Execute the main functionality
nx_stress.main()
```

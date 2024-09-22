# FEM Battery solver

## Installation

Follow the steps in `Dockerfile` to get all required packages.

## Compilation

For users using `icx`, use the following command:

```cmake .. -DDEBUG=OFF -DUSE_MKL=ON -DCMAKE_CXX_COMPILER=icx-cc -DCMAKE_C_COMPILER=icx -DMKL_INTERFACE_FULL=intel_lp64```

For users using `gcc`:

```cmake .. -DDEBUG=OFF -DUSE_MKL=ON -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DMKL_INTERFACE_FULL=intel_lp64```

## Usage

Currently you need to create a `config.json` file. One sample of the file:

```json
{
    "tolerance": 1e-12,
    "dt": 10,
    "step": 350,
    "epsilon_e_an": 0.3,
    "epsilon_e_ca": 0.3,
    "epsilon_e_sep": 1,
    "epsilon_s_an": 0.6,
    "epsilon_s_ca": 0.5,
    "epsilon_s_sep": 0,
    "sigma_an": 100,
    "sigma_ca": 10,
    "sigma_sep": 0,
    "de_an": 2.7877e-10,
    "de_ca": 2.7877e-10,
    "de_sep": 2.7877e-10,
    "ds_an": 3.9e-14,
    "ds_ca": 1e-13,
    "ds_sep": 0,
    "r_p": 10e-6,
    "l_ref": 1e-6,
    "j_ref": 1e-3,
    "c_max_an": 24983,
    "c_max_ca": 51218,
    "c_int_an": 19624,
    "c_int_ca": 20046,
    "k_an": 1e-10,
    "k_ca": 3e-11,
    "ce_int": 1000,
    "bruggeman": 1.5,
    "trans": 0.4
}
```
| Status | Coverage |
| :----: | :----: |
| [![Build Status](https://travis-ci.org/rpgowers/NoisyNeuron.jl.svg?branch=master)](https://travis-ci.org/rpgowers/NoisyNeuron.jl) | [![codecov.io](http://codecov.io/github/rpgowers/NoisyNeuron.jl/coverage.svg?branch=master)](http://codecov.io/github/rpgowers/NoisyNeuron.jl?branch=master) |

# NoisyNeuron
Simulator package of single neurons subject to stochastic drive.

So far, one can only simulate point neurons when subjected to steady synaptic drive, current modulation, and variance modulation.

Cell parameters for the point neuron are entered with PointNeuron(E_L,τ_v).

| Status | Coverage |
| :----: | :----: |
| [![Build Status](https://travis-ci.org/rpgowers/NoisyNeuron.jl.svg?branch=master)](https://travis-ci.org/rpgowers/NoisyNeuron.jl) | [![codecov.io](http://codecov.io/github/rpgowers/NoisyNeuron.jl/coverage.svg?branch=master)](http://codecov.io/github/rpgowers/NoisyNeuron.jl?branch=master) |

# NoisyNeuron
Simulator package of single neurons subject to stochastic drive.

One can simulate point neurons when subjected to steady synaptic drive, current modulation, and variance modulation.

One can also simulate a sealed cable with steady synaptic drive.

Cell parameters for the point neuron are entered with PointNeuron(τ_v) and for a cable with Neurite(τ_v,λ,M,dx).
| Status | Coverage |
| :----: | :----: |
| [![Build Status](https://travis-ci.org/rpgowers/NoisyNeuron.svg?branch=master)](https://travis-ci.org/rpgowers/NoisyNeuron) | [![codecov.io](http://codecov.io/github/rpgowers/NoisyNeuron/coverage.svg?branch=master)](http://codecov.io/github/rpgowers/NoisyNeuron?branch=master) |

# NoisyNeuron
Simulator package of single neurons subject to stochastic drive.

One can simulate point neurons when subjected to steady synaptic drive, current modulation, and variance modulation.

One can also simulate a sealed cable with steady synaptic drive.

Cell parameters for the point neuron are entered with PointNeuron(τ_v) and for a cable with Neurite(τ_v,λ,M,dx).
# A Pre-Communication Mechanism for Evolutionary Multitasking Optimization

This repository contains the MATLAB implementation for the paper: *"A Pre-Communication Mechanism for Evolutionary Multitasking Optimization"*.

This work introduces a novel Pre-Communication Mechanism (PCM) designed to enhance the performance of existing Evolutionary Multitasking Optimization (EMTO) algorithms by leveraging prior information between different tasks during the early stages of evolution.

## About The Project

This project provides the implementation of the PCM framework and the scripts required to reproduce the experimental results presented in the PCM paper. The codebase is built upon the MTO-Platform (https://github.com/intLyc/MTO-Platform), and all necessary files from the platform are included in this repository for convenience.

The proposed PCM-enhanced algorithms are located in the `Algorithms/Multi-task/PCM` directory.

### Prerequisites

* MATLAB (R2021b or later)

## Getting Started

To get a local copy up and running, follow these simple steps.

### Installation

1.  Clone the repository to your local machine:
    ```sh
    git clone [https://github.com/lix1108/PCM-EMTO.git](https://github.com/lix1108/PCM-EMTO.git)
    ```
2.  Open MATLAB.
3.  Change the current directory in MATLAB to the cloned repository's root folder.

## Running the Experiments

To replicate all experiments from the paper, simply run the main script located in the project's root directory.

By running the script `run_pcm.m` in the root directory of the MATLAB project, the computation will be executed automatically, and the results will be saved to the `./Exp_data/` directory.

```matlab
% In MATLAB Command Window
run_pcm.m
```

The `run_pcm.m` script will automatically execute the experiments for all benchmark functions used in the paper, including:
* CEC2017-MTSO
* WCCI2020-MTSO
* The real-world Sensor Coverage Problem (SCP)
* The Different Task Combination Problems (DTCP)

The results will be organized into subfolders within `./Exp_data/` corresponding to each experiment. You can comment out lines in the `run_pcm.m` script if you wish to run only specific experiments.

## Citation

If you use this code or the PCM framework in your research, please cite the papers:

```bibtex
@article{Yang2025PCM,
  title={A pre-communication mechanism for evolutionary multitasking optimization},
  author={Yang, Cuicui and Li, Xiang and Ji, Junzhong and Zhang, Xiaoyu},
  journal={Neural Computing and Applications},
  year={2025},
  doi={10.1007/s00521-025-11534-6},
  publisher={Springer}
}

@misc{li2023mtop,
    title={MToP: A MATLAB Optimization Platform for Evolutionary Multitasking}, 
    author={Yanchi Li and Wenyin Gong and Fei Ming and Tingyu Zhang and Shuijia Li and Qiong Gu},
    year={2023},
    eprint={2312.08134},
    archivePrefix={arXiv},
    primaryClass={cs.NE},
    doi={10.48550/arXiv.2312.08134},
    url={https://arxiv.org/abs/2312.08134}
}
```
## Copyright and License

Copyright (c) 2025, lix1108 and all authors of the paper.

## Acknowledgments

This research code is implemented based on the MTO-Platform. We sincerely thank the original authors for their significant contribution to the field.

https://github.com/intLyc/MTO-Platform

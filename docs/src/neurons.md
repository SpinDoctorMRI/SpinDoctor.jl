# Neuron recipes

The diffusion MRI signal arising from neurons can be numerically simulated by solving the
Bloch-Torrey partial differential equation. In order to facilitate the diffusion MRI
simulation of realistic neurons by the research community, we constructed finite element
meshes for a group of 36 pyramidal neurons and a group of 29 spindle neurons whose
morphological descriptions were found in the publicly available neuron repository
[NeuroMopho.org](http://neuromorpho.org). These finite elements meshes range from having
15163 nodes to 622553 nodes. We also broke the neurons into the soma and dendrite branches
and created finite elements meshes for these cell components. Through the Neuron Module,
these neuron and components finite element meshes can be seamlessly coupled with the
functionalities of SpinDoctor to provide the diffusion MRI signal that can be attributed to
spins inside neurons.

The neuron mesh files are available at
[https://github.com/SpinDoctorMRI/NeuronMeshes](https://github.com/SpinDoctorMRI/NeuronMeshes).

The table below lists the names and the finite element mesh sizes of the group of 36
pyramidal neurons and the group of 29 spindle neurons.

| Neuron ID           | Num. of FE mesh nodes | Neuron ID            | Num. of FE mesh nodes |
| ------------------- | --------------------  | -------------------- | --------------------- |
| 03a\_spindle2aFI    | 38202                 | 02a\_pyramidal2aFI   | 119156                |
| 03a\_spindle6aFI    | 44000                 | 02b\_pyramidal1aACC  | 45216                 |
| 03b\_spindle4aACC   | 17370                 | 02b\_pyramidal1aFI   | 105384                |
| 03b\_spindle5aACC   | 26345                 | 03a\_pyramidal9aFI   | 81530                 |
| 03b\_spindle6aACC   | 26792                 | 03b\_pyramidal2aACC  | 28183                 |
| 03b\_spindle7aACC   | 21618                 | 03b\_pyramidal3aACC  | 27607                 |
| 04b\_spindle3aFI    | 51265                 | 03b\_pyramidal3aFI   | 151362                |
| 05b\_spindle5aFI    | 22457                 | 03b\_pyramidal4aFI   | 96177                 |
| 06b\_spindle8aACC   | 15163                 | 03b\_pyramidal9aFI   | 66162                 |
| 07b\_spindle9aACC   | 54952                 | 04a\_pyramidal4aACC  | 150897                |
| 08a\_spindle13aACC  | 46293                 | 04a\_pyramidal5aACC  | 89256                 |
| 09o\_spindle7aFI    | 38992                 | 04b\_pyramidal5aFI   | 95784                 |
| 09o\_spindle8aFI    | 60755                 | 04b\_pyramidal6aACC  | 87195                 |
| 10a\_spindle18aACC  | 25797                 | 04b\_pyramidal6aFI   | 90482                 |
| 12a\_spindle19aACC  | 31841                 | 04b\_pyramidal7aACC  | 622553                |
| 12o\_spindle9aFI    | 29320                 | 05a\_pyramidal10aACC | 201506                |
| 13o\_spindle10aFI   | 43081                 | 05a\_pyramidal8aACC  | 139975                |
| 15o\_spindle12aFI   | 101548                | 05b\_pyramidal7aFI   | 208203                |
| 16o\_spindle13aFI   | 18266                 | 05b\_pyramidal8aFI   | 124350                |
| 19o\_spindle14aFI   | 25786                 | 05b\_pyramidal9aACC  | 366659                |
| 21o\_spindle15aFI   | 28822                 | 06a\_pyramidal11aACC | 319574                |
| 23o\_spindle16aFI   | 30073                 | 06b\_pyramidal10aFI  | 106808                |
| 25o\_spindle17aFI   | 52919                 | 06b\_pyramidal12aACC | 277718                |
| 26o\_spindle18aFI   | 36239                 | 07a\_pyramidal13aACC | 155854                |
| 27o\_spindle19aFI   | 50807                 | 07b\_pyramidal14aACC | 309789                |
| 28o\_spindle20aFI   | 56036                 | 08o\_pyramidal11aFI  | 419651                |
| 28o\_spindle21aFI   | 17581                 | 10a\_pyramidal15aACC | 56184                 |
| 29o\_spindle22aFI   | 18414                 | 11a\_pyramidal16aACC | 222732                |
| 30o\_spindle23aFI   | 26357                 | 11o\_pyramidal12aFI  | 380293                |
| 22o\_pyramidal16aFI | 389878                | 17o\_pyramidal13aFI  | 326989                |
| 24o\_pyramidal17aFI | 245058                | 18o\_pyramidal14aFI  | 338453                |
| 25o\_pyramidal18aFI | 71209                 | 20o\_pyramidal15aFI  | 247116                |
| 31o\_pyramidal19aFI | 619390                |                      |                       |

The following table shows the morphological characteristics for the neurons.

| Neuron ID            | Brain region       | Avg. dia. ``(\mu m)`` | Height ``(\mu m)`` | Soma vol. ``(\mu m^3)`` | Total vol. ``(\mu m^3)`` |
| -------------------- | ------------------ | --------------------- | ------------------ | ----------------------- | ------------------------ |
| 02a\_pyramidal2aFI   | fronto-insula      | 1.27                  | 404.85             | 19701.05                | 25639.61                 |
| 02b\_pyramidal1aACC  | anterior cingulate | 1.58                  | 363.08             | 9065.56                 | 11579.71                 |
| 02b\_pyramidal1aFI   | fronto-insula      | 1.62                  | 381.56             | 22475.52                | 29804.96                 |
| 03a\_pyramidal9aFI   | fronto-insula      | 2.15                  | 532.30             | 22557.27                | 30189.28                 |
| 03a\_spindle2aFI     | fronto-insula      | 1.74                  | 387.16             | 13406.27                | 17684.23                 |
| 03a\_spindle6aFI     | fronto-insula      | 1.66                  | 501.47             | 33458.19                | 37812.72                 |
| 03b\_pyramidal2aACC  | anterior cingulate | 1.48                  | 189.29             | 2977.21                 | 4487.44                  |
| 03b\_pyramidal3aACC  | anterior cingulate | 1.14                  | 188.45             | 6005.06                 | 6891.04                  |
| 03b\_pyramidal3aFI   | fronto-insula      | 1.84                  | 496.35             | 32510.62                | 46154.08                 |
| 03b\_pyramidal4aFI   | fronto-insula      | 1.33                  | 414.70             | 35253.85                | 39324.87                 |
| 03b\_pyramidal9aFI   | fronto-insula      | 1.92                  | 430.06             | 15263.14                | 20532.57                 |
| 03b\_spindle4aACC    | anterior cingulate | 1.43                  | 336.33             | 3098.39                 | 4070.19                  |
| 03b\_spindle5aACC    | anterior cingulate | 1.49                  | 221.52             | 11925.78                | 13242.53                 |
| 03b\_spindle6aACC    | anterior cingulate | 1.33                  | 398.36             | 4027.74                 | 6058.67                  |
| 03b\_spindle7aACC    | anterior cingulate | 1.18                  | 369.51             | 4982.41                 | 6076.52                  |
| 04a\_pyramidal4aACC  | anterior cingulate | 1.52                  | 705.96             | 5684.55                 | 13637.33                 |
| 04a\_pyramidal5aACC  | anterior cingulate | 1.86                  | 410.59             | 15010.43                | 24648.35                 |
| 04b\_pyramidal5aFI   | fronto-insula      | 1.78                  | 480.13             | 10312.87                | 17184.26                 |
| 04b\_pyramidal6aACC  | anterior cingulate | 1.41                  | 465.88             | 3129.97                 | 7497.12                  |
| 04b\_pyramidal6aFI   | fronto-insula      | 1.56                  | 310.21             | 14718.05                | 21708.21                 |
| 04b\_pyramidal7aACC  | anterior cingulate | 1.3                   | 610.42             | 17060.60                | 28552.49                 |
| 04b\_spindle3aFI     | fronto-insula      | 2.71                  | 391.14             | 22569.99                | 28404.13                 |
| 05a\_pyramidal10aACC | anterior cingulate | 1.74                  | 281.20             | 16604.06                | 21826.41                 |
| 05a\_pyramidal8aACC  | anterior cingulate | 1.29                  | 430.37             | 24709.77                | 29778.79                 |
| 05b\_pyramidal7aFI   | fronto-insula      | 2.18                  | 281.02             | 25720.11                | 32731.05                 |
| 05b\_pyramidal8aFI   | fronto-insula      | 1.66                  | 361.45             | 32527.06                | 44679.46                 |
| 05b\_pyramidal9aACC  | anterior cingulate | 1.56                  | 650.60             | 23948.05                | 40014.54                 |
| 05b\_spindle5aFI     | fronto-insula      | 2.35                  | 381.88             | 15383.08                | 18190.63                 |
| 06a\_pyramidal11aACC | anterior cingulate | 1.46                  | 437.60             | 17222.02                | 29995.02                 |
| 06b\_pyramidal10aFI  | fronto-insula      | 1.92                  | 365.18             | 43127.81                | 52179.53                 |
| 06b\_pyramidal12aACC | anterior cingulate | 1.52                  | 324.94             | 17181.33                | 24931.32                 |
| 06b\_spindle8aACC    | anterior cingulate | 1.92                  | 342.21             | 18237.49                | 19462.92                 |
| 07a\_pyramidal13aACC | anterior cingulate | 1.37                  | 325.73             | 6254.53                 | 8738.01                  |
| 07b\_pyramidal14aACC | anterior cingulate | 1.67                  | 350.40             | 16053.07                | 22772.96                 |
| 07b\_spindle9aACC    | anterior cingulate | 1.75                  | 437.87             | 21344.83                | 27307.48                 |
| 08a\_spindle13aACC   | anterior cingulate | 1.74                  | 814.45             | 9911.07                 | 14113.32                 |
| 08o\_pyramidal11aFI  | fronto-insula      | 1.91                  | 421.68             | 11512.38                | 24326.94                 |
| 09o\_spindle7aFI     | fronto-insula      | 2.90                  | 472.87             | 22052.10                | 27905.89                 |
| 09o\_spindle8aFI     | fronto-insula      | 2.05                  | 376.73             | 11923.76                | 15189.32                 |
| 10a\_pyramidal15aACC | anterior cingulate | 1.40                  | 341.48             | 8522.11                 | 10960.84                 |
| 10a\_spindle18aACC   | anterior cingulate | 1.57                  | 457.90             | 5895.17                 | 7219.28                  |
| 11a\_pyramidal16aACC | anterior cingulate | 1.27                  | 486.31             | 8807.01                 | 12263.84                 |
| 11o\_pyramidal12aFI  | fronto-insula      | 1.91                  | 369.34             | 70786.62                | 79516.92                 |
| 12a\_spindle19aACC   | anterior cingulate | 2.05                  | 431.22             | 12178.08                | 15618.67                 |
| 12o\_spindle9aFI     | fronto-insula      | 3.41                  | 305.31             | 29983.79                | 36678.18                 |
| 13o\_spindle10aFI    | fronto-insula      | 2.69                  | 516.92             | 39866.55                | 46022.15                 |
| 15o\_spindle12aFI    | fronto-insula      | 3.60                  | 604.57             | 53192.65                | 79170.43                 |
| 16o\_spindle13aFI    | fronto-insula      | 2.17                  | 364.66             | 17467.88                | 18888.13                 |
| 17o\_pyramidal13aFI  | fronto-insula      | 1.89                  | 340.77             | 11004.30                | 21167.19                 |
| 18o\_pyramidal14aFI  | fronto-insula      | 1.74                  | 288.41             | 69851.56                | 78999.20                 |
| 19o\_spindle14aFI    | fronto-insula      | 2.18                  | 232.21             | 10507.15                | 12905.43                 |
| 20o\_pyramidal15aFI  | fronto-insula      | 1.82                  | 383.18             | 22344.32                | 27667.19                 |
| 21o\_spindle15aFI    | fronto-insula      | 2.36                  | 286.33             | 17567.69                | 29466.53                 |
| 22o\_pyramidal16aFI  | fronto-insula      | 1.94                  | 585.35             | 18776.05                | 29441.43                 |
| 23o\_spindle16aFI    | fronto-insula      | 1.67                  | 420.05             | 10429.13                | 13482.93                 |
| 24o\_pyramidal17aFI  | fronto-insula      | 2.04                  | 371.99             | 40986.40                | 47377.09                 |
| 25o\_pyramidal18aFI  | fronto-insula      | 1.80                  | 364.05             | 18587.13                | 23572.15                 |
| 25o\_spindle17aFI    | fronto-insula      | 1.79                  | 358.70             | 7897.44                 | 13563.26                 |
| 26o\_spindle18aFI    | fronto-insula      | 2.27                  | 442.65             | 52911.93                | 56084.44                 |
| 27o\_spindle19aFI    | fronto-insula      | 1.73                  | 275.08             | 20640.14                | 25423.96                 |
| 28o\_spindle20aFI    | fronto-insula      | 3.00                  | 520.69             | 35442.59                | 51267.07                 |
| 28o\_spindle21aFI    | fronto-insula      | 2.62                  | 298.57             | 35579.06                | 37783.31                 |
| 29o\_spindle22aFI    | fronto-insula      | 3.52                  | 402.84             | 62928.22                | 83279.12                 |
| 31o\_pyramidal19aFI  | fronto-insula      | 2.26                  | 303.55             | 65950.80                | 86376.72                 |

The neuron models and the measurement data are from (Ascoli 2007) and (Watson 2006).

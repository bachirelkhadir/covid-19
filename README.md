# Introduction

We study the following vector field modeling the transmission of COVID-19.

![f(S,E,I_r,I_u) = \begin{pmatrix}         \beta  \frac S N I_r - \mu  \beta  \frac S N I_u\\         \beta  \frac S N I_r + \mu  \beta  \frac S N I_u - \frac E Z\\         \alpha  \frac E Z - \frac {I_r} D\\         (1-\alpha)  \frac E Z - \frac {I_u} D \end{pmatrix}](https://render.githubusercontent.com/render/math?math=f(S%2CE%2CI_r%2CI_u)%20%3D%20%5Cbegin%7Bpmatrix%7D%20%20%20%20%20%20%20%20%20%5Cbeta%20%20%5Cfrac%20S%20N%20I_r%20-%20%5Cmu%20%20%5Cbeta%20%20%5Cfrac%20S%20N%20I_u%5C%5C%20%20%20%20%20%20%20%20%20%5Cbeta%20%20%5Cfrac%20S%20N%20I_r%20%2B%20%5Cmu%20%20%5Cbeta%20%20%5Cfrac%20S%20N%20I_u%20-%20%5Cfrac%20E%20Z%5C%5C%20%20%20%20%20%20%20%20%20%5Calpha%20%20%5Cfrac%20E%20Z%20-%20%5Cfrac%20%7BI_r%7D%20D%5C%5C%20%20%20%20%20%20%20%20%20(1-%5Calpha)%20%20%5Cfrac%20E%20Z%20-%20%5Cfrac%20%7BI_u%7D%20D%20%5Cend%7Bpmatrix%7D)

![\beta = 1.0, \mu=.8, \alpha=.1](https://render.githubusercontent.com/render/math?math=%5Cbeta%20%3D%201.0%2C%20%5Cmu%3D.8%2C%20%5Calpha%3D.1)
![D = Z = 4 \text{days}](https://render.githubusercontent.com/render/math?math=D%20%3D%20Z%20%3D%204%20%5Ctext%7Bdays%7D)
![$N = 1$](https://render.githubusercontent.com/render/math?math=%24N%20%3D%201%24)


# To run the code
Make sure you already have Julia installed on your system. Clone this repo

```
git clone https://github.com/bachirelkhadir/covid-19.git
```

Move inside the folder and run julia using the following commands.
```
cd covid-19/
julia --project
```

Finally, run the code.
```
include("covid_19.jl")
```


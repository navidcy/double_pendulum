# Double Pendulum

Julia scripts to make double pendulum animations. Dynamics include damping and forcing.

![logo](https://user-images.githubusercontent.com/7112768/122657073-866b8580-d1a3-11eb-88cc-27189a85d37c.gif)


## Instructions

First you need to [install Julia](https://julialang.org/downloads/). We suggest using Julia version 1.6 or later.

Then clone the repository, e.g.,

```
git clone https://github.com/navidcy/double_pendulum.git
```

Enter the directory you've cloned the repository in, e.g., 

```
cd double_pendulum
```

and run

```
julia --project= -e 'using Pkg; Pkg.instantiate()'
```

to download all required dependencies. Afterwards, run

```
julia --project double_pendulum_logo.jl
```

to produce the logo above or 

```
julia --project double_pendulum.jl
```

to produce a series of animations demonstrating the normal modes of the double pendulum without forcing nor dissipation and how we can excite these normal modes in the presence of forcing and dissipation.

# Final Project

Andrew Boyacigiller, Cecily Chung, Connor Flick 

## How to Run

The project runs using the same format as Projects 3 and 4 for this course. JSON scene files have extra features, listed below.

Global:
- `skybox` is added to specify a texture used to create a skysphere. Consumes a filepath.
- `textureU` and `textureV` are added to specify how many times `skybox` should be tiled in the sphere. Consumes a float.

Material:
- `solid` indicates whether an object is a solid object or a non-solid object, as a boolean.
- `density` indicates how dense a non-solid object should be. Consumes a float. Densities above 100 are treated as solid objects. Best results are with a density between 0 and 2. 


## Outstanding Issues

- The star generation consumes a significant amount of memory in some circumstances on some platforms, particularly on a Windows release build. 
- No black hole :(
- Edges of solid objects may not always be "captured" by raymarching, making them uneven.

## Andrew Boyacigiller

### Running

The current version of my code is on the black-hole branch. The latest commit is using a non physically-based approach however previous commits are using the full pipeline. The pipeline is loosely based on the designs outlined in [this](https://20k.github.io/c++/2024/05/31/schwarzschild.html) document.

### Sources

1. https://profoundphysics.com/christoffel-symbols-a-complete-guide-with-examples/
1. https://arxiv.org/pdf/1511.06025
1. https://iopscience.iop.org/article/10.3847/1538-4365/aac9ca/pdf
1. https://iopscience.iop.org/article/10.1088/0264-9381/32/6/065001/pdf
1. https://en.wikipedia.org/wiki/Kerr_metric
1. https://20k.github.io/c++/2024/05/31/schwarzschild.html
1. https://www2.mpia-hd.mpg.de/homes/tmueller/pdfs/catalogue_2014-05-21.pdf
1. https://arxiv.org/pdf/gr-qc/0204066
1. https://en.wikipedia.org/wiki/Geodesics_in_general_relativity
1. https://physics.stackexchange.com/questions/733433/christoffel-symbols-for-schwarzschild-metric
1. https://en.wikipedia.org/wiki/Schwarzschild_metric

## Connor Flick

Much of my focus for this project was on raymarching. I used the following sources, along with what was discussed in class, to guide my thinking:

### Sources

1. https://www.scratchapixel.com/lessons/3d-basic-rendering/volume-rendering-for-developers/intro-volume-rendering.html
1. https://cglearn.eu/pub/advanced-computer-graphics/volumetric-rendering

### Assets

I pulled the space skybox asset (`spaceagain.png`) from Adobe Stock under their licensing terms. The asset is labelled as AI generated. 


# Flintknapping Simulation
This project is intended to simulate the process of flint-knapping (forming/shaping certain types of stone to manufacture tools and other objects). This is primarily done by striking a prepared “core” at certain angles to produce flakes of stone. The flaking process follows a largely predictable conchoidal fracture phenomenon known as a Hertzian Cone, which will be the focus of this project. The final goal is to simulate the generation of a flake from a sample core using a Hertzian Cone calculated from a variable input point, angle, and force.

## Team
Allen Briggs

## Introduction
Flint-knapping is the process of creating stone tools through breaking flakes off of a stone core. Stones frequently used include flint, chert, obsidian, and others. This is done in 3 major ways - hard hammer percussion flaking (for large pieces), soft hammer percussion flaking (for smaller, more controllable flakes that still remove a great deal of material), and pressure flaking (for minute details). In ancient times this technique was widely used for manufacturing stone tools, and it was used up until the mid-twentieth century for manufacturing gun flints. In modern times it is the domain of hobbyists and archeologists studing stone technology.

Due to the predictable nature of the Hertzian Cone, which is relatively easy to generate geometrically based on the position and angle of a given strike, the actual simulation can be considered as a geometric Boolean operation. A 3D Boolean operation can be a Difference, a Union, or an Intersection, each of which is a type of geometry calculated when two other geometries intersect in 3D space. The two I focussed on for this are the Intersection and Difference. The Intersection of two geometries is the geometry representing where the two objects overlap - it can be thought of as the central oval of a venn diagram. With regards to this project, the Intersection of the Hertzian Cone with the stone Core is the flake produced by the strike. The Difference is the modification of one of the original geometries to include only the space that does not overlap with the other geometry. In this case, the Core after the flake has been removed is the result of a Difference of the Core by the Hertzian Cone. 

### Project Use Cases
The primary use of this project is in visualizing and teaching the process to novices. The inspiriation came from taking an Experimental Archeology class in which we learned the process. Visualizing the outcome of a particular strike was often difficult for me, and it would be useful to be able to interactively see how the flaking works without breaking real rock.

## Implementation
To cut past the work of setting up an environment in which to simulate an interactive 3D simulation I chose to use the Unity Engine as the basic environment for my project. Unity provides a simple development environment, robust renderer, and easy access to internal mesh data, which is essential to the chosen algorithm. Efficient booleans on triangle meshes (most 3D renderers, including Unity, work in triangle representations) are still an area of active research, and existing methods can be complex and computationally expensive. I followed the approach outlined by Mei and Tipper in [Simple and Robust Boolean Operations for Triangulated Surfaces](https://arxiv.org/pdf/1308.4434.pdf). This algorithm proceeds with the following series of general steps, operating on two input meshes:

1. For each triangle in the Core mesh, if it possibly intersects with a triangle in the Cone mesh, compute the line of intersection between the two triangles in 3D space. Keep a record of each triangle associated to the lines of intersection that apply to it (so one intersection calculation result would be applied to both triangles involved). For this step I used Moller's Fast Triangle-Triangle Intersection test, following his own adaptations for computing the intersection line.
2. Using the intersection lines, divide each triangle into smaller polygons. Triangulate these newly formed polygons into smaller triangles. If both meshes were properly made, there should be no open polygons (which would be invalid for triangulation). Replace each triangle in the original meshes with the list of triangles created during this process.
3. Calculate subsurfaces of the geometries by separating them at edges where they begin or end an overlap. These subsurfaces are then combined into the geometry for any boolean type. At this time, I take the Difference subsections to output the model of the Core after the flake has been removed, since this is the simplest of the two outputs to calculate.

## Demo
Provided is a demonstration build of the simulation in its current state. Work is intended to continue on this project, focusing immediately in the areas of user interaction (allowing the user to set up the strike the way they want it), followed by generating the flake in conjunction with the modified Core. As it stands now, there is only a static input for the strike, and that can only be triggered once. The static input is configured to cut a slice out of the bottom of the Core, which is the large, conical object that solely inhabits the simulation when launched.

The Demo can be found [here](https://people.rit.edu/arb7109/csci716/project/build/).

### Controls
1. The camera can be controlled with standard WASD keys, rotated by holding down the right mouse button and dragging, and moved up and down with the Spacebar and Left Shift keys, respectively.
2. The Strike can be triggered by pressing the 'I' key. A good indication of its success is that part of the right of the Core will have hard edges differently, instead of being smoothly shaded. The bottom will also have changed to have a conical dent/cutout, though with Unity's default shading that is difficult to notice. While a wireframe shading would be more ideal, I was not able to get one working to apply to this project.

## Final Presentation
[A presentation on the results of this project may be found here.](https://www.youtube.com/watch?v=Izvm0YxuIx0)

## References
[Knapping - Wikipedia](https://en.wikipedia.org/wiki/Knapping)

[Hertzian Cone - Wikipedia](https://en.wikipedia.org/wiki/Hertzian_cone)

[Flint-knapping with the Hertzian Cone](https://www.researchgate.net/figure/a-Flint-knapping-utilising-the-Hertzian-cone-phenomenon-to-remove-long-sharp-pieces-of_fig2_326838071)

[The Formation of Flakes](https://www-jstor-org.ezproxy.rit.edu/stable/281)

[Generating Mesh Geometry in Unity](https://docs.unity3d.com/Manual/GeneratingMeshGeometryProcedurally.html)

[Simple and Robust Boolean Operations for Triangulated Surfaces](https://arxiv.org/pdf/1308.4434.pdf)

[A Fast Triangle-Triangle Intersection Test](https://web.stanford.edu/class/cs277/resources/papers/Moller1997b.pdf)

[Triangulation By Ear-Clipping]()

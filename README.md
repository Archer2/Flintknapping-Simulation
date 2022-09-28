# Flintknapping Simulation
Unity simulation of basic percussion flintknapping. Uses an approximation of a Hertzian Cone to generate a flake from a core, with the simulated material being based off of obsidian.
## Team
Allen Briggs
## Description
This project is intended to simulate the process of flint-knapping (forming/shaping certain types of stone to manufacture tools and other objects). This is primarily done by striking a prepared “core” at certain angles to produce flakes of stone. The flaking process follows a largely predictable conchoidal fracture phenomenon known as a Hertzian Cone, which will be the focus of this project. The final goal is to simulate the generation of a flake from a sample core using a Hertzian Cone calculated from a variable input point, angle, and force.
### Flint-knapping
Flint-knapping is the process of creating stone tools through breaking flakes off of a stone core. Stones frequently used include flint, chert, obsidian, and others. This is done in 3 major ways - hard hammer percussion flaking (for large pieces), soft hammer percussion flaking (for smaller, more controllable flakes that still remove a great deal of material), and pressure flaking (for minute details). In ancient times this technique was widely used for manufacturing stone tools, and it was used up until the mid-twentieth century for manufacturing gun flints. In modern times it is the domain of hobbyists and archeologists studing stone technology.
### Project Use Cases
The primary use of this project is in visualizing and teaching the process to novices. The inspiriation came from taking an Experimental Archeology class in which we learned the process. Visualizing the outcome of a particular strike was difficult for me, and it would be useful to be able to interactively see how the flaking works without breaking real rock.
## Project Timeline

| Week | Tasks |
| --- | --- |
| 5 | Background research |
| 6 | More research, project setup, and creation of initial "Core" model |
| 7 | Extract the Core's geometrical representation into a useable format |
| 8 | Generate the Hertzian Cone structure and generate geometry from input |
| 9 | Develop boolean algorithm to “cut” the Core |
| 10 | Midway Presentation - complete algorithm implementation |
| 11 | Refine functionality and begin implementing user input |
| 12 | Finish implementing user input |
| 13 | Develop project webpage |
| 14 | Finish webpage, final presentation |

## References
[Knapping - Wikipedia](https://en.wikipedia.org/wiki/Knapping)

[Hertzian Cone - Wikipedia](https://en.wikipedia.org/wiki/Hertzian_cone)

[Flint-knapping with the Hertzian Cone](https://www.researchgate.net/figure/a-Flint-knapping-utilising-the-Hertzian-cone-phenomenon-to-remove-long-sharp-pieces-of_fig2_326838071)

[The Formation of Flakes](https://www-jstor-org.ezproxy.rit.edu/stable/281)

[Generating Mesh Geometry in Unity](https://docs.unity3d.com/Manual/GeneratingMeshGeometryProcedurally.html)

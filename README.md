# Im2mesh (2D or 3D image to finite element mesh)



**Im2mesh** is an open-source MATLAB/Octave package for generating finite element mesh based on 2D or 3D multi-phase image. It provides a robust workflow capable of processing various input images, such as microstructure images of engineering materials. Due to its generalized framework, Im2mesh can handle segmented image with more than 10 phases.  Im2mesh was originally released on [MathWorks File Exchange](https://www.mathworks.com/matlabcentral/fileexchange/71772-im2mesh-2d-image-to-finite-element-mesh) in 2019. 

Im2mesh can also be used as a mesh generation interface for MATLAB 2D multi-part geometry, aka multi-domain or multi-phase geometry (see demo12-18).

<p align="center">
  <img src = "https://mjx888.github.io/im2mesh_demo_html/cover_2602.jpg" height="220"> 
</p>


**Downloads:**

- [Im2mesh package](https://github.com/mjx888/im2mesh/releases)
- [GUI version](https://mjx888.github.io/others/Im2mesh_GUI.mlappinstall) (MATLAB app)
- [GUI version](https://mjx888.github.io/others/Installer_Im2mesh_GUI.zip) (standalone desktop application & no need to install MATLAB)



**News:**

- Version 2.60 can generate tetrahedral mesh based on 3D voxel image! See [demo20](https://mjx888.github.io/im2mesh_demo_html/demo20.html).
- Version 2.45 can export image boundaries as `dxf` file (CAD).
- Version 2.2.0 can use Gmsh as mesh generator (quadrilateral mesh).



**Features (for 2d):**

- Accurately preserve the contact details between different phases. 
- Incorporates polyline smoothing and simplification
- Able to edit polygonal boundary before mesh generation.
- Support phase selection and local mesh refinement.
- 4 mesh generators are available for selection: [MESH2D](https://github.com/dengwirda/mesh2d), [generateMesh](https://www.mathworks.com/help/pde/ug/pde.pdemodel.generatemesh.html), [Gmsh](https://gmsh.info/), and pixelMesh.
- Graphical user interface (GUI) version is available as a MATLAB app and as a standalone desktop application.


<p align="center">
  <img src = "https://mjx888.github.io/im2mesh_demo_html/GUI.png" height="300"> 
</p>



**Generated mesh can be exported as:** 

- `inp` file with boundary node set (Abaqus)
- `bdf` file (Nastran bulk data, compatible with COMSOL)
- `msh` file (Gmsh mesh format)
- `stl` file
- For other formats, you can import the generated `msh` file into software Gmsh and then export.



## Dependencies

- When running demo02 and demo18 in Im2mesh package, you need to install MATLAB Partial Differential Equation Toolbox. For other demos, no need to install any MATLAB toolboxes.
- When running demo20 in Im2mesh package, you need to install fTetWild. 
- When using Im2mesh_GUI as a standalone desktop application, there is no need to install MATLAB or any MATLAB toolboxes. 



## Version compatibility

- Im2mesh_GUI: MATLAB R2017b or later; version higher than R2018b is preferred.
- Im2mesh package: MATLAB R2017b or later. GNU Octave 9.3.0 or later.
- Gmsh: tested with version 4.13 and 4.14.




## How to start

After downloading Im2mesh package, I suggest you start with [Im2mesh_GUI app](https://github.com/mjx888/im2mesh/tree/main/Im2mesh_GUI%20app) in the folder, which will help you understand the workflow and parameters of Im2mesh. A detailed tutorial is provided in [Im2mesh_GUI Tutorial.pdf](https://github.com/mjx888/im2mesh/blob/main/Im2mesh_GUI%20Tutorial.pdf). **Note that Im2mesh_GUI is for 2D images.** 

Then, you can learn to use Im2mesh package in the folder "Im2mesh_Matlab" or "Im2mesh_Octave". 20 examples are provided. demo01-18 are for 2D images. demo19-20 are for 3D voxel images.

- If you're using MATLAB,  examples are live script `mlx` files (`demo01.mlx` ~ `demo20.mlx`). If you find some text in the `mlx` file is missing, please read the `html` file instead.
- If you're using Octave,  examples are `m` files (`demo01.m` ~ `demo10.m`).
- Examples are also available as `html` files in the folder "demo_html".
- You can skip demo04 and demo09. These demos are kept for historical reason.
- If you're only interested in 3D voxel images, you can skip demo01-18.

**Examples:**

- [demo01](https://mjx888.github.io/im2mesh_demo_html/demo01.html) - Demonstrate function `im2mesh`, which use `MESH2D` as mesh generator.
- [demo02](https://mjx888.github.io/im2mesh_demo_html/demo02.html) - Demonstrate function `im2meshBuiltIn`, which use MATLAB built-in function `generateMesh` as mesh generator.
- [demo03](https://mjx888.github.io/im2mesh_demo_html/demo03.html) - Export: save mesh as `inp`, `bdf`, `msh` or, `stl` file; save image boundary as `dxf` file, Gmsh `geo` file, or PSLG data.
- [demo04](https://mjx888.github.io/im2mesh_demo_html/demo04.html) - What is inside function `im2mesh`
- [demo05](https://mjx888.github.io/im2mesh_demo_html/demo05.html) - Thresholds in polyline smoothing
- [demo06](https://mjx888.github.io/im2mesh_demo_html/demo06.html) - Parameter `hmax` and `grad_limit` in mesh generation
- [demo07](https://mjx888.github.io/im2mesh_demo_html/demo07.html) - Function `plotMeshes`
- [demo08](https://mjx888.github.io/im2mesh_demo_html/demo08.html) - Select phases for meshing
- [demo09](https://mjx888.github.io/im2mesh_demo_html/demo09.html) - Find node sets at the interface and boundary
- [demo10](https://mjx888.github.io/im2mesh_demo_html/demo10.html) - Function `pixelMesh` (pixel-based quadrilateral mesh)
- [demo11](https://mjx888.github.io/im2mesh_demo_html/demo11.html) - Use `Gmsh` as mesh generator
- [demo12](https://mjx888.github.io/im2mesh_demo_html/demo12.html) - Use polyshape to define geometry
- [demo13](https://mjx888.github.io/im2mesh_demo_html/demo13.html) - 2D mesh for periodic boundary conditions
- [demo14](https://mjx888.github.io/im2mesh_demo_html/demo14.html) - Edit polygonal boundaries before meshing
- [demo15](https://mjx888.github.io/im2mesh_demo_html/demo15.html) - Generate mesh based on 2D contours
- [demo16](https://mjx888.github.io/im2mesh_demo_html/demo16.html) - Add mesh seeds/nodes
- [demo17](https://mjx888.github.io/im2mesh_demo_html/demo17.html) - Refine mesh
- [demo18](https://mjx888.github.io/im2mesh_demo_html/demo18.html) - 2D image to tetrahedral mesh
- [demo19](https://mjx888.github.io/im2mesh_demo_html/demo19.html) - Function `voxelMesh` (voxel-based hexahedral mesh)
- [demo20](https://mjx888.github.io/im2mesh_demo_html/demo20.html) - 3D voxel image to tetrahedral mesh (via fTetWild)



## Author

Jiexian Ma

## Cite as

If you use Im2mesh, please cite it as follows. You can click the DOI link below to see other citation styles.

Ma, J., & Li, Y. (2025). Im2mesh: A MATLAB/Octave package for generating finite element mesh based on 2D multi-phase image (2.1.5). Zenodo. https://doi.org/10.5281/zenodo.14847059

Once my paper is published, I will update a new DOI here.

## Acknowledgments

Many thanks to Dr. Yang Lu for providing valuable suggestions and testing of export formats. 

This project incorporates code from the following open-source projects. I appreciate the contributions of the original authors. Each incorporated code retains its original copyright.

- [MESH2D](https://github.com/dengwirda/mesh2d) by Darren Engwirda
- [dpsimplify](https://www.mathworks.com/matlabcentral/fileexchange/21132-line-simplification) by Wolfgang Schwanghart
- [p_poly_dist](https://www.mathworks.com/matlabcentral/fileexchange/12744-distance-from-points-to-polyline-or-polygon) by Michael Yoshpe
- [MeshQualityQuads](https://www.mathworks.com/matlabcentral/fileexchange/33108-unstructured-quadrilateral-mesh-quality-assessment) by Allan Peter Engsig-Karup
- [XtalMesh](https://github.com/jonathanhestroffer/XtalMesh) by Jonathan Hestroffer



## Other related projects

- [writeMesh (write mesh to inp, bdf, and msh files)](https://github.com/mjx888/writeMesh)


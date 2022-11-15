using System;
using System.Collections;
using System.Collections.Generic;
using TMPro.EditorUtilities;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.UIElements;

public class Striker : MonoBehaviour
{
    protected MeshFilter coneMeshFilter;
    [SerializeField] protected int tesselationCount = 10;
    [SerializeField] protected GameObject coreModel;
    protected Mesh coreMesh;

    // Start is called before the first frame update
    void Start()
    {
        // Find any random vertex on the bottom of the core, and move to that vertex location
        if(coreModel != null && coreMesh != null)
        {
            Vector3 min = new Vector3(0, float.MaxValue, 0);
            for(int i = 0; i < coreMesh.vertices.Length; i++)
            {
                if(coreMesh.vertices[i].y < min.y)
                {
                    min = coreMesh.vertices[i];
                }
            }

            transform.position = coreModel.transform.TransformPoint(min);
            transform.rotation = Quaternion.Euler(0.0f, 0.0f, 31.0f); // Hardcoded to cut a chunk out of the bottom of the default cone
        }
    }

    // Awake is called when the script instance is being loaded
    void Awake()
    {
        // Get a reference to the MeshFilter's Mesh
        coneMeshFilter = GetComponent<MeshFilter>();
        if (coneMeshFilter == null)
        {
            coneMeshFilter = gameObject.AddComponent<MeshFilter>();
        }
        coneMeshFilter.mesh.Clear(); // Reset any existing mesh data

        // Make sure there is a MeshRenderer
        if (GetComponent<MeshRenderer>() == null)
        {
            gameObject.AddComponent<MeshRenderer>();
        }

        // If we don't have a Core set, try and find one
        if (coreModel == null)
        {
            coreModel = GameObject.Find("Basic_Core");
        }
        // A bit wacky syntax, but the null state could change inside the previous if{}
        if(coreModel != null)
        {
            MeshFilter core = coreModel.GetComponent<MeshFilter>();
            if(core != null)
            {
                coreMesh = core.mesh; // Could still be null afterwards - CHECK BEFORE USING (TODO: Default to a simple cube?)
            }
        }
    }

    // Update is called once per frame
    void Update()
    {
        // Sample input action for now
        if (Input.GetKeyDown(KeyCode.I) && coneMeshFilter.mesh.name == gameObject.name) // Check that the name is not the default name assigned (even when mesh is set to null)
        {
            Mesh cone = GenerateCone();
            cone.name = "Hertzian_Cone";
            coneMeshFilter.mesh = cone;

            //SutherlandHodgmanAlgorithm/*CullOnly*/(coreMesh, cone);
        }
    }

    // Generates a simple cone mesh (the Hertzian Cone) from tesellation data. Returns
    // the created mesh for use in the algorithm
    //  - TODO: Implement Impact Force instead of a basic cone height
    //  - TODO: Instead of generating the whole mesh, probably only vertices are needed.
    //          Mesh fully formed for now for visualization
    Mesh GenerateCone()
    {
        Mesh cone = new Mesh();

        List<Vector3> vertexList = new List<Vector3>();
        List<Vector3> normalList = new List<Vector3>();
        List<Vector2> uvList = new List<Vector2>();
        List<int> indexList = new List<int>();

        // Designed for solid-color material (uvs are 00, 10, 11)
        Action<Vector3, Vector3, Vector3, Vector3> AddTriangle = 
            (Vector3 vert1, Vector3 vert2, Vector3 vert3, Vector3 norm) =>
        {
            vertexList.Add(vert1); vertexList.Add(vert2); vertexList.Add(vert3);
            uvList.Add(new Vector2(0, 0)); uvList.Add(new Vector2(1, 0)); uvList.Add(new Vector2(1, 1));
            
            // Increment indices 3 times, add the same normal to each vertex
            for(int i = 0; i < 3; i++)
            {
                indexList.Add(indexList.Count);
                normalList.Add(norm);
            }
        };

        // Create Cone
        if (tesselationCount < 3) tesselationCount = 3;

        float tanSixty = Mathf.Tan(Mathf.PI / 3f); // 60 degree angle
        float height = 1.8f; // Constant for now, must be greater than 1 to go through sample Core - TODO: Depend on Force --- IF THIS CHANGES MUST CHANGE CONSTANT IN ALGORITHM'S IsPointWithinCone() HELPER!
        float radius = height * tanSixty;
        Vector3 peak = new Vector3(0, 0, 0); // Cone origin is at this position, and expands downwards
        Vector3 valley = new Vector3(0, -height, 0);

        float angle = 0;
        for (int i = 0; i < tesselationCount; i++)
        {
            float nextAngle = angle + Mathf.PI * 2f / tesselationCount;

            float cosAngle = Mathf.Cos(angle); // Chache for normal use
            float sinAngle = Mathf.Sin(angle); // Cache for normal use
            Vector3 p1 = new Vector3(radius * cosAngle, valley.y, radius * sinAngle);
            Vector3 p2 = new Vector3(radius * Mathf.Cos(nextAngle), valley.y, radius * Mathf.Sin(nextAngle));

            Vector3 normal = Vector3.Normalize(new Vector3(cosAngle, tanSixty, sinAngle)); // Requires normalizing with just x and z first, but basic cos/sin is already normalized

            AddTriangle(peak, p2, p1, normal);
            AddTriangle(valley, p1, p2, Vector3.down);

            angle = nextAngle;
        }

        // Assign data to Mesh
        cone.vertices = vertexList.ToArray();
        cone.normals = normalList.ToArray();
        cone.uv = uvList.ToArray();
        cone.triangles = indexList.ToArray();

        cone.RecalculateNormals();
        cone.RecalculateBounds();

        return cone;
    }

    // "Normalizes" or Unpacks a Mesh's data from efficient indexed arrays of data to basic
    // arrays with duplicated data, where the "triangles" array is simply a linearly increasing
    // value. While this is horribly inefficient and creates potentially several times the
    // number of vertices, it makes solving the Boolean algorithms much easier
    //  - This will be removed if possible, but at the moment I'm not sure I can actually solve
    //    the problem without the data unpacked
    void NormalizeMesh(Mesh mesh)
    {
        List<Vector3> newVertices = new List<Vector3>();
        List<Vector3> newNormals = new List<Vector3>();
        List<Vector4> newTangents = new List<Vector4>();
        List<Vector2> newUVs = new List<Vector2>();
        List<int> newIndices = new List<int>();

        // For every index/triangle, extract its data to new array, duplicating any data optimized by an index buffer
        foreach(int index in mesh.triangles)
        {
            newVertices.Add(mesh.vertices[index]);
            newNormals.Add(mesh.normals[index]);
            newTangents.Add(mesh.tangents[index]);
            newUVs.Add(mesh.uv[index]);
            newIndices.Add(newIndices.Count); // Incremental
        }

        mesh.vertices = newVertices.ToArray();
        mesh.normals = newNormals.ToArray();
        mesh.tangents = newTangents.ToArray();
        mesh.uv = newUVs.ToArray();
        mesh.triangles = newIndices.ToArray();
    }

    // Performs an implementation of the Sutherland-Hodgman Clip Algorithm in 3D space
    // using planes instead of lines. The "clipped" polygon is the Boolean Intersection, which
    // is then differenced from the Subject (Core) Mesh.
    // This algorithm was chosen for 2 main reasons: 1: It looked fairly simple in 3D compared to
    // listed alternative algorithms (and was stated to be simpler); and 2: It's input is just a simple
    // list of vertices, which is fairly trivial to obtain from a Mesh
    //  - Currently this ONLY functions if the ClipMesh is a Cone generated by GenerateCone()
    //  - Other implementations of better algorithms will be considered if this fails
    //      - Greiner-Hormann is the next bet, but does not handle vertex intersections or common edges
    //  - subjectMesh: core mesh to operate on - this mesh will be clipped
    //  - clipMesh: mesh to clip the subject with - this MUST be generated by GenerateCone() at this time
    void SutherlandHodgmanAlgorithm(Mesh subjectMesh, Mesh clipMesh)
    {
        NormalizeMesh(subjectMesh); // Expand Mesh data from indexed to facilitate algorithm (yes, this is horribly inefficient)

        List<Vector3> outputVertices = new List<Vector3>();
        subjectMesh.GetVertices(outputVertices); // Starting output vertices are the subject's vertices

        // Checks if a given cone-space point is within the volume of the cone defined by GenerateCone()
        //  - This fully assumes the clipMesh is that Cone, because a simple bounds check is not sufficient
        //  - This uses constant height and radius from the GenerateCone() method - TODO: Make those class-level member vars
        //  - Code is taken from https://stackoverflow.com/questions/12826117/
        Func<Vector3, bool> IsPointWithinCone = (Vector3 coneSpaceVert) =>
        {
            Vector3 dir = Vector3.down; // Direction from Cone Peak to Base
            Vector3 coneSpacePeak = Vector3.zero; // Local space base of the Cone was created at origin
            float distanceFromPeak = Vector3.Dot(coneSpaceVert - coneSpacePeak, dir); // Dot Peak->Point with Peak->Base
            if (distanceFromPeak < 0 || distanceFromPeak > 1.8f) // 1.8 IS THE HEIGHT CONSTANT FROM CONE GENERATION - MUST BE CHANGED IF THAT IS
            {
                return false; // Not within vertical range of cone
            }
            float baseRadius = 1.8f * Mathf.Tan(Mathf.PI / 3f); // 1.8 IS THE HEIGHT CONSTANT FROM CONE GENERATION - MUST BE CHANGED IF THAT IS
            float radius = distanceFromPeak / 1.8f * baseRadius; // 1.8 IS THE HEIGHT CONSTANT FROM CONE GENERATION - MUST BE CHANGED IF THAT IS
            float distanceFromCenter = Vector3.Magnitude((coneSpaceVert - coneSpacePeak) - distanceFromPeak * dir);

            return distanceFromCenter < radius;
        };

        // Test if a Vertex is within the Triangle defined by 3 vertices
        // Code taken from https://stackoverflow.com/questions/2049582/ (C# answer by Glenn Slayden 2013, 2021) --- Because the test Vertex has already been projected onto the plane, y can be removed and it can become a test of the projection onto (x,z)
        // Everything in world space
        // This just flat out does not return well for anything except the peak vertex, which is ON THE BOUNDING TRIANGLE
        Func<Vector3, Vector3, Vector3, Vector3, bool> IsPointWithinBoundingTriangle =
            (Vector3 vertexOnPlane, Vector3 boundingVert0, Vector3 boundingVert1, Vector3 boundingVert2) =>
            {
                //Debug.Log("Testing Point to Triangle:");
                //Debug.Log(vertexOnPlane);
                //Debug.Log(boundingVert0);
                //Debug.Log(boundingVert1);
                //Debug.Log(boundingVert2);
                float s = (boundingVert0.x - boundingVert2.x) * (vertexOnPlane.z - boundingVert2.z) - (boundingVert0.z - boundingVert2.z) * (vertexOnPlane.x - boundingVert2.x);
                float t = (boundingVert1.x - boundingVert0.x) * (vertexOnPlane.z - boundingVert0.z) - (boundingVert1.z - boundingVert0.z) * (vertexOnPlane.x - boundingVert0.x);

                if ((s < 0) != (t < 0) && s != 0 && t != 0)
                {
                    //Debug.Log(false);
                    //Debug.Log("----------------------------");
                    return false;
                }

                float d = (boundingVert2.x - boundingVert1.x) * (vertexOnPlane.z - boundingVert1.z) - (boundingVert2.z - boundingVert1.z) * (vertexOnPlane.x - boundingVert1.x);
                bool ret = d == 0 || (d < 0) == (s + t <= 0);

                //Debug.Log(ret);
                //Debug.Log("----------------------------");
                return ret;
            };

        // Iterate over every second Triangle of the Clip Mesh - all bottom
        // triangles have even indices and the same normal plane - check last.
        // This should be every triangle around the cone
        for (int i = 0; i < clipMesh.vertices.Length; i += 6)
        {
            Vector3 peakPoint = transform.TransformPoint(clipMesh.vertices[i]); // ClipMesh is the Cone, so use this transform to get to World Space
            Vector3 planePoint1 = transform.TransformPoint(clipMesh.vertices[i + 1]); // The Triangle's other Vertices are points in the Plane
            Vector3 planePoint2 = transform.TransformPoint(clipMesh.vertices[i + 2]); // The Triangle's other Vertices are points in the Plane

            Plane clipPlane = new Plane(peakPoint, planePoint1, planePoint2); // Construct plane

            // Iterate over all still-valid points on the subjectMesh, removing any points that are not on the correct side of the clip Plane
            List<Vector3> newOutputVertices = new List<Vector3>();
            for (int j = 0; j < outputVertices.Count; j += 3)  // 3 Vertices for a single Triangle
            {
                List<Vector3> validPoints = new List<Vector3>();
                List<Vector3> invalidPoints = new List<Vector3>();

                // Go through Triangle and mark each of the 3 Vertices as Valid or Invalid, based on whether or not they are
                // on the normal side of the Triangle and within it's infinite vertical extents
                for (int k = 0; k < 3; k++)
                {
                    Vector3 vertex = coreModel.transform.TransformPoint(outputVertices[j+k]); // outputVertices comes from subjectMesh, which is the Core - FOR THIS TO WORK SUBJECTMESH MUST = COREMESH

                    // If this vertex is NOT on the side that the normal faces, remove it. To account for normal Planes extending
                    // in every direction, also verify that the point is within the clip Cone. A Point is only invalid if it is
                    // both on the negative side of the current Plane AND if it is within the Cone's volume.
                    if (!clipPlane.GetSide(vertex) // If we are on the negative side of the Plane
                        && IsPointWithinCone(transform.InverseTransformPoint(vertex))) // Old version
                        //&& IsPointWithinBoundingTriangle(clipPlane.ClosestPointOnPlane(vertex), peakPoint, transform.TransformPoint(clipMesh.vertices[i+1]), transform.TransformPoint(clipMesh.vertices[i+2]))) // ClipMesh is the Cone, so use this transform to get from World to Clip Local Space
                    {
                        invalidPoints.Add(coreModel.transform.InverseTransformPoint(vertex)); // Convert back to Core space now
                    }
                    else
                    {
                        validPoints.Add(coreModel.transform.InverseTransformPoint(vertex)); // Convert back to Core space now
                    }
                }

                // Perform each possible Case of vertex replacement based on the number of invalid Vertices
                switch(invalidPoints.Count)
                {
                    case 0: // All 3 Vertices are valid
                        if(validPoints.Count != 3)
                        {
                            Debug.Log($"So the invalid count is 0, but the valid count is {validPoints.Count} instead of 3?");
                        }
                        newOutputVertices.AddRange(validPoints);
                        break;
                    case 1: // 2 Vertices valid, 1 invalid - must create 2 new Triangles <----------------------------------------- CASE VERY BROKEN!!!!!!!!!!!!!!!!!!!
                        //newOutputVertices.AddRange(validPoints); // <------ Placeholder Visualization
                        //newOutputVertices.Add(invalidPoints[0]); // <------ Placeholder Visualization
                        //break;

                        List<Vector3> verticesToAdd = new List<Vector3>(4); // Will hold 2 old valid points and 2 new points. Renormalization means order does not matter. Indices 0 and 2 will be duplicated between Tris
                        verticesToAdd.AddRange(validPoints); // Add valid Points

                        // For each valid Vertex, cast towards the invalid Vertex. The intersection of that Ray with the
                        // Triangle (if present) is the location of the new Vertex in that direction. A new Vertex should
                        // be found for each Valid Vertex, totalling 4
                        foreach(Vector3 vertex in validPoints)
                        {
                            Vector3 dir = Vector3.Normalize(vertex - invalidPoints[0]); // Direction doesn't care about core vs. world space
                            Ray castToValid = new Ray(coreModel.transform.TransformPoint(invalidPoints[0]), dir); // World Space for Plane
                            float distance = 0.0f;
                            clipPlane.Raycast(castToValid, out distance);
                            Vector3 worldIntersectPoint = castToValid.GetPoint(distance);

                            // Check if the new vertex is within the bounds of the Triangle (works because it should already be on the plane that all 3 Vertices are on)
                            if (//IsPointWithinBoundingTriangle(worldIntersectPoint, peakPoint, transform.TransformPoint(clipMesh.vertices[i + 1]), transform.TransformPoint(clipMesh.vertices[i + 2])) // Function does not work?
                                /*&& */distance <= (vertex - invalidPoints[0]).magnitude)
                            {
                                verticesToAdd.Add(coreModel.transform.InverseTransformPoint(worldIntersectPoint)); // Cast World to local Core transform, as it is now a Core Vertex
                            }
                        }

                        // If neither test replacement Vertex is valid, just use the bad one (TODO: DON'T?)
                        if(verticesToAdd.Count == 2)
                        {
                            verticesToAdd.Add(invalidPoints[0]);
                        }

                        // If there is only 1 new Triangle just add that one - More work needs to be done
                        // In reality this means that the Clip Triangle only partially intersected this one, and 1 of the
                        // "Invalid" Vertices is really still Valid, so more work needs to be done to figure out where
                        // to connect that to the New Valid Vertex
                        if(verticesToAdd.Count == 3)
                        {
                            newOutputVertices.AddRange(verticesToAdd);
                        }
                        else if(verticesToAdd.Count == 4) // If there are 2 new Triangles then 2 vertices are shared between them
                        {
                            // Triangle 1
                            newOutputVertices.Add(verticesToAdd[0]);
                            newOutputVertices.Add(verticesToAdd[1]);
                            newOutputVertices.Add(verticesToAdd[2]);

                            // Triangle 2
                            newOutputVertices.Add(verticesToAdd[0]);
                            newOutputVertices.Add(verticesToAdd[2]);
                            newOutputVertices.Add(verticesToAdd[3]);
                        }
                        else
                        {
                            Debug.Log("What? There should be 3 vertices, or 4 if 2 triangles were created");
                        }
                        break;
                    case 2: // 1 Vertex valid, 2 invalid - still output 1 Triangle <------------------------------------- CASE WORKING (ish? No clue what those bottom Triangles are doing)
                        newOutputVertices.Add(validPoints[0]); // Add valid point

                        // For each invalid Vertex move it up to a point on the Clip Plane. Ignore if
                        // the new point would fall outside the Triangle causing this Clip (it would
                        // not be on the clip mesh)
                        foreach (Vector3 vertex in invalidPoints) // vertex is in core space
                        {
                            Vector3 dir = Vector3.Normalize(vertex - validPoints[0]); // Direction doesn't care about core vs. world space
                            Ray castToInvalid = new Ray(coreModel.transform.TransformPoint(validPoints[0]), dir); // World Space for Plane
                            float distance = 0.0f;
                            clipPlane.Raycast(castToInvalid, out distance);
                            Vector3 worldIntersectPoint = castToInvalid.GetPoint(distance); // World Space

                            // Check if the new vertex is within the bounds of the Triangle (works because it should already be on the plane that all 3 Vertices are on)
                            if (IsPointWithinBoundingTriangle(worldIntersectPoint, peakPoint, transform.TransformPoint(clipMesh.vertices[i + 1]), transform.TransformPoint(clipMesh.vertices[i + 2])))
                            {
                                newOutputVertices.Add(coreModel.transform.InverseTransformPoint(worldIntersectPoint)); // Cast World to local Core transform, as it is now a Core Vertex
                            }
                            else
                            {
                                newOutputVertices.Add(vertex); // TODO: DON'T DO THIS
                            }
                        }

                        // Swap winding order to face outward if required
                        int size = newOutputVertices.Count;
                        Vector3 newNormal = Vector3.Normalize(Vector3.Cross(newOutputVertices[size-2] - newOutputVertices[size-3], newOutputVertices[size-1] - newOutputVertices[size-3]));
                        Vector3 toCoreCenter = coreModel.transform.position - ((newOutputVertices[size-3] + newOutputVertices[size-2] + newOutputVertices[size-1]) / 3.0f); // Core position - Triangle centroid
                        if(Vector3.Dot(newNormal, toCoreCenter) > 0) // If the dot product is above 0 then the normal is facing somewhat towards the center of the mesh - which it should not be
                        {
                            Vector3 tempVert = newOutputVertices[size-1];
                            newOutputVertices[size-1] = newOutputVertices[size-2];
                            newOutputVertices[size-2] = tempVert;
                        }
                        break;
                    case 3: // All 3 Vertices are invalid
                        // Placeholder - This should add 1-2 Triangles along the Plane instead of it not existing at all?
                        break;
                    default:
                        Debug.Log($"How the hell are there not 0-3 Vertices in the Invalid Array? ({invalidPoints.Count})");
                        break;
                }

                // No more operations for this Triangle. It has been added with all modifications needed
            }

            outputVertices = newOutputVertices; // Set output for next cycle
        }
        // TODO - Check bottom face!!!!

        // Reset Indices and other data for the subject mesh
        List<int> newIndices = new List<int>();
        for(int i = 0; i < outputVertices.Count; i++)
        {
            newIndices.Add(i);
        }
        subjectMesh.triangles = newIndices.ToArray();
        subjectMesh.vertices = outputVertices.ToArray();
        subjectMesh.RecalculateNormals();
        subjectMesh.RecalculateTangents();
    }

    // Performs an implementation of the Sutherland-Hodgman Clip Algorithm in 3D space
    // using planes instead of lines. The "clipped" polygon is the Boolean Intersection, which
    // is then differenced from the Subject (Core) Mesh.
    // This algorithm was chosen for 2 main reasons: 1: It looked fairly simple in 3D compared to
    // listed alternative algorithms (and was stated to be simpler); and 2: It's input is just a simple
    // list of vertices, which is fairly trivial to obtain from a Mesh
    //  - Currently this ONLY functions if the ClipMesh is a Cone generated by GenerateCone()
    //  - Other implementations of better algorithms will be considered if this fails
    //      - Greiner-Hormann is the next bet, but does not handle vertex intersections or common edges
    //  - subjectMesh: core mesh to operate on - this mesh will be clipped
    //  - clipMesh: mesh to clip the subject with - this MUST be generated by GenerateCone() at this time
    void SutherlandHodgmanAlgorithmCullOnly(Mesh subjectMesh, Mesh clipMesh)
    {
        NormalizeMesh(subjectMesh); // Expand Mesh data from indexed to facilitate algorithm (yes this is horribly inefficient)

        // The Triangles array is just indexes into the vertex array, so by culling Triangles (groups of 3 indices) from
        // this array the vertices of the mesh are effectively unaltered
        //  Ex: Removing (3, 4, 5) from (0, 1, 2, 3, 4, 5, 6, 7, 8) will make the Mesh have 2 triangles, even though
        //      enough vertices exist for 3
        List<int> outputIndices = new List<int>();
        subjectMesh.GetTriangles(outputIndices, 0); // Start with output as the entirety of the subjectMesh

        // Checks if a given cone-space point is within the volume of the cone defined by GenerateCone()
        //  - This fully assumes the clipMesh is that Cone, because a simple bounds check is not sufficient
        //  - This uses constant height and radius from the GenerateCone() method - TODO: Make those class-level member vars
        //  - Code is taken from https://stackoverflow.com/questions/12826117/
        Func<Vector3, bool> IsPointWithinCone = (Vector3 coneSpaceVert) =>
        {
            Vector3 dir = Vector3.down; // Direction from Cone Peak to Base
            Vector3 coneSpacePeak = Vector3.zero; // Local space base of the Cone was created at origin
            float distanceFromPeak = Vector3.Dot(coneSpaceVert - coneSpacePeak, dir); // Dot Peak->Point with Peak->Base
            if (distanceFromPeak < 0 || distanceFromPeak > 1.8f) // 1.8 IS THE HEIGHT CONSTANT FROM CONE GENERATION - MUST BE CHANGED IF THAT IS
            {
                return false; // Not within vertical range of cone
            }
            float baseRadius = 1.8f * Mathf.Tan(Mathf.PI / 3f); // 1.8 IS THE HEIGHT CONSTANT FROM CONE GENERATION - MUST BE CHANGED IF THAT IS
            float radius = distanceFromPeak / 1.8f * baseRadius; // 1.8 IS THE HEIGHT CONSTANT FROM CONE GENERATION - MUST BE CHANGED IF THAT IS
            float distanceFromCenter = Vector3.Magnitude((coneSpaceVert - coneSpacePeak) - distanceFromPeak * dir);

            return distanceFromCenter < radius;
        };

        // Iterate over every second Triangle of the Clip Mesh - all bottom
        // triangles have even indices and the same normal plane - check last.
        // This should be every triangle around the cone
        for (int i = 0; i < clipMesh.vertices.Length; i += 6)
        {
            Vector3 peakPoint = clipMesh.vertices[i]; // Should be same for all planes, which is fine since normal is different
            peakPoint = transform.TransformPoint(peakPoint); // ClipMesh is the Cone, so use this transform to get to World Space
            Vector3 planeNormal = clipMesh.normals[i]; // Should be = for i+1 and i+2

            Plane clipPlane = new Plane(planeNormal, peakPoint);

            // Iterate over all still-valid points on the subjectMesh, removing any points that are not on the correct side of the clip Plane
            for(int j = 0; j < outputIndices.Count; j++)
            {
                Vector3 vertex = coreModel.transform.TransformPoint(subjectMesh.vertices[outputIndices[j]]); // outputVertices comes from subjectMesh, which is the Core - FOR THIS TO WORK SUBJECTMESH MUST = COREMESH

                // If this vertex is NOT on the side that the normal faces, remove it. To account for normal Planes extending
                // in every direction, also verify that the point is within the clip Cone. A Point is only invalid if it is
                // both on the negative side of the current Plane AND if it is within the Cone's volume.
                if(!clipPlane.GetSide(vertex) && IsPointWithinCone(transform.InverseTransformPoint(vertex))) // ClipMesh is the Cone, so use this transform to get from World to Clip Local Space
                {
                    int k = (j / 3) * 3; // Truncate to let k = the first index of the triangle containing j (j = 3 -> k = 3, j = 4 -> k = 3, j = 5 -> k = 3, j = 6 -> k = 6, etc.)
                    outputIndices.RemoveAt(k);
                    outputIndices.RemoveAt(k); // Removal subtracts 1 from all subsequent indices, so the formerly "next" index is the current index
                    outputIndices.RemoveAt(k); // Removal subtracts 1 from all subsequent indices, so the formerly "next" index is the current index
                    j -= (j % 3) + 1; // Revert to index before the triangle that was deleted (j = 0 -> j = -1, j = 2 -> j = -1, j = 3 -> j = 0, etc.)
                }
            }
        }
        // TODO - Check bottom face!!!!

        subjectMesh.triangles = outputIndices.ToArray();
    }
}

// Try this:
// 1: Normals and Tangents can be Recalculated with Mesh functions, so erase completely
//    UVs can be ignored - solid color material (too bad, uv calculations are WAY out of scope)
//    This leaves only Vertices and Triangles/Indices to deal with
//      a: This is absolutely essential, because modifying indices allows the removal of existing
//         values, but does NOT allow for the addition of new vertices
// 2: Foreach Clip Plane, start with an empty list of VERTICES (the new ouput), and an old output
//    list of VERTICES from the previous stage (inital is all vertices in the Subject Mesh).
//    Indices do not matter here, but the Mesh must have been normalized so that every 3 Vertices
//    forms a single Triangle. Process every Subject Triangle (set of 3 Vertices), finding which points
//    lie on the negative (!GetSide(v)) of the Clip Plane. There are 4 possible cases:
//      a: 0 Invalid Points: This case is easy, keep the entire Triangle. Add all 3 Vertices to
//         the new ouput list
//      b: 3 Invalid Points: Also easy, but the reverse of Case A. The entire triangle is within
//         the clip polygon, so remove it. Add no Vertices to the new output list
//      c: 2 Invalid Points: Find the Valid Point, and Raycast (twice) along the directions of that   <--- Plane.Raycast() accepts a Ray (origin and direction) and returns
//         towards each invalid Point. The intersection of that Raycast with the Clip Plane           <--- the distance along the Ray from its origin to the Plane
//         determines the replacement Vertex in that direction. Since 2 Vertices are removed and
//         2 are added, there is still just 1 Triangle. Add the Valid Vertex and the new Vertices
//         to the new output list in the same winding order
//      d: 1 Invalid Point: Find the Invalid Point, and Raycast from that Point to each valid Point.
//         The intersection of that Ray and the Clip Plane provides 2 new Vertices. These new
//         Vertices and the 2 valid Vertices form 2 new Triangles. Organize the 4 Vertices into 2
//         Triangles (accounting for winding order), and add those 6 Vertices to the new output list.
//         2 Vertices will be duplicated between the Triangles - 1 of these will be a new Vertex and
//         1 will be an old, valid Vertex
//  3: After the main run is completed the final iterations "new output" will consist of the Vertices
//     of the output Mesh. Run through this array, incrementing an Index counter at every step for the
//     Mesh's new Triangle array (start from scratch at empty and 0). Since so many Vertices and Triangles
//     are added and removed during this process, it is best to redo it from scratch. Once complete call
//     Mesh.RecalculateNormals() and Mesh.RecalculateTangents(). Iterate across a new array for UVs and
//     set them all to some arbitrary value (0.5 and 0.5 works). This should not matter with solid-color
//     materials, and performing UV alterations seems EXTREMELY complex (maybe just lerp them with the new
//     Triangles? Doesn't seem that simple)
//    
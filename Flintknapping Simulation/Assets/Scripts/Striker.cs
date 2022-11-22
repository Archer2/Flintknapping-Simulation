using System;
using System.Collections;
using System.Collections.Generic;
using TMPro.EditorUtilities;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.UIElements;

// Struct Representing all necessary Vertex Data to be tracked for individual
// Vertices in the Striking operation. Saves needing many separate arrays for
// data that should be grouped, added, altered, and removed together
struct Vertex
{
    public Vector3 Position;
    public Vector3 Normal; // Should end up unused - too complex to fragment into more Triangles
    public Vector3 Tangent; // Should end up unused - too complex to fragment into more Triangles
    public Vector2 UV; // Should end up unused - too complex to fragment into more Triangles
}

// Struct Representing an undirected Edge between two Points in space. This can be stored with
// a Triangle to determine how it fragments during the intersection (step 2 of Mei and Tipper).
// Edge Loops can be formed where for any Edge, P1 is P0 of the next Edge in the Loop
//  - Directions should be implemented, but are difficult to do with Moller
struct Edge
{
    public Vector3 P0; // Head, if being used as Directed
    public Vector3 P1; // Tail, if being used as Directed
}

// Struct Representing a single Mesh Triangle. This consists of 3 Vertices, and
// supports Triangle operations like planar normal and intersections
struct Triangle
{
    // 3 Vertices to make a triangle
    public Vertex V0; public Vertex V1; public Vertex V2;

    // Calculate the Triangle's Normal - this is DIFFERENT from each Vertex Normal,
    // it is the Normal of the Plane the Triangle lies on
    public Vector3 CalculateNormal()
    {
        Vector3 v1 = V1.Position - V0.Position;
        Vector3 v2 = V2.Position - V0.Position;
        return Vector3.Cross(v1, v2).normalized; // Assumes Winding Order in Clockwise
    }

    // Performs Moller's(97) Fast Intersection Test, returning the Edge
    // upon which the Triangles Intersect
    // This assumes both Triangles are within the same coordinate system - 
    // ideally pre-transformed into World Units
    //  - Performing other.Intersect(this) returns an Edge with the same points
    //    in reverse order
    //
    // Returns a boolean indicating whether or not there is an intersecion. If an intersection
    // is found, o_intersection is also set to an edge representing this intersection
    public bool Intersect(Triangle other, ref Edge o_intersection)
    {
        // -- Step 1 - Compute Plane equation of Other Triangle (T2)
        Vector3 N2 = other.CalculateNormal();
        float d2 = Vector3.Dot(-N2, other.V0.Position);
        Plane plane2 = new Plane(N2, d2); // Create a Unity Plane for access to functions like GetDistanceToPoint, etc.

        // Distances of T1 (this) vertices to T2 Plane (other)
        // Naming convention is d2v1i
        //  - d - Distance to
        //  - 2 - Triangle 2
        //  - v1 - Vertex of Triangle 1
        //  - i - Vertex # of Triangle 1
        float d2v10 = plane2.GetDistanceToPoint(V0.Position);
        float d2v11 = plane2.GetDistanceToPoint(V1.Position);
        float d2v12 = plane2.GetDistanceToPoint(V2.Position);

        // -- Step 2 - Reject as trivial if all points of This (T1) are on the same side
        if (d2v10 != 0 && d2v11 != 0 && d2v12 != 0) // All distances are not 0 (no points on the plane)
        {
            float d2v10Sign = Mathf.Sign(d2v10);
            if (d2v10Sign == Mathf.Sign(d2v11) && d2v10Sign == Mathf.Sign(d2v12)) // If all share the same sign
            {
                Debug.Log("All Vertices of Triangle 1 are on the same side of the Plane defined by Triangle 2");
                return false;
            }
        }
        else if (d2v10 == 0 && d2v11 == 0 && d2v12 == 0) // Distance from plane to all points is 0 - Triangles are co-planar
        {
            Debug.Log("Triangles are co-planar. Handle case separately");
            return false; // Placeholder
        }

        // Reorganize Vertices and floats (d2v1i) to have the first and third be on the same side of the P2
        // Modifying Vertex ordering will not have any negative effects, so long as the float distances are updated too
        {
            float d2v10Sign = Mathf.Sign(d2v10);
            float d2v11Sign = Mathf.Sign(d2v11);
            float d2v12Sign = Mathf.Sign(d2v12);

            // If the first and third signs are not equal, something must be swapped
            if (d2v10Sign != d2v12Sign)
            {
                if (d2v11Sign == d2v12Sign) // second and third are on same side, so swap first and second
                {
                    Vertex tmpV = V0; // Swap Vertices
                    V0 = V1;
                    V1 = tmpV;

                    float tmpF = d2v10; // Swap distances calculated from the Vertices
                    d2v10 = d2v11;
                    d2v11 = tmpF;
                }
                else if (d2v11Sign == d2v10Sign) // first and second are on same side, so put third in second
                {
                    Vertex tmpV = V2; // Swap Vertices
                    V2 = V1;
                    V1 = tmpV;

                    float tmpF = d2v12; // Swap distances calculated from the Vertices
                    d2v12 = d2v11;
                    d2v11 = tmpF;
                }
            }
        }

        // -- Step 3 - Compute Plane equation of This Triangle (T1)
        Vector3 N1 = CalculateNormal();
        float d1 = Vector3.Dot(-N1, V0.Position);
        Plane plane1 = new Plane(N1, d1);

        // Distances of T2 (other) vertices to T1 Plane (this)
        // Naming convention is d1v2i
        //  - d1 - Distance to Triangle 1
        //  - v2 - Vertex of Triangle 2
        //  - i - Vertex # of Triangle 2
        float d1v20 = plane1.GetDistanceToPoint(other.V0.Position);
        float d1v21 = plane1.GetDistanceToPoint(other.V1.Position);
        float d1v22 = plane1.GetDistanceToPoint(other.V2.Position);

        // -- Step 4 - Reject as trivial if all points of Other (T2) are on the same side
        if (d1v20 != 0 && d1v21 != 0 && d1v22 != 0) // All distances are not 0 (no points on the plane)
        {
            float d1v20Sign = Mathf.Sign(d1v20);
            if (d1v20Sign == Mathf.Sign(d1v21) && d1v20Sign == Mathf.Sign(d1v22)) // If all share the same sign
            {
                Debug.Log("All Vertices of Triangle 2 are on the same side of the Plane defined by Triangle 2");
                return false;
            }
        }
        // No need to check co-planar again - it would have been caught beforehand

        // Reorganize Vertices and floats (d2v1i) to have the first and third be on the same side of the P2
        // Modifying Vertex ordering will not have any negative effects, so long as the float distances are updated too
        {
            float d1v20Sign = Mathf.Sign(d1v20);
            float d1v21Sign = Mathf.Sign(d1v21);
            float d1v22Sign = Mathf.Sign(d1v22);

            // If the first and third signs are not equal, something must be swapped
            if (d1v20Sign != d1v22Sign)
            {
                if (d1v21Sign == d1v22Sign) // second and third are on same side, so swap first and second
                {
                    Vertex tmpV = other.V0; // Swap Vertices
                    other.V0 = other.V1;
                    other.V1 = tmpV;

                    float tmpF = d1v20; // Swap distances calculated from the Vertices
                    d1v20 = d1v21;
                    d1v21 = tmpF;
                }
                else if (d1v21Sign == d1v20Sign) // first and second are on same side, so put third in second
                {
                    Vertex tmpV = other.V2; // Swap Vertices
                    other.V2 = other.V1;
                    other.V1 = tmpV;

                    float tmpF = d1v22; // Swap distances calculated from the Vertices
                    d1v22 = d1v21;
                    d1v21 = tmpF;
                }
            }
        }

        // -- Step 5 - Compute intersection line and project onto largest axis - Projection onto Largest axis is not done,
        // as I do not fully understand the notation described and how the simplification works
        Vector3 D = Vector3.Cross(N1, N2); // D in the equation L = O + tD (O is ignored, since Translation to Origin does not affect result)

        // Project Triangle 1 (this) points
        // Naming convention is pv1i
        //  - pv1 - Projection Vertex of Triangle 1
        //  - i - Vertex # of Triangle 1
        float pv10 = Vector3.Dot(D, V0.Position);
        float pv11 = Vector3.Dot(D, V1.Position);
        float pv12 = Vector3.Dot(D, V2.Position);

        // Repeat for Triangle 2 (other) points, same convention
        float pv20 = Vector3.Dot(D, other.V0.Position);
        float pv21 = Vector3.Dot(D, other.V1.Position);
        float pv22 = Vector3.Dot(D, other.V2.Position);

        // -- Step 6 - Compute intervals for each Triangle
        float t11 = (pv11 - pv10) * (d2v10 / (d2v10 - d2v11)) + pv10;
        float t12 = (pv11 - pv12) * (d2v12 / (d2v12 - d2v11)) + pv12;

        float t21 = (pv21 - pv20) * (d1v20 / (d1v20 - d1v21)) + pv20;
        float t22 = (pv21 - pv22) * (d1v22 / (d1v22 - d1v21)) + pv22;

        Debug.Log($"\nT11: {t11}, T12: {t12}\nT21: {t21}, T22: {t22}");

        // -- Step 7 - Intersect the intervals
        // Ensure t11 < t12 and t21 < t22
        float tmp;
        if (t11 > t12)
        {
            tmp = t11;
            t11 = t12;
            t12 = tmp;
        }
        if (t21 > t22)
        {
            tmp = t21;
            t21 = t22;
            t22 = tmp;
        }

        if (t12 < t21 || t22 < t11)
        {
            Debug.Log("Distance between shortest and longest is greater than combined interval length, so they can't overlap");
            return false;
        }

        // Find actual intersection edges
        List<float> lengths = new List<float>(4);

        // Insert each tij value and sort in increasing order
        lengths.Add(t11);
        lengths.Add(t12);
        lengths.Add(t21);
        lengths.Add(t22);
        lengths.Sort();
        foreach(float v in lengths)
        {
            Debug.Log(v);
        }
        Debug.Log(D);

        o_intersection.P0 = D * lengths[1]; // Actual intersection is from index 1 to 2
        o_intersection.P1 = D * lengths[2];
        Debug.Log("INTERSECTION FOUND");
        return true;
    }
}

// Actual Unity Script
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
            StrikeMesh(cone);
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


    // Performs an implementation of the Boolean Operations on Triangulated Surfaces described by
    // Mei and Tipper (2013)
    // The general steps (as laid out in their paper) are as follows:
    //  1 - Take 2 Surface Mesh inputs (clipMesh and coreMesh, clipMesh MUST have been generated by GenerateCone())
    //  2 - Search for pairs of intersecting Triangles
    //  3 - Compute Intersection and re-triangulate
    //  4 - Merge and Update Surface Meshes (in this case only the coreMesh)
    //  5 - Form Intersection Loops (???)
    //  6 - Create sub-surfaces (???)
    //  7 - Assemble and distinguish sub-blocks (???)
    //  8 - Output sub-surfaces and/or sub-blocks (in this case just 1, but I'm not sure which - I believe an Intersection sub-block)
    void StrikeMesh(Mesh clipMesh)
    {
        Vertex v0 = new Vertex(), v1 = new Vertex(), v2 = new Vertex();
        v0.Position = new Vector3(2.0f, -2.0f, 0.0f);
        v1.Position = new Vector3(-2.0f, -2.0f, 0.0f);
        v2.Position = new Vector3(0.0f, 2.0f, 0.0f);
        Triangle t1 = new Triangle();
        t1.V0 = v0;
        t1.V1 = v1;
        t1.V2 = v2;

        Vertex u0 = new Vertex(), u1 = new Vertex(), u2 = new Vertex();
        u0.Position = new Vector3(0.0f, 0.0f, 1.0f);
        u1.Position = new Vector3(1.0f, 0.0f, -1.0f);
        u2.Position = new Vector3(-1.0f, 0.0f, -1.0f);
        Triangle t2 = new Triangle();
        t2.V0 = u0;
        t2.V1 = u1;
        t2.V2 = u2;

        Edge e = new Edge();
        t1.Intersect(t2, ref e);
        Debug.Log($"Intersection Edge: {e.P0}, {e.P1}");
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
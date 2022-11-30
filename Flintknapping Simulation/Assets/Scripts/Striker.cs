using System;
using System.Collections;
using System.Collections.Generic;
using TMPro.EditorUtilities;
using Unity.VisualScripting;
using UnityEditor.Profiling.Memory.Experimental;
using UnityEngine;
using UnityEngine.UI;
using UnityEngine.UIElements;

// Struct Representing all necessary Vertex Data to be tracked for individual
// Vertices in the Striking operation. Saves needing many separate arrays for
// data that should be grouped, added, altered, and removed together
struct Vertex
{
    public Vector3 Position;
    public Vector3 Normal; // Should end up unused - too complex to fragment into more Triangles
    public Vector4 Tangent; // Should end up unused - too complex to fragment into more Triangles
    public Vector2 UV; // Should end up unused - too complex to fragment into more Triangles
}

// Struct Representing an undirected Edge between two Points in space. This can be stored with
// a Triangle to determine how it fragments during the intersection (step 2 of Mei and Tipper).
// Edge Loops can be formed where for any Edge, P1 is P0 of the next Edge in the Loop
//  - Directions should be implemented, but are difficult to do with Moller
class Edge
{
    public Vector3 P0; // Head, if being used as Directed
    public Vector3 P1; // Tail, if being used as Directed

    // Checks if a given Point is on this Edge. Returns false if the
    // Point is an Endpoint
    public bool IsPointOnEdge(Vector3 testPoint)
    {
        //if (testPoint == P0 || testPoint == P1) return false;
        //double alpha = (testPoint.x - P0.x) / (P1.x - P0.x);
        //return (alpha < 0) ? false : (alpha > 1) ? false : true;

        Vector3 p0p1 = P1 - P0;
        Vector3 p0test = testPoint - P0;

        Vector3 cross = Vector3.Cross(p0p1, p0test);
        if (cross.sqrMagnitude != 0) return false; // Cross Product magnitude of 0 - points are colliniar

        float dotEdge = Vector3.Dot(p0p1, p0p1);
        float dotTest = Vector3.Dot(p0p1, p0test);

        return (0 < dotTest && dotTest < dotEdge); // Return False if 0 or dotEdge - if on an Endpoint
    } 
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
        // Use source code translated to C# from Moller, since mine Will Not Work
        Vector3 p1 = new Vector3(), p2 = new Vector3();
        bool intersected = TriTriOverlap.TriTriIntersect(V0.Position, V1.Position, V2.Position, 
            other.V0.Position, other.V1.Position, other.V2.Position,
            ref p1, ref p2);
        
        if (intersected)
        {
            o_intersection.P0 = p1;
            o_intersection.P1 = p2;
            Debug.DrawLine(p1, p2, Color.magenta, 10000);
        }
        return intersected;

        // -------------------------- My own from-scratch implementation of Moller's 1997 paper is inaccurate
        //                            detecting the actual intersection line, and attempts to use a method other
        //                            than the one he implements himself that seems to be inaccurate
        //                              - I attempt to (presumably incorrectly) calculate an O for the equation L = O + tD

        // -- Step 0: Copy data into operatorTris, since the order of Vertices may be manipulated.
        //            This is extremely dirty, but its better than coding in a set of equations for
        //            each case (4-6 different sets of equations, nested within several control blocks,
        //            is quite dirty)
        Triangle operandThis = new Triangle(), operandOther = new Triangle();
        operandThis.V0 = V0;
        operandThis.V1 = V1;
        operandThis.V2 = V2;
        operandOther.V0 = other.V0;
        operandOther.V1 = other.V1;
        operandOther.V2 = other.V2;

        // -- Step 1 - Compute Plane equation of Other Triangle (T2)
        Vector3 N2 = other.CalculateNormal();
        float d2 = Vector3.Dot(-N2, operandOther.V0.Position);
        Plane plane2 = new Plane(N2, d2); // Create a Unity Plane for access to functions like GetDistanceToPoint, etc.

        // Distances of T1 (this) vertices to T2 Plane (other)
        // Naming convention is d2v1i
        //  - d - Distance to
        //  - 2 - Triangle 2
        //  - v1 - Vertex of Triangle 1
        //  - i - Vertex # of Triangle 1
        float d2v10 = plane2.GetDistanceToPoint(operandThis.V0.Position);
        float d2v11 = plane2.GetDistanceToPoint(operandThis.V1.Position);
        float d2v12 = plane2.GetDistanceToPoint(operandThis.V2.Position);

        // -- Step 2 - Reject as trivial if all points of This (T1) are on the same side
        if (d2v10 != 0 && d2v11 != 0 && d2v12 != 0) // All distances are not 0 (no points on the plane)
        {
            float d2v10Sign = Mathf.Sign(d2v10);
            if (d2v10Sign == Mathf.Sign(d2v11) && d2v10Sign == Mathf.Sign(d2v12)) // If all share the same sign
            {
                //Debug.Log("All Vertices of Triangle 1 are on the same side of the Plane defined by Triangle 2");
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
                    Vertex tmpV = operandThis.V0; // Swap Vertices
                    operandThis.V0 = operandThis.V1;
                    operandThis.V1 = tmpV;

                    float tmpF = d2v10; // Swap distances calculated from the Vertices
                    d2v10 = d2v11;
                    d2v11 = tmpF;
                }
                else if (d2v11Sign == d2v10Sign) // first and second are on same side, so put third in second
                {
                    Vertex tmpV = operandThis.V2; // Swap Vertices
                    operandThis.V2 = operandThis.V1;
                    operandThis.V1 = tmpV;

                    float tmpF = d2v12; // Swap distances calculated from the Vertices
                    d2v12 = d2v11;
                    d2v11 = tmpF;
                }
            }
        }

        // -- Step 3 - Compute Plane equation of This Triangle (T1)
        Vector3 N1 = CalculateNormal();
        float d1 = Vector3.Dot(-N1, operandThis.V0.Position);
        Plane plane1 = new Plane(N1, d1);

        // Distances of T2 (other) vertices to T1 Plane (this)
        // Naming convention is d1v2i
        //  - d1 - Distance to Triangle 1
        //  - v2 - Vertex of Triangle 2
        //  - i - Vertex # of Triangle 2
        float d1v20 = plane1.GetDistanceToPoint(operandOther.V0.Position);
        float d1v21 = plane1.GetDistanceToPoint(operandOther.V1.Position);
        float d1v22 = plane1.GetDistanceToPoint(operandOther.V2.Position);

        // -- Step 4 - Reject as trivial if all points of Other (T2) are on the same side
        if (d1v20 != 0 && d1v21 != 0 && d1v22 != 0) // All distances are not 0 (no points on the plane)
        {
            float d1v20Sign = Mathf.Sign(d1v20);
            if (d1v20Sign == Mathf.Sign(d1v21) && d1v20Sign == Mathf.Sign(d1v22)) // If all share the same sign
            {
                //Debug.Log("All Vertices of Triangle 2 are on the same side of the Plane defined by Triangle 2");
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
                    Vertex tmpV = operandOther.V0; // Swap Vertices
                    operandOther.V0 = operandOther.V1;
                    operandOther.V1 = tmpV;

                    float tmpF = d1v20; // Swap distances calculated from the Vertices
                    d1v20 = d1v21;
                    d1v21 = tmpF;
                }
                else if (d1v21Sign == d1v20Sign) // first and second are on same side, so put third in second
                {
                    Vertex tmpV = operandOther.V2; // Swap Vertices
                    operandOther.V2 = operandOther.V1;
                    operandOther.V1 = tmpV;

                    float tmpF = d1v22; // Swap distances calculated from the Vertices
                    d1v22 = d1v21;
                    d1v21 = tmpF;
                }
            }
        }

        // -- Step 5 - Compute intersection line and project onto largest axis - Projection onto Largest axis is not done,
        // as I do not fully understand the notation described and how the simplification works
        Vector3 D = Vector3.Cross(N1, N2); // D in the equation L = O + tD - direction of the intersection line
        Vector3 O = new Vector3(0.0f, 0.0f, 0.0f); // Any point on line L, must be calculated to record Edge of Intersection
        // To calculate O, start with X = 0, and solve system of equations in planar general form
        // N1.Y*y + N1.Z*z + d1 = 0, N2.Y*y + N2.Z*z + d2 = 0
        // If D parallel to X, do Y, if parallel to Y too the Line is a problem
        // z = (N1.Y * y + d1) / -N1.Z
        // y = (N2.Z(N1.Y*y+d1)/-N1.Z + d2) / N2.Y
        // -N1.Z*N2.Y*y + N2.Z*N1.Y*y + N2.Z*d1 + -N1.Z*d2 = 0
        // N2.Z*N1.y*y - N1.Z*N2.Y*y = N1.Z*d2 - N2.Z*d1
        // y = (N1.Z*d2 - N2.Z*d1) / (N2.Z*N1.Y - N1.Z*N2.Y)
        O.y = ((N1.z * d2) - (N2.z * d1)) / ((N2.z * N1.y) - N1.z * N2.y);
        O.z = ((N1.y * O.y) + d1) / -N1.z;
        O = (Vector3.Dot(N2, operandThis.V0.Position)*operandThis.V1.Position - Vector3.Dot(N2, operandThis.V1.Position)*operandThis.V0.Position + d2*(operandThis.V1.Position - operandThis.V0.Position)) 
            / (Vector3.Dot(N2, operandThis.V0.Position) - Vector3.Dot(N2, operandThis.V1.Position));

        // Project Triangle 1 (this) points
        // Naming convention is pv1i
        //  - pv1 - Projection Vertex of Triangle 1
        //  - i - Vertex # of Triangle 1
        float pv10 = Vector3.Dot(D, operandThis.V0.Position - O);
        float pv11 = Vector3.Dot(D, operandThis.V1.Position - O);
        float pv12 = Vector3.Dot(D, operandThis.V2.Position - O);

        // Repeat for Triangle 2 (other) points, same convention
        float pv20 = Vector3.Dot(D, operandOther.V0.Position - O);
        float pv21 = Vector3.Dot(D, operandOther.V1.Position - O);
        float pv22 = Vector3.Dot(D, operandOther.V2.Position - O);

        // -- Step 6 - Compute intervals for each Triangle
        float t11 = (pv11 - pv10) * (d2v10 / (d2v10 - d2v11)) + pv10;
        float t12 = (pv11 - pv12) * (d2v12 / (d2v12 - d2v11)) + pv12;

        float t21 = (pv21 - pv20) * (d1v20 / (d1v20 - d1v21)) + pv20;
        float t22 = (pv21 - pv22) * (d1v22 / (d1v22 - d1v21)) + pv22;

        //Debug.Log($"\nT11: {t11}, T12: {t12}\nT21: {t21}, T22: {t22}");

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
            //Debug.Log("Distance between shortest and longest is greater than combined interval length, so they can't overlap");
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
        //foreach(float v in lengths)
        //{
        //    Debug.Log(v);
        //}
        //Debug.Log(D);
        o_intersection.P0 = O + D * lengths[1]; // Actual intersection is from index 1 to 2 - ISSUE: These numbers
        o_intersection.P1 = O + D * lengths[2]; // are local to the intersection, NOT World or even Local
        //Debug.DrawLine(O + D * lengths[0], O + D * lengths[3], Color.cyan, 1000);
        //Debug.DrawLine(o_intersection.P0, o_intersection.P1, Color.magenta, 1000);
        //Debug.Log("INTERSECTION FOUND");
        return true;
    }
}

// A Loop represents a Loop of connected Edges, which can either be Open or Closed. A Closed Loop
// represents a simple Polygon
class Loop
{
    private List<Edge> Edges;
    public List<Edge> GetEdges() { return Edges; }
    //private List<Vector3> Vertices;

    // Adds an Edge to the Polygon. The new Edge must start at the endpoint of the last Edge
    // in the Loop OR end at the start of the first Edge. If the Loop is Closed, no new Edges
    // may be added
    //  - This is very restrictive, but simplifies the construction process while still adding
    //    one Edge at a time
    public bool AddEdge(Edge e)
    {
        if (Edges == null) Edges = new List<Edge>();

        if (Edges.Count == 0) // Just Add as the first, but both Vertices are added
        {
            Edges.Add(e);
            //Vertices.Add(e.P0);
            //Vertices.Add(e.P1);
            return true;
        }
        else if (Closed()) // If the Loop is Closed, don't Add
        {
            return false;
        }

        if (e.P0 == Edges[Edges.Count-1].P1) // Edge added to the Front
        {
            Edges.Add(e);
            //Vertices.Add(e.P1);
            return true;
        }
        else if (e.P1 == Edges[0].P0) // Edge added to the Back
        {
            Edges.Insert(0, e);
            //Vertices.Insert(0, e.P0);
            return true;
        }
        
        return false;
    }

    // If the Polygon is not Closed, and already has a working Loop, add an Edge to the specified Point
    public bool AddPoint(Vector3 p)
    {
        if (Edges == null || Closed()) return false;

        Edge e = new Edge();
        e.P0 = Edges[Edges.Count-1].P1;
        e.P1 = p;
        Edges.Add(e);
        //Vertices.Add(e.P1);
        return true;   
    }

    // Returns true if the Polygon Loop is Closed (last Edge ends at first Edge's start)
    public bool Closed()
    {
        if (Edges == null || Edges.Count == 0) return false; // Can't be Closed if there is no Loop

        return Edges[Edges.Count - 1].P1 == Edges[0].P0;
    }

    // Returns true if the given Edge is in the Loop
    public bool ContainsEdge(Edge e)
    {
        return Edges.Contains(e);
    }

    // Retrieves the first point in the Loop. If there is no Loop, returns false and does not set the point
    public bool FirstPoint(ref Vector3 point)
    {
        if (Edges == null || Edges.Count == 0) return false;

        point = Edges[0].P0;
        return true;
    }

    // Retrieves the last pont in the Loop. If there is no Loop, returns false and does not set point
    public bool LastPoint(ref Vector3 point)
    {
        if (Edges == null || Edges.Count == 0) return false;

        point = Edges[Edges.Count - 1].P1;
        return true;
    }
}

// A TriangleMesh is a representation of raw Mesh using whole Triangles. This
// simplifies the management and code of the algorithm, converting Unity Mesh
// objects into StrikeMesh objects and back simplifies the number of parts that
// must be kept track of as the algorithm runs
struct TriangleMesh
{
    public List<Triangle> Triangles;

    // Source need not be normalized into non-indexed storage
    public void Create(Mesh source)
    {
        if (Triangles == null)
        {
            Triangles = new List<Triangle>();
        }
        else
        {
            Triangles.Clear();
        }

        for(int i = 0; i < source.triangles.Length; i += 3)
        {
            Vertex v0 = new Vertex(), v1 = new Vertex(), v2 = new Vertex();
            int t0 = source.triangles[i], t1 = source.triangles[i + 1], t2 = source.triangles[i + 2];
            
            v0.Position = source.vertices[t0];
            v1.Position = source.vertices[t1];
            v2.Position = source.vertices[t2];

            v0.Normal = source.normals[t0];
            v1.Normal = source.normals[t1];
            v2.Normal = source.normals[t2];

            if (source.tangents.Length > t0) v0.Tangent = source.tangents[t0];
            if (source.tangents.Length > t1) v1.Tangent = source.tangents[t1];
            if (source.tangents.Length > t2) v2.Tangent = source.tangents[t2];

            v0.UV = source.uv[t0];
            v1.UV = source.uv[t1];
            v2.UV = source.uv[t2];

            Triangle tri = new Triangle();
            tri.V0 = v0;
            tri.V1 = v1;
            tri.V2 = v2;
            Triangles.Add(tri);
        }
    }

    public Mesh Extract()
    {
        List<Vector3> vertices = new List<Vector3>();
        List<Vector3> normals = new List<Vector3>();
        List<Vector4> tangents = new List<Vector4>();
        List<Vector2> uvs = new List<Vector2>();
        List<int> indices = new List<int>();

        // Convert each Triangle into Mesh data
        foreach (Triangle tri in Triangles)
        {
            vertices.Add(tri.V0.Position);
            vertices.Add(tri.V1.Position);
            vertices.Add(tri.V2.Position);

            normals.Add(tri.V0.Normal);
            normals.Add(tri.V1.Normal);
            normals.Add(tri.V2.Normal);

            tangents.Add(tri.V0.Tangent);
            tangents.Add(tri.V1.Tangent);
            tangents.Add(tri.V2.Tangent);

            uvs.Add(tri.V0.UV);
            uvs.Add(tri.V1.UV);
            uvs.Add(tri.V2.UV);

            for(int i = 0; i < 3; i++)
            {
                indices.Add(indices.Count);
            }
        }

        // Create Unity Mesh and assign data
        Mesh mesh = new Mesh();
        mesh.vertices = vertices.ToArray();
        mesh.normals = normals.ToArray();
        mesh.tangents = tangents.ToArray();
        mesh.uv = uvs.ToArray();
        mesh.triangles = indices.ToArray();
        return mesh;
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
        if (coreModel != null && coreMesh != null)
        {
            Vector3 min = new Vector3(0, float.MaxValue, 0);
            for (int i = 0; i < coreMesh.vertices.Length; i++)
            {
                if (coreMesh.vertices[i].y < min.y)
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
        if (coreModel != null)
        {
            MeshFilter core = coreModel.GetComponent<MeshFilter>();
            if (core != null)
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
            for (int i = 0; i < 3; i++)
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
        foreach (int index in mesh.triangles)
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
        // ----------------------- TEST TRIANGLE INTERSECTION CODE ------------------------
        /*
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

        Edge ed = new Edge();
        t1.Intersect(t2, ref ed);
        Debug.Log($"Intersection Edge: {ed.P0}, {ed.P1}");
        return;
        */
        // --------------------- END TEST TRIANGLE INTERSECTION CODE ----------------------

        // - Step 1: Validate input
        if (coreMesh == null) { Debug.Log("No Core Mesh Found"); return; }
        NormalizeMesh(coreMesh);
        //NormalizeMesh(clipMesh);
        TriangleMesh core = new TriangleMesh(), clip = new TriangleMesh();
        core.Create(coreMesh);
        clip.Create(clipMesh);

        // Convert Triangle Mesh vertices to World Position
        for (int i = 0; i < core.Triangles.Count; i++)
        {
            Triangle coreTri = core.Triangles[i];
            coreTri.V0.Position = coreModel.transform.TransformPoint(coreTri.V0.Position);
            coreTri.V1.Position = coreModel.transform.TransformPoint(coreTri.V1.Position);
            coreTri.V2.Position = coreModel.transform.TransformPoint(coreTri.V2.Position);
            core.Triangles[i] = coreTri;
        }
        for (int i = 0; i < clip.Triangles.Count; i++)
        {
            Triangle clipTri = clip.Triangles[i];
            clipTri.V0.Position = transform.TransformPoint(clipTri.V0.Position);
            clipTri.V1.Position = transform.TransformPoint(clipTri.V1.Position);
            clipTri.V2.Position = transform.TransformPoint(clipTri.V2.Position);
            clip.Triangles[i] = clipTri;
        }

        // - Step 2: Search for pairs of intersecting Triangles
        // This step is primarily for identifying Triangles that may potentially
        // intersect each other using Octrees or some other type of spatial partitioning.
        // Since I am working with a relatively small number of Triangles, in the
        // interest of time this step is skipped, and an O(nh) loop is used


        // - Step 3: Intersect and re-triangulate

        // This dictionary records all Edges within a Triangle of intersections with other triangles.
        // From it Edge loops can be constructed (since any edge sharing a point MUST be in the same loop)
        //  - It is easier to record all Edges, then construct Loops, since that avoids the case of 2 Loops
        //    beginning separately and then being joined into one (ie: 4-5-6-7 and 8-9-10 and adding 7->8 to either of them)
        Dictionary<Triangle, List<Edge>> triIntersectionEdges = new Dictionary<Triangle, List<Edge>>();
        
        // For each triangle of the Core Mesh, record Edges of intersections with
        // the Clip Mesh
        foreach (Triangle coreTri in core.Triangles)
        {
            // Add the Triangle to the list if it is not already present (it never should be here)
            // Also add it's bounds as Edges for the fre-Triangulation
            if (!triIntersectionEdges.ContainsKey(coreTri))
            {
                triIntersectionEdges[coreTri] = new List<Edge>();
                //Edge v0v1 = new Edge(), v1v2 = new Edge(), v2v0 = new Edge();
                //v0v1.P0 = coreTri.V0.Position;
                //v0v1.P1 = coreTri.V1.Position;
                //v1v2.P0 = coreTri.V1.Position;
                //v1v2.P1 = coreTri.V2.Position;
                //v2v0.P0 = coreTri.V2.Position;
                //v2v0.P1 = coreTri.V0.Position;
                //triIntersectionEdges[coreTri].Add(v0v1);
                //triIntersectionEdges[coreTri].Add(v1v2);
                //triIntersectionEdges[coreTri].Add(v2v0);
            }

            // Check against every Clipping Triangle for intersections
            foreach (Triangle clipTri in clip.Triangles)
            {          
                // If this Clip Triangle is not recorded yet, add it (this should only happen on the first loop iteration)
                if (!triIntersectionEdges.ContainsKey(clipTri))
                {
                    triIntersectionEdges[clipTri] = new List<Edge>();
                    //Edge v0v1 = new Edge(), v1v2 = new Edge(), v2v0 = new Edge();
                    //v0v1.P0 = clipTri.V0.Position;
                    //v0v1.P1 = clipTri.V1.Position;
                    //v1v2.P0 = clipTri.V1.Position;
                    //v1v2.P1 = clipTri.V2.Position;
                    //v2v0.P0 = clipTri.V2.Position;
                    //v2v0.P1 = clipTri.V0.Position;
                    //triIntersectionEdges[clipTri].Add(v0v1);
                    //triIntersectionEdges[clipTri].Add(v1v2);
                    //triIntersectionEdges[clipTri].Add(v2v0);
                }

                // Intersect Core to Clip, adding the edge in both entries.
                //  - Potentially, the Edge entered into the Clip Tri should be reversed
                Edge intersection = new Edge();
                if (coreTri.Intersect(clipTri, ref intersection)) // Only add Edge if there was an intersection
                {
                    // Occasionally this is false, which is likely due to floating point error. Regardless,
                    // an Edge of a point to itself is considered irrelevant
                    if (intersection.P0 != intersection.P1)
                    {
                        // Intersection now holds the intersection Edge, but before it is added
                        // it must be checked for having a point on another Edge (it would then
                        // need to split that Edge in 2. 2 Edges should NEVER cross naturally)
                        //  - UPDATE: This shouldn't happen at all, since with well-formed meshes
                        //            the only crossing intersections should be on the edge of a Tri,
                        //            and the original Tri is not currently in the Edge List

                        // Check coreTri for Edge-splitting on Intersection Endpoints
                        triIntersectionEdges[coreTri].Add(intersection);
                        /*foreach (Edge e in triIntersectionEdges[coreTri])
                        {
                            // If either point is on the Edge, remove this Edge and replace it with 2 others
                            // OR statement ignores collinear Edges ------------------------------ TODO: POSSIBLE ERROR LOCATION
                            if (e.IsPointOnEdge(intersection.P0))
                            {
                                triIntersectionEdges[coreTri].Remove(e);
                                Edge e1 = new Edge(), e2 = new Edge();
                                e1.P0 = e.P0;
                                e1.P1 = intersection.P0;
                                e2.P0 = intersection.P0;
                                e2.P1 = e.P1;
                                triIntersectionEdges[coreTri].Add(e1);
                                triIntersectionEdges[coreTri].Add(e2);

                            }
                            else if (e.IsPointOnEdge(intersection.P1))
                            {
                                triIntersectionEdges[coreTri].Remove(e);
                                Edge e1 = new Edge(), e2 = new Edge();
                                e1.P0 = e.P0;
                                e1.P1 = intersection.P1;
                                e2.P0 = intersection.P1;
                                e2.P1 = e.P1;
                                triIntersectionEdges[coreTri].Add(e1);
                                triIntersectionEdges[coreTri].Add(e2);
                            }
                        }*/

                        // Check clipTri for Edge-splitting on Intersection Endpoints
                        triIntersectionEdges[clipTri].Add(intersection);
                        /*foreach (Edge e in triIntersectionEdges[clipTri])
                        {
                            // If either point is on the Edge, remove this Edge and replace it with 2 others
                            // OR statement ignores collinear Edges ------------------------------ TODO: POSSIBLE ERROR LOCATION
                            if (e.IsPointOnEdge(intersection.P0))
                            {
                                triIntersectionEdges[clipTri].Remove(e);
                                Edge e1 = new Edge(), e2 = new Edge();
                                e1.P0 = e.P0;
                                e1.P1 = intersection.P0;
                                e2.P0 = intersection.P0;
                                e2.P1 = e.P1;
                                triIntersectionEdges[clipTri].Add(e1);
                                triIntersectionEdges[clipTri].Add(e2);

                            }
                            else if (e.IsPointOnEdge(intersection.P1))
                            {
                                triIntersectionEdges[clipTri].Remove(e);
                                Edge e1 = new Edge(), e2 = new Edge();
                                e1.P0 = e.P0;
                                e1.P1 = intersection.P1;
                                e2.P0 = intersection.P1;
                                e2.P1 = e.P1;
                                triIntersectionEdges[clipTri].Add(e1);
                                triIntersectionEdges[clipTri].Add(e2);
                            }
                        }*/
                    }
                }
            }
        }

        // Triangulate takes a Triangle and a series of Edge Loops within that Triangle, splinters that Triangle into Polygons,
        // and Triangulates those Polygons. The output is a list of all new Triangles that will replace the origin Triangle
        Func<Triangle, List<Loop>, List<Triangle>> Triangulate = (Triangle a_origin, List<Loop> a_intersections) =>
        {
            Debug.Log("Triangulate Invoked");
            // Cut Loops off of the Triangle
            List<Loop> polys = new List<Loop>();
            Loop initialTri = new Loop();
            // Create initial Triangle polygon. Local braces added to avoid recording of variables
            {
                Edge e = new Edge();
                e.P0 = a_origin.V0.Position;
                e.P1 = a_origin.V1.Position;
                initialTri.AddEdge(e);
                initialTri.AddPoint(a_origin.V2.Position);
                initialTri.AddPoint(a_origin.V0.Position); // Close Polygon
            }
            polys.Add(initialTri);

            // Iteratively remove Loops
            foreach (Loop loop in a_intersections)
            {
                // If the Loop is closed it is located within one of the final Polygons, and severly complicates
                // Triangulation, but should not happen with well-formed meshes
                //  - In some error cases (likely floating point math), a Loop of 1 Edge with identical points
                //    will be found. This is culled before adding to the Triangle's intersection list
                if (loop.Closed())
                {
                    Debug.Log($"Closed Loop Intersection Found - Edge Count: {loop.GetEdges().Count}");
                    continue; // Closed Loops involve fractionation and are not supported
                }

                // Get Endpoints of the Loop
                Vector3 loopStart = new Vector3(), loopEnd = new Vector3();
                if (!loop.FirstPoint(ref loopStart) || !loop.LastPoint(ref loopEnd))
                {
                    Debug.Log("Loop does not exist?");
                    continue;
                }

                Loop cutPolygon = null; // Polygon that this Loop will cut

                // Find the Polygon in the Triangle to be split in 2 by this Loop
                foreach (Loop poly in polys)
                {
                    Edge[] loopEdges = new Edge[2]; // Index 0 is the starting Edge, Index 1 is the ending Edge

                    // Test to find Edges that hold the ends of the Loop
                    foreach(Edge polyEdge in poly.GetEdges())
                    {
                        // Check each endpoint separately, but they may be on the same Edge
                        if (polyEdge.IsPointOnEdge(loopStart))
                        {
                            loopEdges[0] = polyEdge;       
                        }
                        if (polyEdge.IsPointOnEdge(loopEnd))
                        {
                            loopEdges[1] = polyEdge;
                        }

                        // Cut out early if both Edges are found
                        if (loopEdges[0] != null && loopEdges[1] != null) break;
                    }

                    // If both endpoints are on this Polygon, this one will be split by the Loop
                    //  - With well-formed meshes, the only possible Polygon Edges that should have
                    //    both endpoints are Edges stemming from the original Triangle, so the found
                    //    polygon will have no duplicate on the opposite side of the Edge ------------- TODO: True?
                    if (loopEdges[0] != null && loopEdges[1] != null) 
                    {
                        cutPolygon = poly;
                        break; // The Loop cannot split more than one polygon, so do not keep searching
                    }
                }

                // If there was no Polygon found for this Loop to divide, there must have been an error
                if (cutPolygon == null)
                {
                    Debug.Log("Error: No Polygon within the Triangle found to split by this Edge Loop");
                    continue;
                }
            }

            List<Triangle> tris = new List<Triangle>();
            // Perform Triangulation

            return tris;
        };

        // Foreach Dictionary entry, form Edge Loops and Triangulate. If a point on an Edge does not
        // have a corresponding point in another Edge, then it is on one of the 3 original Tri edges
        foreach (var pair in triIntersectionEdges) // Key is the Triangle, Value is the list of intersecting Edges
        {
            if (pair.Value.Count == 0) continue; // If there are no intersections, no need to retriangulate

            // Accumulate all Intersection Edges into Loops (Open or Closed). These will then be added
            // to the source Triangle iteratively to splinter it for Re-Triangulation
            List<Loop> intersectionLoops = new List<Loop>();
            // While there are still Edges, Edge "Loops" can still be created, even if it is just one Edge
            while (pair.Value.Count != 0)
            {
                // Add first Edge
                Loop loop = new Loop();
                Edge initial = pair.Value[0];
                loop.AddEdge(initial);
                bool bEdgeAdded = true;

                // Iterate through other Edges, attempting to add to the Loop. This nested loop
                // is quite ugly, but it resolves a case where Edges in a loop are ignored because
                // they appear in different orderings in the List
                while (bEdgeAdded)
                {
                    bEdgeAdded = false;
                    foreach(Edge e in pair.Value)
                    {
                        if (loop.ContainsEdge(e)) continue; // Quick Exit

                        // Attempt to Add Edge
                        if(loop.AddEdge(e))
                        {
                            bEdgeAdded = true;
                            break; // Only add one Edge at a time
                        }
                    }
                }

                // No Edges were added, so remove the Loop from the list
                foreach (Edge e in loop.GetEdges())
                {
                    pair.Value.Remove(e);
                }
                intersectionLoops.Add(loop);
            }

            if (intersectionLoops.Count == 0) continue; // Do not retriangulate if no Loops are found
            List<Triangle> triangulation = Triangulate(pair.Key, intersectionLoops);
        }
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
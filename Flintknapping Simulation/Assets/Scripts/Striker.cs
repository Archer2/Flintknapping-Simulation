using System;
using System.Collections;
using System.Collections.Generic;
using TMPro.EditorUtilities;
using Unity.VisualScripting;
using UnityEngine;

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
        float height = 1; // Constant for now - TODO: Depend on Force
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
}

using System.Collections;
using System.Collections.Generic;
using UnityEngine;

/// <summary>
/// Basic Input script utilizing the New Input System
/// to capture basic Camera Movement input
/// </summary>
public class CameraInput : MonoBehaviour
{
    private CameraInputActions InputActions; // Actual Input object that captures data
    [SerializeField] private float Speed = 1;

    // On Start make sure the cursor is locked
    void Start()
    {
        Cursor.lockState = CursorLockMode.Locked;
    }

    // When this script awakens set the new InputActions
    void Awake()
    {
        InputActions = new CameraInputActions();
    }

    // Enable Input too
    void OnEnable()
    {
        InputActions.CameraMap.Enable();
    }

    // Disable Input too
    void OnDisable()
    {
        InputActions.CameraMap.Disable();
    }

    // Apparently the New Input System works better when accessed once per frame in Late Update
    //  - Particularly MouseDelta rotations
    private void LateUpdate()
    {
        // Translation
        transform.Translate(InputActions.CameraMap.Movement.ReadValue<Vector3>() * Time.deltaTime * Speed);

        // Rotation
        Vector2 mouseDelta = InputActions.CameraMap.Rotation.ReadValue<Vector2>() * Mathf.Rad2Deg * Time.deltaTime;
        transform.rotation *= Quaternion.AngleAxis(mouseDelta.y, Vector3.left);
        transform.rotation = Quaternion.Euler(transform.eulerAngles.x, transform.eulerAngles.y + mouseDelta.x, transform.eulerAngles.z);
    }
}

#!/usr/bin/env python3
import bpy
import json
import os
import math
import mathutils

def clear_scene():
    """Clear all objects in the scene"""
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete()
    
    # Clear all collections except the default one
    for collection in bpy.data.collections:
        bpy.data.collections.remove(collection)
    
    # Create collections for obstacles and robots
    obstacles_collection = bpy.data.collections.new("Obstacles")
    robots_collection = bpy.data.collections.new("Robots")
    bpy.context.scene.collection.children.link(obstacles_collection)
    bpy.context.scene.collection.children.link(robots_collection)

def create_material(name, color):
    """Create a new material with the given color"""
    material = bpy.data.materials.new(name=name)
    material.use_nodes = True
    material.node_tree.nodes["Principled BSDF"].inputs[0].default_value = color
    return material

def create_box(name, size, position, rotation, collection_name, material):
    """Create a box with the given parameters"""
    bpy.ops.mesh.primitive_cube_add(size=1, location=position)
    box = bpy.context.active_object
    box.name = name
    box.scale = (size[0], size[1], size[2])
    
    # Apply rotation as quaternion (w, x, y, z)
    box.rotation_mode = 'QUATERNION'
    box.rotation_quaternion = (rotation[0], rotation[1], rotation[2], rotation[3])
    
    # Link to collection
    bpy.data.collections[collection_name].objects.link(box)
    bpy.context.collection.objects.unlink(box)
    
    # Apply material
    if box.data.materials:
        box.data.materials[0] = material
    else:
        box.data.materials.append(material)
    
    return box

def create_sphere(name, radius, position, rotation, collection_name, material):
    """Create a sphere with the given parameters"""
    bpy.ops.mesh.primitive_uv_sphere_add(radius=radius, location=position, segments=32, ring_count=16)
    sphere = bpy.context.active_object
    sphere.name = name
    
    # Apply rotation as quaternion (w, x, y, z)
    sphere.rotation_mode = 'QUATERNION'
    sphere.rotation_quaternion = (rotation[0], rotation[1], rotation[2], rotation[3])
    
    # Link to collection
    bpy.data.collections[collection_name].objects.link(sphere)
    bpy.context.collection.objects.unlink(sphere)
    
    # Apply material
    if sphere.data.materials:
        sphere.data.materials[0] = material
    else:
        sphere.data.materials.append(material)
    
    return sphere

def create_capsule(name, radius, height, position, rotation, collection_name, material):
    """Create a capsule with the given parameters"""
    # Create a cylinder
    bpy.ops.mesh.primitive_cylinder_add(radius=radius, depth=height, location=position, vertices=32)
    capsule = bpy.context.active_object
    capsule.name = name
    
    # Add sphere at the top
    bpy.ops.mesh.primitive_uv_sphere_add(radius=radius, location=(position[0], position[1], position[2] + height/2))
    top_sphere = bpy.context.active_object
    
    # Add sphere at the bottom
    bpy.ops.mesh.primitive_uv_sphere_add(radius=radius, location=(position[0], position[1], position[2] - height/2))
    bottom_sphere = bpy.context.active_object
    
    # Join all objects
    bpy.ops.object.select_all(action='DESELECT')
    capsule.select_set(True)
    top_sphere.select_set(True)
    bottom_sphere.select_set(True)
    bpy.context.view_layer.objects.active = capsule
    bpy.ops.object.join()
    
    # Apply rotation as quaternion (w, x, y, z)
    capsule.rotation_mode = 'QUATERNION'
    capsule.rotation_quaternion = (rotation[0], rotation[1], rotation[2], rotation[3])
    
    # Link to collection
    bpy.data.collections[collection_name].objects.link(capsule)
    bpy.context.collection.objects.unlink(capsule)
    
    # Apply material
    if capsule.data.materials:
        capsule.data.materials[0] = material
    else:
        capsule.data.materials.append(material)
    
    return capsule

def create_cylinder(name, radius, height, position, rotation, collection_name, material):
    """Create a cylinder with the given parameters"""
    bpy.ops.mesh.primitive_cylinder_add(radius=radius, depth=height, location=position, vertices=32)
    cylinder = bpy.context.active_object
    cylinder.name = name
    
    # Apply rotation as quaternion (w, x, y, z)
    cylinder.rotation_mode = 'QUATERNION'
    cylinder.rotation_quaternion = (rotation[0], rotation[1], rotation[2], rotation[3])
    
    # Link to collection
    bpy.data.collections[collection_name].objects.link(cylinder)
    bpy.context.collection.objects.unlink(cylinder)
    
    # Apply material
    if cylinder.data.materials:
        cylinder.data.materials[0] = material
    else:
        cylinder.data.materials.append(material)
    
    return cylinder

def create_object_from_json(obj_data, collection_name, material):
    """Create an object based on the JSON data"""
    name = obj_data["name"]
    position = obj_data["position"]
    rotation = obj_data["rotation"]
    
    # Create different types of geometry based on the geometry type
    geometry = obj_data["geometry"]
    geometry_type = geometry["type"]
    
    if geometry_type == "box":
        size = geometry["size"]
        return create_box(name, size, position, rotation, collection_name, material)
    
    elif geometry_type == "sphere":
        radius = geometry["radius"]
        return create_sphere(name, radius, position, rotation, collection_name, material)
    
    elif geometry_type == "capsule":
        radius = geometry["radius"]
        height = geometry["height"]
        return create_capsule(name, radius, height, position, rotation, collection_name, material)
    
    elif geometry_type == "cylinder":
        radius = geometry["radius"]
        height = geometry["height"]
        return create_cylinder(name, radius, height, position, rotation, collection_name, material)
    
    elif geometry_type == "convex":
        # Create a placeholder for convex shapes (not fully implemented)
        bpy.ops.mesh.primitive_ico_sphere_add(radius=0.2, location=position)
        obj = bpy.context.active_object
        obj.name = name + "_convex_placeholder"
        
        # Link to collection
        bpy.data.collections[collection_name].objects.link(obj)
        bpy.context.collection.objects.unlink(obj)
        
        # Apply material
        if obj.data.materials:
            obj.data.materials[0] = material
        else:
            obj.data.materials.append(material)
        
        return obj
    
    return None

def setup_animation(obstacle_time_map, num_frames):
    """Setup animation for obstacles based on their time values"""
    # Set scene frame end
    bpy.context.scene.frame_end = num_frames
    
    # For each time frame, make objects from that time visible only in that frame
    for time, objects in obstacle_time_map.items():
        frame = int(time)
        
        # Set keyframes for all objects
        for obj in objects:
            # Make object invisible by default
            obj.hide_render = True
            obj.hide_viewport = True
            obj.keyframe_insert(data_path="hide_render", frame=0)
            obj.keyframe_insert(data_path="hide_viewport", frame=0)
            
            # Make object visible only for its time frame
            obj.hide_render = False
            obj.hide_viewport = False
            obj.keyframe_insert(data_path="hide_render", frame=frame)
            obj.keyframe_insert(data_path="hide_viewport", frame=frame)
            
            # Make object invisible again after its frame
            obj.hide_render = True
            obj.hide_viewport = True
            obj.keyframe_insert(data_path="hide_render", frame=frame+1)
            obj.keyframe_insert(data_path="hide_viewport", frame=frame+1)

def import_fcl_scene(json_path):
    """Import FCL scene from JSON file"""
    # Clear the scene first
    clear_scene()
    
    # Load JSON data
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    # Create materials
    obstacle_material = create_material("ObstacleMaterial", (0.8, 0.2, 0.2, 1.0))  # Red
    robot_material = create_material("RobotMaterial", (0.2, 0.2, 0.8, 1.0))        # Blue
    
    # Track obstacles by time for animation
    obstacle_time_map = {}
    max_time = 0
    
    # Create obstacles
    for obstacle_data in data["obstacles"]:
        time = obstacle_data["time"]
        max_time = max(max_time, int(time))
        
        # Create the object
        obj = create_object_from_json(obstacle_data, "Obstacles", obstacle_material)
        
        # Add to time map for animation
        if time not in obstacle_time_map:
            obstacle_time_map[time] = []
        obstacle_time_map[time].append(obj)
    
    # Create robots
    for robot_data in data["robots"]:
        create_object_from_json(robot_data, "Robots", robot_material)
    
    # Setup animation for obstacles
    setup_animation(obstacle_time_map, max_time + 1)
    
    print(f"Imported FCL scene with {len(data['obstacles'])} obstacles and {len(data['robots'])} robot links")
    print(f"Animation set up with {max_time + 1} frames")

# Main execution - can be run from Blender's Text Editor
def main():
    # Default path - can be changed as needed
    json_path = os.path.join(os.path.dirname(bpy.data.filepath), "fcl_export", "scene.json")
    
    # If run from command line with arguments
    import sys
    argv = sys.argv
    
    if "--" in argv:
        argv = argv[argv.index("--") + 1:]
        if len(argv) > 0:
            json_path = argv[0]
    
    import_fcl_scene(json_path)

if __name__ == "__main__":
    main() 
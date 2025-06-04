#!/usr/bin/env python3
import json
import os
import sys
import math
import argparse

def extract_joint_angles(json_file):
    """
    Extract joint angles from the JSON state export
    and convert to the format expected by robot_animator.py
    """
    # Load the JSON data
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Dictionary to store joint angles
    joint_angles = {}
    
    # Dictionary to map qIndex to joint name
    qindex_to_joint = {}
    
    # First pass: map qIndices to joint names
    for robot in data.get('robots', []):
        if 'joint' in robot and 'qIndex' in robot['joint']:
            joint_name = robot['name']
            qindex = robot['joint']['qIndex']
            
            # Store mapping of qIndex to joint name
            if qindex >= 0:  # Only store valid indices
                qindex_to_joint[qindex] = joint_name
    
    # Second pass: collect joint angles
    for robot in data.get('robots', []):
        if 'joint' in robot and 'angle' in robot['joint']:
            joint_name = robot['name']
            angle = robot['joint']['angle']
            qindex = robot['joint']['qIndex']
            
            # Store the angle using the joint name
            # First try numeric joint naming (joint1, joint2, etc.)
            if qindex >= 0 and qindex < 10:  # Assume single-digit joint indices
                joint_angles[f"joint{qindex+1}"] = angle
            
            # Also store with actual joint name for completeness
            if qindex >= 0:
                # Clean up the joint name - sometimes we need the suffix only
                clean_name = joint_name.split('_')[-1] if '_' in joint_name else joint_name
                joint_angles[clean_name] = angle
    
    return joint_angles

def convert_to_robot_animator_format(joint_angles, output_file=None, frame=1):
    """
    Convert the joint angles to the format expected by robot_animator.py
    
    Args:
        joint_angles: Dictionary of joint angles
        output_file: Path to output JSON file (if None, return the data)
        frame: Frame number for the animation
    
    Returns:
        JSON data structure if output_file is None, otherwise None
    """
    # Create the robot animator format
    animation_data = [
        {
            "frame": frame,
            "angles": joint_angles
        }
    ]
    
    # If output file is specified, write to file
    if output_file:
        with open(output_file, 'w') as f:
            json.dump(animation_data, f, indent=2)
        print(f"Converted joint angles written to {output_file}")
        return None
    
    # Otherwise return the data
    return animation_data

def batch_convert_directory(input_dir, output_dir=None):
    """
    Convert all JSON files in a directory
    
    Args:
        input_dir: Directory containing JSON state files
        output_dir: Directory to write converted files (if None, use input_dir)
    """
    if output_dir is None:
        output_dir = input_dir
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all JSON files
    json_files = [f for f in os.listdir(input_dir) if f.endswith('.json')]
    
    for i, json_file in enumerate(sorted(json_files)):
        input_path = os.path.join(input_dir, json_file)
        output_name = f"robot_animation_{i+1}.json"
        output_path = os.path.join(output_dir, output_name)
        
        try:
            # Extract joint angles
            joint_angles = extract_joint_angles(input_path)
            
            # Convert to robot animator format
            convert_to_robot_animator_format(joint_angles, output_path, frame=i+1)
            
            print(f"Converted {json_file} to {output_name}")
        except Exception as e:
            print(f"Error converting {json_file}: {e}")
    
    # Now create a combined animation file with all frames
    combine_animation_files(output_dir, os.path.join(output_dir, "combined_animation.json"))

def combine_animation_files(input_dir, output_file):
    """
    Combine multiple single-frame animation files into one multi-frame file
    
    Args:
        input_dir: Directory containing single-frame JSON files
        output_file: Path to output combined JSON file
    """
    # Find all JSON files
    json_files = [f for f in os.listdir(input_dir) if f.startswith('robot_animation_') and f.endswith('.json')]
    
    # Read and combine all frames
    combined_frames = []
    for json_file in sorted(json_files, key=lambda x: int(x.split('_')[2].split('.')[0])):
        input_path = os.path.join(input_dir, json_file)
        with open(input_path, 'r') as f:
            data = json.load(f)
            if data and len(data) > 0:
                combined_frames.extend(data)
    
    # Write combined frames to output file
    with open(output_file, 'w') as f:
        json.dump(combined_frames, f, indent=2)
    
    print(f"Combined animation written to {output_file}")

def convert_fcl_json_to_animation(fcl_json, output_file=None, frame_step=10):
    """
    Convert an FCL JSON export to a robot animation file
    
    Args:
        fcl_json: Path to FCL JSON file
        output_file: Path to output animation JSON file
        frame_step: Number of frames between each pose
    """
    try:
        # Load the FCL JSON data
        with open(fcl_json, 'r') as f:
            data = json.load(f)
        
        # Create dictionary of robots with their positions and rotations
        robots = {}
        for robot in data.get('robots', []):
            name = robot['name']
            robots[name] = {
                'position': robot.get('position', [0, 0, 0]),
                'rotation': robot.get('rotation', [1, 0, 0, 0])
            }
        
        # Extract joint angles if available
        joint_angles = extract_joint_angles(fcl_json)
        
        # Create animation data
        if joint_angles:
            animation_data = [
                {
                    "frame": 1,
                    "angles": joint_angles
                }
            ]
            
            # Write to output file
            if output_file:
                with open(output_file, 'w') as f:
                    json.dump(animation_data, f, indent=2)
                print(f"Converted animation written to {output_file}")
            
            return animation_data
        else:
            print("No joint angles found in the FCL JSON file")
            return None
    
    except Exception as e:
        print(f"Error converting FCL JSON to animation: {e}")
        return None

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Convert FCL JSON state export to robot animation')
    parser.add_argument('input', help='Input JSON file or directory')
    parser.add_argument('--output', '-o', help='Output file or directory')
    parser.add_argument('--batch', '-b', action='store_true', help='Process all JSON files in input directory')
    parser.add_argument('--frame', '-f', type=int, default=1, help='Frame number for single file conversion')
    
    args = parser.parse_args()
    
    if args.batch:
        # Batch process directory
        batch_convert_directory(args.input, args.output)
    else:
        # Process single file
        if os.path.isfile(args.input):
            joint_angles = extract_joint_angles(args.input)
            if joint_angles:
                convert_to_robot_animator_format(joint_angles, args.output, args.frame)
            else:
                print("No joint angles found in the JSON file")
        else:
            print(f"Input file not found: {args.input}")

if __name__ == "__main__":
    main() 
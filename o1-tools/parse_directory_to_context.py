#!/usr/bin/env -S conda run --live-stream -n scientific python

import os
import json
import sys
import xerox # pip install xerox

def build_structure_dict(root_dir):
    """
    Recursively walk through `root_dir` and build a dictionary
    capturing the structure and contents of `.R` files.
    """
    structure = []

    for entry in os.scandir(root_dir):
        if entry.is_dir():
            # If entry is a directory, recurse into it.
            structure.append({
                "type": "directory",
                "name": entry.name,
                "contents": build_structure_dict(entry.path)
            })
        else:
            # If entry is a file, check if it's an R file. If so, read content.
            if entry.name.endswith(".R"):
                with open(entry.path, "r", encoding="utf-8") as f:
                    file_contents = f.read()
                structure.append({
                    "type": "file",
                    "name": entry.name,
                    "relative_path": os.path.relpath(entry.path, root_dir),
                    "content": file_contents
                })
            else:
                # If you also want to keep track of non-R files, you can include them here.
                # Otherwise, you can omit them from the structure.
                structure.append({
                    "type": "file",
                    "name": entry.name,
                    "relative_path": os.path.relpath(entry.path, root_dir),
                    "content": None  # or "skipped" / remove key entirely
                })

    # Sort entries alphabetically by name for a tidy structure (optional)
    structure.sort(key=lambda x: x["name"])
    return structure

def parse_R_directory(directory_path):
    """
    Main driver function. Takes in a directory path, parses the structure, and
    returns a JSON string representing that structure (with all `.R` file contents).
    """
    # Build the directory/content structure.
    dir_structure = {
        "root_directory": os.path.abspath(directory_path),
        "structure": build_structure_dict(directory_path)
    }

    # Convert to JSON with indentation for readability.
    return json.dumps(dir_structure, indent=2)

if __name__ == "__main__":
    # If you pass the directory path as a CLI argument, e.g.:
    #   python parse_r_dir.py /path/to/shiny_app
    # or set a default here for testing:

    if len(sys.argv) > 1:
        root_dir = sys.argv[1]
    else:
        # Fallback: use current directory if none is provided
        root_dir = "."

    result_json = parse_R_directory(root_dir)
    print(result_json)

    xerox.copy(result_json)

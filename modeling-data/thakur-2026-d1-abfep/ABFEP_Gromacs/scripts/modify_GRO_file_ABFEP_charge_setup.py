#!/usr/bin/env python3
"""
Script to modify solvated.gro file for ligand charge compensation
"""

def modify_gro_file(input_file, output_file):
    """
    Modify a .gro file by:
    - Subtracting 2 from the atom count
    - Converting second-to-last TP3 to Na+ (remove H atoms)
    - Converting last TP3 to ION (O -> cl, H -> du)
    
    Preserves exact column formatting.
    """
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Line 0: Title
    # Line 1: Number of atoms
    # Lines 2 to -1: Atom lines
    # Last line: Box vectors
    
    title = lines[0]
    atom_count = int(lines[1].strip())
    box_line = lines[-1]
    atom_lines = lines[2:-1]
    
    # Decrease atom count by 2
    new_atom_count = atom_count - 2
    
    # Find all TP3 molecules by checking columns 5-10 (residue name)
    tp3_indices = []
    for i, line in enumerate(atom_lines):
        if len(line) >= 10:
            residue_name = line[5:10].strip()
            if residue_name == 'TP3':
                tp3_indices.append(i)
    
    if len(tp3_indices) < 6:
        print(f"Warning: Found only {len(tp3_indices)} TP3 atoms, need at least 6!")
        return
    
    # Get indices for the last two TP3 molecules (each has 3 atoms)
    last_tp3_indices = tp3_indices[-3:]  # Last 3 atoms
    second_last_tp3_indices = tp3_indices[-6:-3]  # Second-to-last 3 atoms
    
    # Process each atom line
    modified_lines = []
    lines_to_remove = set()
    
    # Mark H atoms from second-to-last TP3 for removal
    lines_to_remove.add(second_last_tp3_indices[1])  # H1
    lines_to_remove.add(second_last_tp3_indices[2])  # H2
    
    for i, line in enumerate(atom_lines):
        if i in lines_to_remove:
            continue  # Skip H atoms from second-to-last TP3
        
        # Second-to-last TP3: O atom -> Na+
        if i == second_last_tp3_indices[0]:
            # Only change columns 5-10 (residue) and 10-15 (atom name)
            modified_line = line[:5] + 'Na+  ' + '  Na+' + line[15:]
            modified_lines.append(modified_line)
        
        # Last TP3: O atom -> cl
        elif i == last_tp3_indices[0]:
            # Change residue to ION and atom to cl
            modified_line = line[:5] + 'ION  ' + '   cl' + line[15:]
            modified_lines.append(modified_line)
        
        # Last TP3: H1 -> du1
        elif i == last_tp3_indices[1]:
            modified_line = line[:5] + 'ION  ' + '  du1' + line[15:]
            modified_lines.append(modified_line)
        
        # Last TP3: H2 -> du2
        elif i == last_tp3_indices[2]:
            modified_line = line[:5] + 'ION  ' + '  du2' + line[15:]
            modified_lines.append(modified_line)
        
        else:
            modified_lines.append(line)
    
    # Write output file
    with open(output_file, 'w') as f:
        f.write(title)
        f.write(f"{new_atom_count}\n")
        for line in modified_lines:
            f.write(line)
        f.write(box_line)
    
    print(f"Modified .gro file written to {output_file}")
    print(f"Original atom count: {atom_count}")
    print(f"New atom count: {new_atom_count}")
    print(f"Found {len(tp3_indices)//3} TP3 molecules")


if __name__ == "__main__":
    import shutil
    
    # Always use solvated.gro
    input_file = "solvated.gro"
    temp_file = "solvated.gro.tmp"
    
    print(f"Modifying {input_file} in place...")
    
    # Modify and save to temporary file
    modify_gro_file(input_file, temp_file)
    
    # Replace original with modified version
    shutil.move(temp_file, input_file)
    
    print(f"Done! {input_file} has been modified.")

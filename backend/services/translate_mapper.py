from Bio.Data import CodonTable
import json
from pathlib import Path

def create_translations_from_biopython():
    """
    Create translations.json using Biopython's comprehensive genetic codes.
    """
    translations = []
    
    # Biopython has all NCBI genetic code tables
    genetic_code_names = {
        1: "Standard",
        2: "Vertebrate mitochondrial",
        3: "Yeast mitochondrial",
        4: "Mold, protozoan, and coelenterate mitochondrial + Mycoplasma/Spiroplasma",
        5: "Invertebrate mitochondrial",
        6: "Ciliate, dasycladacean and Hexamita nuclear",
        9: "Echinoderm and flatworm mitochondrial",
        10: "Euplotid nuclear",
        11: "Bacterial, archaeal and plant plastid",
        12: "Alternative yeast nuclear",
        13: "Ascidian mitochondrial",
        14: "Alternative flatworm mitochondrial",
        15: "Blepharisma nuclear",
        16: "Chlorophycean mitochondrial",
        21: "Trematode mitochondrial",
        22: "Scenedesmus obliquus mitochondrial",
        23: "Thraustochytrium mitochondrial",
        24: "Pterobranchia mitochondrial",
        25: "Candidate Division SR1 and Gracilibacteria",
    }
    
    for code_id, code_name in genetic_code_names.items():
        try:
            table = CodonTable.unambiguous_dna_by_id[code_id]
            
            # Forward table (codon -> amino acid)
            for codon, aa in table.forward_table.items():
                translations.append({
                    "whichGroup": code_name,
                    "Codon": codon,
                    "Translation": aa
                })
            
            # Stop codons
            for codon in table.stop_codons:
                translations.append({
                    "whichGroup": code_name,
                    "Codon": codon,
                    "Translation": "*"
                })
        
        except KeyError:
            print(f"Warning: Genetic code {code_id} not found")
            continue
    
    # Save to fileiokll.
    output_dir = Path("data")
    output_dir.mkdir(exist_ok=True)
    
    output_file = output_dir / "translations.json"
    with open(output_file, 'w') as f:
        json.dump(translations, f, indent=2)
    
    print(f"Created {output_file} with {len(translations)} entries")
    print(f"Genetic codes included: {len(genetic_code_names)}")

# if __name__ == "__main__":
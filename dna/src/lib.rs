use std::fmt;
use std::result;

#[derive(Debug, Clone)]
pub enum DNANucleotide {
    Adenine,
    Guanine,
    Cytosine,
    Thymine,
}

impl DNANucleotide {
    pub fn new(c: char) -> Self {
        match c {
            'A' => DNANucleotide::Adenine,
            'G' => DNANucleotide::Guanine,
            'C' => DNANucleotide::Cytosine,
            'T' => DNANucleotide::Thymine,
            _   => panic!("Invalid Base Character")
        }
    }

    pub fn get_partner(n: &DNANucleotide) -> DNANucleotide {
        match n {
            DNANucleotide::Adenine => DNANucleotide::Thymine,
            DNANucleotide::Thymine => DNANucleotide::Adenine,
            DNANucleotide::Cytosine => DNANucleotide::Guanine,
            DNANucleotide::Guanine => DNANucleotide::Cytosine,
        }
    }
}

impl fmt::Display for DNANucleotide {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> result::Result<(), fmt::Error> { 
        match self {
            DNANucleotide::Adenine => write!(f, "A"),
            DNANucleotide::Thymine => write!(f, "T"),
            DNANucleotide::Cytosine => write!(f, "C"),
            DNANucleotide::Guanine => write!(f, "G"),
        }
    }
}

#[derive(Debug, Clone)]
pub enum RNANucleotide {
    Adenine,
    Guanine,
    Cytosine,
    Uracil
}

impl RNANucleotide {
    pub fn new(c: char) -> Self {
        match c {
            'A' => RNANucleotide::Adenine,
            'G' => RNANucleotide::Guanine,
            'C' => RNANucleotide::Cytosine,
            'U' => RNANucleotide::Uracil,
            _   => panic!("Invalid Base Character")
        }
    }

    pub fn get_partner(n: &DNANucleotide) -> RNANucleotide {
        match n {
            DNANucleotide::Adenine => RNANucleotide::Uracil,
            DNANucleotide::Thymine => RNANucleotide::Adenine,
            DNANucleotide::Cytosine => RNANucleotide::Guanine,
            DNANucleotide::Guanine => RNANucleotide::Cytosine,
        }
    }
}

impl fmt::Display for RNANucleotide {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> result::Result<(), fmt::Error> { 
        match self {
            RNANucleotide::Adenine => write!(f, "A"),
            RNANucleotide::Uracil => write!(f, "U"),
            RNANucleotide::Cytosine => write!(f, "C"),
            RNANucleotide::Guanine => write!(f, "G"),
        }
    }
}

#[derive(Debug)]
pub enum AminoAcid {
    // Neutral, Non-polar amino acids
    Glycine, Alanine, Valine, Isoleucine, Tryptophan, 
    Proline, Methionine, Leucine, PhenylAlanine,

    // Acidic amino acids
    AsparticAcid, GlutamicAcid,

    // Basic amino acids
    Lysine, Arginine, Histidine,

    // Neutral Polar, amino acids
    Serine, Threonine, Tyrosine, Asparagine, Glutamine,
    Cysteine,

    STOP
}

impl fmt::Display for AminoAcid {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> result::Result<(), fmt::Error> { 
        match self {
            AminoAcid::Glycine => write!(f, "Gly"),
            AminoAcid::Alanine => write!(f, "Ala"),
            AminoAcid::Valine => write!(f, "Val"),
            AminoAcid::Isoleucine => write!(f, "IIe"),
            AminoAcid::Tryptophan => write!(f, "Try"),
            AminoAcid::Proline => write!(f, "Pro"),
            AminoAcid::Methionine => write!(f, "Met"),
            AminoAcid::Leucine => write!(f, "Leu"),
            AminoAcid::PhenylAlanine => write!(f, "Phe"),
            AminoAcid::Asparagine => write!(f, "Asn"),
            AminoAcid::AsparticAcid => write!(f, "Asp"),
            AminoAcid::GlutamicAcid => write!(f, "Glu"),
            AminoAcid::Lysine => write!(f, "Lys"),
            AminoAcid::Arginine => write!(f, "Arg"),
            AminoAcid::Serine => write!(f, "Ser"),
            AminoAcid::Histidine => write!(f, "His"),
            AminoAcid::Threonine => write!(f, "Thr"),
            AminoAcid::Tyrosine => write!(f, "Tyr"),
            AminoAcid::Glutamine => write!(f, "Glu"),
            AminoAcid::Cysteine => write!(f, "Cys"),
            AminoAcid::STOP => write!(f, "#")
        }
    }
}


/// A Gene is a sequence of DNA Nucleotides.
/// The nucleotides are arranged with
/// the 3' end at the first index and the 
/// 5' end at the last index
#[derive(Debug)]
pub struct DNASequence {
    sequence: Vec<DNANucleotide>
}

#[derive(Debug)]
pub struct RNASequence {
    sequence: Vec<RNANucleotide>
}

#[derive(Debug)]
pub struct AminoAcidSequence {
    sequence: Vec<AminoAcid>
}

#[derive(Debug)]
pub struct Codon {
    codon : [RNANucleotide; 3]
}

impl  fmt::Display for DNASequence {
     fn fmt(&self, f: &mut fmt::Formatter<'_>) -> result::Result<(), fmt::Error> {
         let mut string_seq = String::new();
         for n in &self.sequence {
            if let Some(c) = n.to_string().chars().next() {
                string_seq.push(c);
             }
         }
         write!(f, "{}", string_seq)
     }
}

impl  fmt::Display for RNASequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> result::Result<(), fmt::Error> {
        let mut string_seq = String::new();
        for n in &self.sequence {
           if let Some(c) = n.to_string().chars().next() {
               string_seq.push(c);
            }
        }
        write!(f, "{}", string_seq)
    }
}

impl fmt::Display for Codon {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> result::Result<(), fmt::Error> {
        write!(f, "({}, {}, {})", self.codon[0], self.codon[1], self.codon[2])
    }
}

impl DNASequence {
    pub fn new(s: String) -> Self {
        DNASequence {
            sequence: s.chars()
                      .map(|c| DNANucleotide::new(c))
                      .collect()
        }
    }

    pub fn get_complement(&self) -> Self {
        DNASequence {
            sequence: self.sequence
                          .clone().iter().rev()
                          .map(|n| DNANucleotide::get_partner(n))
                          .collect()
        }
    }

    pub fn transcribe(&self) -> RNASequence {
        RNASequence {
            sequence: self.sequence
                          .clone().iter().rev()
                          .map(|n| RNANucleotide::get_partner(n))
                          .collect()
        }
    }
}

impl RNASequence {
    pub fn new(s: String) -> Self {
        RNASequence {
            sequence: s.chars()
                      .map(|c| RNANucleotide::new(c))
                      .collect()
        }
    }

    pub fn translate(&self) -> AminoAcidSequence {
        unimplemented!()
    }
}

impl Codon {
    pub fn from_rna_sequence(rna: &RNASequence, start: usize, end: usize) -> Self {
        unimplemented!()
    }

    pub fn get_amino_acid(&self) -> AminoAcid {
        use RNANucleotide::*;
        use AminoAcid::*;
        match self.codon {
            [Uracil, Guanine, Guanine] => Tryptophan,
            [Uracil, Uracil, Uracil] | [Uracil, Uracil, Cytosine] => PhenylAlanine,
            [Uracil, Uracil, Adenine] | [Uracil, Uracil, Guanine] => Leucine,
            [Adenine, Uracil, Uracil] | [Adenine, Uracil, Cytosine] => Isoleucine,
            [Adenine, Uracil, Adenine] | [Adenine, Uracil, Guanine] => Methionine,
            [Adenine, Adenine, Uracil] | [Adenine, Adenine, Cytosine] => Asparagine,
            [Adenine, Adenine, Adenine] | [Adenine, Adenine, Guanine] => Lysine,
            [Adenine, Guanine, Uracil] | [Adenine, Guanine, Cytosine] => Serine,
            [Adenine, Guanine, Adenine] | [Adenine, Guanine, Guanine] => Arginine,
            [Cytosine, Adenine, Uracil] | [Cytosine, Adenine, Cytosine] => Histidine,
            [Cytosine, Adenine, Adenine] | [Cytosine, Adenine, Guanine] => Glycine,
            [Guanine, Adenine, Uracil] | [Guanine, Adenine, Cytosine] => AsparticAcid,
            [Uracil, Adenine, Uracil] | [Uracil, Adenine, Cytosine] => Tyrosine,
            [Uracil, Guanine, Uracil] | [Uracil, Guanine, Cytosine] => Cysteine,
            [Uracil, Adenine, Adenine] |
            [Uracil, Adenine, Guanine] | [Uracil, Guanine, Adenine] => STOP,
            [Uracil, Cytosine, Uracil] | [Uracil, Cytosine, Cytosine] |
            [Uracil, Cytosine, Adenine] | [Uracil, Cytosine, Guanine] => Serine,
            [Cytosine, Uracil, Uracil] | [Cytosine, Uracil, Cytosine] |
            [Cytosine, Uracil, Adenine] | [Cytosine, Uracil, Guanine] => Leucine,
            [Cytosine, Cytosine, Uracil] | [Cytosine, Cytosine, Cytosine] |
            [Cytosine, Cytosine, Adenine] | [Cytosine, Cytosine, Guanine] => Proline,
            [Adenine, Cytosine, Uracil] | [Adenine, Cytosine, Cytosine] |
            [Adenine, Cytosine, Adenine] | [Adenine, Cytosine, Guanine] => Threonine,
            [Guanine, Uracil, Uracil] | [Guanine, Uracil, Cytosine] |
            [Guanine, Uracil, Adenine] | [Guanine, Uracil, Guanine]  => Valine,
            [Guanine, Cytosine, Uracil] | [Guanine, Cytosine, Cytosine] |
            [Guanine, Cytosine, Adenine] | [Guanine, Cytosine, Guanine] => Alanine,
            _ => panic!(format!("{}: Not a valid Codon.", self))
        }
    }
}

#[cfg(test)]
mod test {
    use super::{Codon, DNANucleotide, DNASequence, RNANucleotide, RNASequence};
    #[test]
    fn test_gene() {
        let gene = DNASequence::new(String::from("ATGCCCAACGGCATGCT"));
        assert_eq!(gene.get_complement().to_string(), "AGCATGCCGTTGGGCAT");
    }
}




from Bio.SeqUtils.ProtParam import ProteinAnalysis

class Iterator:

    def __init__(self, dna_sequence):
        self.length = len(dna_sequence)
        self.sequence = dna_sequence
        self.it = 0
        self.data = {'A': 'T', 'T': 'A', 'C': 'G', 'G':'C'}

    def __iter__(self):
        return self

    def __next__(self):
        if self.length == self.it:
            raise StopIteration
        base = self.sequence[self.it]
        self.it += 1
        return self.data.get(base,'N')

class Sequence:
  
    def __init__(self, dna_sequence, accession_number):
       self.name = accession_number
       self.sequence = dna_sequence

    def content(self):
      dictionary={}
      sequence_length = len(self.sequence)
      content_A = self.sequence.count('A')
      dictionary['A'] = round(content_A/sequence_length*100,2)
      content_T = self.sequence.count('T')
      dictionary['T'] = round(content_T/sequence_length*100,2)
      content_C = self.sequence.count('C')
      dictionary['C'] = round(content_C/sequence_length*100,2)
      content_G = self.sequence.count('G')
      dictionary['G'] = round(content_G/sequence_length*100,2)
      return dictionary

    def length_of_sequence(self):
      return len(self.sequence)

class DNA(Sequence):

    def dna_to_cdna(self):
      cdna=''
      for letter in Iterator(self.sequence):
         cdna += letter
      return cdna

    def dna_to_rna(self):
        rna = self.sequence.replace('T','U')
        return rna

    def reverse(self):
        return self.sequence[::-1]

    def reverse_complement(self):
        cdna=''
        for letter in Iterator(self.sequence):
           cdna += letter
        return cdna[::-1]

    def translate(self):
        table_amino_acids = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
        }
        protein=''
        end = len(self.sequence) - (len(self.sequence) %3) - 1
        for i in range(0,end,3):
            codon = self.sequence[i:i+3]
            protein += table_amino_acids[codon]
        return protein 

    def mutate(self):
        import random
        mutation = random.randint(0,len(self.sequence)-1)
        mutation_in_dna = self.sequence[:mutation] + random.choice(['A', 'T', 'C', 'G']) + self.sequence [mutation+1:]
        return mutation_in_dna
      
class Protein:

      def __init__(self, translate_protein):
         self.protein = translate_protein

      def molar_mass(self):
          weight = ProteinAnalysis(self.protein)
          return round(weight.molecular_weight(),3)
      
      def aromaticity_aminoacids(self):
         aromatic = ProteinAnalysis(self.protein)
         return round(aromatic.aromaticity(),3)
       
      def pI (self):
         IP = ProteinAnalysis(self.protein)
         return round(IP.isoelectric_point(),2)

      def instability(self):
         instable = ProteinAnalysis(self.protein)
         return round(instable.instability_index(),2)

      def hydrophobicity (self):
         hydrophobic = ProteinAnalysis(self.protein)
         return round(hydrophobic.gravy(),2)


if __name__=='__main__':
    fh = open("DNA.txt")
    dna_from_file=''
    accession_number =''
    for line in fh:
        if line.startswith('>'):
            line_1 = line.lstrip('>').split()
            accession_number = line_1[0]
        else:
            dna_from_file += line
            dna_sequence = dna_from_file.replace('\n','')
    fh.close() 

    Sequence_object=Sequence(dna_sequence, accession_number)
    print (f'ACCESSION NUMBER: {Sequence_object.name}\n')
    print (f'DNA SEQUENCE: {Sequence_object.sequence} \n')
    dictionary = Sequence_object.content()
    print(f'LENGTH OF SEQUENCE: {Sequence_object.length_of_sequence()} bases\n')
    print (f'PERCENTAGE CONTENT OF BASES:')
    for key in dictionary:
         print(f'{key} {dictionary[key]} %')

    DNA_object = DNA(dna_sequence, accession_number)
    print (f'\nCOMPLEMENT: {DNA_object.dna_to_cdna()}')
    print (f'\nREVERSE: {DNA_object.reverse()}')
    print (f'\nCOMPLEMENT AND REVERSE: {DNA_object.reverse_complement()}')
    print (f'\nRNA SEQUENCE: {DNA_object.dna_to_rna()}')
    #print (f'\nDNA WITH ONE RANDOM MUTATION: {DNA_object.mutate()}')
    protein = DNA_object.translate()
    print (f'\nPROTEIN: {protein}')

    Protein_object = Protein(protein)
    print (f'\nPROTEIN MOLECULAR WEIGHT: {Protein_object.molar_mass()}')
    print (f'\nPROTEIN AROMATICITY: {Protein_object.aromaticity_aminoacids()}')
    print (f'\nPROTEIN ISOELECTRIC POINT: {Protein_object.pI()}')
    print (f'\nPROTEIN INSTABILITY: {Protein_object.instability()}')
    print (f'\nPROTEIN GRAVY: {Protein_object.hydrophobicity ()}')

    oh = open('PROTEIN.txt','w') 
    oh.write (f'PROTEIN: {protein}')
    oh.write (f'\nPROTEIN MOLECULAR WEIGHT: {Protein_object.molar_mass()}')
    oh.write (f'\nPROTEIN AROMATICITY: {Protein_object.aromaticity_aminoacids()}')
    oh.write (f'\nPROTEIN ISOELECTRIC POINT: {Protein_object.pI()}')
    oh.write (f'\nPROTEIN INSTABILITY: {Protein_object.instability()}')
    oh.write (f'\nPROTEIN GRAVY: {Protein_object.hydrophobicity ()}')
    oh.close()


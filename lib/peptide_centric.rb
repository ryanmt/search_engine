require 'mspire/fasta'


module PeptideCentricDB
  
  # returns a hash pointing from peptide to protein ID.
  # yields a protein sequence and expects an array of peptides
  #
  #    hash = PeptideCentricDB.create_peptide_to_protein_hash("yeast.fasta") do |aaseq| 
  #      digestor.digest(aaseq)
  #    end
  def self.create_peptide_to_protein_hash(fasta_file, &block)
    pep_to_prot_id = Hash.new {|h,k| h[k] = [] }
    Mspire::Fasta.open(fasta_file) do |fasta|
      fasta.each do |entry|
        prot_id = entry.header.split(' ',2).first
        peptides = block.call(entry.sequence)
        peptides.each do |pep|
          pep_to_prot_id[pep] << prot_id
        end
      end
    end
    pep_to_prot_id
  end
end

require 'yaml'
require 'lib/tools/digestor'
require 'peptide_centric'
module MS
  module SearchEngine
    class Searcher
      # Default Options Hash: Tolerances in ppm for now
      Options = {decoy: true, precursor_mass_tolerance: 10, ms2_mass_tolerance: 300, missed_cleavages: 2 } 
      def initialize(file)
        if File.extname(file) == ".fasta"
          @pep_hash = create_theoretical_spectra_file(file)
          @theoretical_spectra =  #TODO
        else 
          @theoretical_spectra = load_theoretical_spectra_from_file(file)
        end
      end
      def create_theoretical_spectra_file(fasta_file)
        outfile_name = File.join(File.dir_name(fasta_file),File.basename(fasta_file).sub('.fasta', 'changethis'))
        d = Digestor.new('trypsin')
        pep_centric_hash = PeptideCentricDB.create_peptide_to_protein_hash(fasta_file) do |aaseq|
          d.digest(aaseq, :missed_cleavages => Options[:missed_cleavages], :minimum_length => 6 )
        end
        # Write out the pepcentric db...
        pepcent_out = outfile_name.sub('changethis', ".pepcentric_#{Time.now.to_i}.yml") 
        putsv "Writing #{pepcent_out}"
        File.open(pepcent_out, 'w') do |out| 
          YAML.dump(pep_centric_hash, out) 
        end


        # Write out the theoretical spectra...
        theospec_out = outfile_name.sub('changethis', '.spectra.yml')
        putsv "Writing #{theospec_out}"
        File.open(theospec_out, 'w') do |out|
          YAML.dump(theoretical_spectra, out)
        end
      end
      def load_theoretical_spectra_from_file(file)
        #Load the YAML file
      end
    end
  end
end

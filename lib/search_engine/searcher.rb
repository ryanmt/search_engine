require 'json'
require 'tools/digestor'
require 'peptide_centric'
require 'tools/fragmenter'

module MS
  module SearchEngine
    class Searcher
      def initialize(file)
        if File.extname(file) == ".fasta"
          putsv "Creating: #create_theoretical_spectra"
          @theoretical_spectra = create_theoretical_spectra(file)
        else 
          putsv "Loading: #load_theoretical_spectra_from_file"
          @theoretical_spectra = load_theoretical_spectra_from_file(file)
        end
        self
      end
      def create_theoretical_spectra(fasta_file)
        outfile_name = File.join(File.dirname(fasta_file),File.basename(fasta_file).sub('.fasta', 'changethis'))
        d = Digestor.new('trypsin')
        pep_centric_hash = PeptideCentricDB.create_peptide_to_protein_hash(fasta_file) do |aaseq|
          d.digest(aaseq, :missed_cleavages => Options[:missed_cleavages], :minimum_length => 6 )
        end
        if MS::SearchEngine::Options[:write_pepdb]
          # Write out the pepcentric db...
          pepcent_out = outfile_name.sub('changethis', ".pepcentric_#{Time.now.to_i}.yml") 
          putsv "Writing #{pepcent_out}"
          File.open(pepcent_out, 'w') do |out| 
            out.print pep_centric_hash.to_json
          end
        end # Options
        # Generate theoretical spectra
        fragger = MS::Fragmenter.new
        theoretical_spectra = pep_centric_hash.keys.map do |aaseq|
          mzs = fragger.fragment(aaseq)
          MS::DataStructs::TheoreticalSpectrum.new(mzs, Array.new(mzs.size, 100.0), aaseq)
        end
        if MS::SearchEngine::Options[:write_spectradb] 
          # Write out the theoretical spectra...
          theospec_out = outfile_name.sub('changethis', ".spectra_#{Time.now.to_i.to_s[3,7]}.yml")
          putsv "Writing #{theospec_out}"
          File.open(theospec_out, 'w') do |out|
            out.print theoretical_spectra.to_json
          end
        end # Option
        theoretical_spectra
      end
      def load_theoretical_spectra_from_file(file)
        #Load the YAML file
        File.open(file, 'r') do |io|
          JSON.parse(io)
        end
        putsv "Loaded the file"
      end
      def search_mzml(mzml_file)
        #Open mzML
        #Each spectrum:
          # Normalize??? XCORR is robust to normalizations...
          # Bin
          # Compare to fragmentation scans within selected mass window (ONLY CALL THIS ONCE, as the experimental scan will be preprocessed then compared to each theoretical scan)
      end

    end
  end
end

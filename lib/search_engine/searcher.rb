require 'tools/digestor'
require 'peptide_centric'
require 'tools/fragmenter'
require 'bin'

require 'mspire/mzml'

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
          MS::DataStructs::TheoreticalSpectrum.new(mzs, Array.new(mzs.size, 100.0), MS::DataStructs::Peptide.new(aaseq))
        end
        if MS::SearchEngine::Options[:write_spectradb] 
          # Write out the theoretical spectra...
          theospec_out = outfile_name.sub('changethis', ".spectra_#{Time.now.to_i.to_s[3,7]}.yml")
          putsv "Writing #{theospec_out}"
          File.open(theospec_out, 'w') do |out|
            Marshal.dump(theoretical_spectra, out)
          end
        end # Option
        theoretical_spectra
      end
      def load_theoretical_spectra_from_file(file)
        #Load the YAML file
        tmp = nil
        File.open(file, 'r') do |io|
          tmp = Marshal::load(io)
        end
        putsv "Loaded the file"
        tmp
      end
      def search_mzml(mzml_file)
        filter_theoretical_spectra_by_precursor_mass
        experimental_spectra = []
        Mspire::Mzml.open(mzml_file) do |mzml|
          mzml.each do |spec|
            next unless spec.ms_level == 2 
            experimental_spectra << MS::DataStructs::ExperimentalSpectrum.new(spec.mzs, spec.intensities.map{|v| Math.sqrt(v) }, spec.precursor_mz, spec.precursor_charge)
          end
        end
        #Each spectrum:
        experimental_spectra.each do |spectrum|
          spectrum.bins = MS::Bin.create_bins_by_ppm(spectrum.mzs.min,spectrum.mzs.max  , MS::SearchEngine::Options[:ms2_mass_tolerance])
          MS::Bin.bin(spectrum.bins, spectrum.mzs_and_intensities, &:first)
          spectrum.bins.each do |bin|
            bin.data= bin.data.map(&:last).inject(:+)
          end
          spectrum.bins = normalize_by_section(spectrum.bins)
          range = PPM.mass_to_window(spectrum.precursor_neutral_mass, MS::SearchEngine::Options[:precursor_mass_tolerance])
          search_space = @theoretical_spectra.select {|a| range.include? a.precursor_mass}
          
          
# NOW I need to bin the search_space arrays to the same bins and then process them
          #search_space.map{|a| a.bins.map(&:data)}
          p XCorr.compare_to(search_space, spectrum.bins.map(&:data))
        end
       # TODO Finish out the search 
      end
      def normalize_by_section(bins)
        end_point = bins.map(&:end).max
        section_top = bins.map(&:begin).min
        incrementer = MS::SearchEngine::Options[:normalization_window_width]
        sections = []
        bins.each do |bin|
          bin.data = bin.data.nil? ? 0 : bin.data
          if bin.end > section_top
            break if section_top == end_point
            section_top += incrementer
            sections << []
          end
          sections.last << bin
        end
        sections.each do |section|
          max = section.map(&:data).max
          section.map {|bin| bin.data = bin.data / max * 100 }
        end
        sections.flatten
      end
      def filter_theoretical_spectra_by_precursor_mass
        @theoretical_spectra.sort_by! {|spectrum| spectrum.precursor_neutral_mass }
      end
    end
  end
end

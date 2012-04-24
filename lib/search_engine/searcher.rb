require 'tools/digestor'
require 'peptide_centric'
require 'tools/fragmenter'
require 'bin'
require 'tools/decoy'

require 'mspire/mzml'

module MS
  module SearchEngine
    class Searcher
      def initialize(file)
        if File.extname(file) == ".fasta"
          decoy_file = create_decoy_file(file) 
          putsv "Creating: #create_theoretical_spectra"
          @theoretical_spectra = create_theoretical_spectra(file)
          @decoy_theoretical_spectra = create_theoretical_spectra(decoy_file)
        else 
          putsv "Loading: #load_theoretical_spectra_from_file"
          @theoretical_spectra = load_theoretical_spectra_from_file(file)
#          @decoy_theoretical_spectra = load_theoretical_spectra(decoy_file)
        end
        @experimental_matched_spectra = []
        @decoy_matched_spectra = []
        self
      end
      def search(file)
        @experimental_matched_spectra = search_mzml(file, @theoretical_spectra)
        @decoy_matched_spectra = search_mzml(file, @decoy_theoretical_spectra)
        qvalues(@experimental_matched_spectra, @decoy_matched_spectra)
        output( @experimental_matched_spectra)
      end

      def qvalues(experiment, decoy)
        [experiment, decoy].map {|spectra| spectra.sort_by! {|spectrum| spectrum.match.nil? ? 0 : spectrum.match.xcorr} }
        decoy_xcorrs = decoy.map{|i| i.match.nil? ? 0 : i.match.xcorr}
        real_xcorrs = experiment.map{|i| i.match.nil? ? 0 : i.match.xcorr}
        decoy_size = decoy.size
        real_size = experiment.size
        experiment.map do |spectra|
          decoy_count = decoy_xcorrs.map{|i| i >= spectra.match.xcorr ? true : nil }.compact.size
          real_count = real_xcorrs.map{|i| i >= spectra.match.xcorr ? true : nil }.compact.size
          spectra.match.qvalue = decoy_count/(real_count+decoy_count)
        end
      end

      def output(spectra)
        puts %w{ID qvalue Xcorr ppm Peptide Proteins}.join("\t\t") 
        puts "spectra.first.match.theoretical_spectrum.peptide: #{spectra.first.match.theoretical_spectrum.peptide}"
        spectra.each_with_index do |spectrum, i|
          puts ["scan=#{i+1}",spectrum.match.qvalue, spectrum.match.xcorr, spectrum.match.ppm, spectrum.match.theoretical_spectrum.peptide, @pep_centric_hash[spectrum.match.theoretical_spectrum.peptide]].join("\t\t")
        end
      end
      def create_decoy_file(file)
        MS::Decoy.new.generate(file)
      end
      def create_theoretical_spectra(fasta_file)
        outfile_name = File.join(File.dirname(fasta_file),File.basename(fasta_file).sub('.fasta', 'changethis'))
        d = Digestor.new('trypsin')
        @pep_centric_hash = PeptideCentricDB.create_peptide_to_protein_hash(fasta_file) do |aaseq|
          d.digest(aaseq, :missed_cleavages => Options[:missed_cleavages], :minimum_length => 6 )
        end
        if MS::SearchEngine::Options[:write_pepdb]
          # Write out the pepcentric db...
          pepcent_out = outfile_name.sub('changethis', ".pepcentric_#{Time.now.to_i}.yml") 
          putsv "Writing #{pepcent_out}"
          File.open(pepcent_out, 'w') do |out| 
            out.print @pep_centric_hash.to_yaml
          end
        end # Options
        # Generate theoretical spectra
        fragger = MS::Fragmenter.new
        theoretical_spectra = @pep_centric_hash.keys.map do |aaseq|
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
      def search_mzml(mzml_file, theoretical_spectra_set)
        matched_spectra = []
        filter_theoretical_spectra_by_precursor_mass(theoretical_spectra_set)
        experimental_spectra = []
        Mspire::Mzml.open(mzml_file) do |mzml|
          mzml.each do |spec|
            next unless spec.ms_level == 2 
            experimental_spectra << MS::DataStructs::ExperimentalSpectrum.new(spec.mzs, spec.intensities.map{|v| Math.sqrt(v) }, spec.precursor_mz, spec.precursor_charge)
          end
        end
        #Each spectrum:
        experimental_spectra.each do |spectrum|
          bins = MS::Bin.create_bins_by_ppm(spectrum.mzs.min,spectrum.mzs.max  , MS::SearchEngine::Options[:ms2_mass_tolerance])
          spectrum.bins = Marshal.load(Marshal.dump(bins))
          MS::Bin.bin(spectrum.bins, spectrum.mzs_and_intensities, &:first)
          spectrum.bins.each do |bin|
            bin.data= bin.data.map(&:last).inject(:+)
          end
          spectrum.bins = normalize_by_section(spectrum.bins)
          range = PPM.mass_to_window(spectrum.precursor_neutral_mass, MS::SearchEngine::Options[:precursor_mass_tolerance])
          search_space = theoretical_spectra_set.select {|a| range.include? a.precursor_mass}
          
# NOW I need to bin the search_space arrays to the same bins and then process them
          search_space.map do |theo_spectrum|
            theo_spectrum.bins = Marshal.load(Marshal.dump(bins))
            MS::Bin.bin(theo_spectrum.bins, theo_spectrum.mzs_and_intensities.sort_by(&:first), &:first)
            theo_spectrum.bins.each do |bin| 
              bin.data = bin.data.empty? ? 0 : bin.data.map(&:last).inject(:+)
                #puts "bin.data.each: #{a}"
            end
          end
          #search_space.map{|a| a.bins.map(&:data)}
          xcorrs = XCorr.compare_to(search_space.map{|a| a.bins.map(&:data)}, spectrum.bins.map(&:data))
          
          best_xcorr_and_index = xcorrs.each_with_index.max
          spectrum.match = MS::DataStructs::Match.new(best_xcorr_and_index.first, search_space[best_xcorr_and_index.last], PPM.calculate(search_space[best_xcorr_and_index.last].precursor_mass, spectrum.precursor_neutral_mass)) unless search_space.empty?
          matched_spectra << spectrum
        end
        matched_spectra
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
      def filter_theoretical_spectra_by_precursor_mass(theoretical_spectra_set)
        theoretical_spectra_set.sort_by! {|spectrum| spectrum.precursor_neutral_mass }
      end
    end
  end
end

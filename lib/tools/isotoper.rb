
require 'tools/fragmenter'
require 'rserve/simpler'
require 'mspire/isotope/aa'
r = Rserve::Simpler.new
module MS
  class Isotoper
    Neutron = 1.00866491600
    TableEntry = Struct.new(:ion, :seq, :mass, :charge, :composition_arr)
    Defaults = { charge_state: 1, graph: false, normalize: :first_peak}
    def initialize(options = {})
      @opts = Defaults.merge(options)
    end
# Control FXN for this group... 
    def generate_spectra(sequence)
      @sequence = sequence
      #TODO Generate a ms1 
      ms1_entries = normalize(isotopic_distribution(sequence))
      #TODO Generate a ms2
      #ms2_entries = Fragmenter.new(fragmenter_opts).fragment(sequence, fragment_opts)
      #spectrum_to_mgf(ms1_entries)
      ms1_entries
    end
#OUTPUT FXNS
# This fxn will take the list of TableEntry values and send them to mgf.
    def spectrum_to_mgf(data)
      output_arr = []
      output_arr << %Q{COM=Project: In-silico Isotoper\nBEGIN IONS\nCHARGE=#{@charge_state}+\nTITLE=Label: Sequence is #{@sequence} and options were #{@opts}}
      data.each do |entry|
        #	TableEntry = Struct.new(:ion, :seq, :mass, :charge)
        output_arr << "#{"%.5f" % entry.first }\t#{"%.5f" % entry.last}"
      end
      output_arr << "END IONS"
      File.open("#{@sequence}.mgf", 'w') {|o| o.print output_arr.join("\n") }
      output_arr.join("\n")
    end
# This fxn will take a list of masses and graph them.
    def graph(data = nil)
      robj = Rserve::Simpler.new
      hash = {}
      hash["mass"] = data.map(&:first)
      hash["intensity"] = data.map(&:last)
      robj.converse( masses: hash.to_dataframe) do 
        "attach(masses)"
      end
      #robj.converse( data: Rserve::DataFrame.from_structs(list))
      robj.converse "setwd('#{Dir.pwd}')"
      output_file_name = "#{@sequence}_spectra.png"
      robj.converse do 
        %Q{png(file='#{output_file_name}')
          plot(masses$mass, masses$intensity, type='h')
          dev.off()
        }
      end	
      output_file_name
    end
    def graph_crap(list)
      robj = Rserve::Simpler.new
      hash = {}
      hash["intensity"] = list
      hash["mass"] = list.each_index.map {|i| i} # Hacky mass definition
      robj.converse( masses: hash.to_dataframe) do 
        "attach(masses)"
      end
      #robj.converse( data: Rserve::DataFrame.from_structs(list))
      robj.converse "setwd('#{Dir.pwd}')"
      output_file_name = "#{@sequence}_spectra.png"
      robj.converse do 
        %Q{png(file='#{output_file_name}')
          plot(masses$mass, masses$intensity, type='h')
          dev.off()
      }
      end	
      output_file_name
    end        
    def monoisotopic_mass 
      @m.monoisotopic
    end
    private
# Private tool functions
    def normalize(fft_resp)
      case @opts[:normalize]
        when :none
          fft_resp
        when :first_peak
          first_intensity = fft_resp.first[1]
          fft_resp.map {|arr| arr[1] = arr[1]/first_intensity*100}
        when :sum
          sum = fft_resp.map(&:last).inject(:+)
          fft_resp.map {|arr| arr[1] = arr[1]/sum*100}
      end
      fft_resp
    end
    def self.precursor_mass(sequence)
      @atom_counts = Mspire::Isotope::AA::ATOM_COUNTS
      resp = sequence.each_char.map {|aa| @atom_counts[aa] }
      pep_counts = Hash.new {|h,k| h[k] = 0 }
      resp.map do |h| 
        h.keys.each {|k| pep_counts[k] += h[k] }
      end
      pep_counts[:h] +=2
      pep_counts[:o] +=1
      require_relative 'isotoper/molecule'
      @m = Molecule.new(pep_counts, charge_state: 0, name: sequence)
      @m.mass(@m.counts_to_formula)
    end
# This function would generate the ms1 spectra masses for a given peptide sequence.
    def isotopic_distribution(sequence)
      @atom_counts = Mspire::Isotope::AA::ATOM_COUNTS
      resp = sequence.each_char.map do |aa|
        @atom_counts[aa]
      end
      pep_counts = Hash.new {|h,k| h[k] = 0 }
      resp.each do |h|
        h.keys.each {|k| pep_counts[k] += h[k]}
      end
      pep_counts[:h] += 2
      pep_counts[:o] += 1
      require_relative 'isotoper/molecule'
      @m = Molecule.new(pep_counts, charge_state: @opts[:charge_state], name: sequence)
      @m.isotopic_distribution
    end
  end
end

if __FILE__ == $0
=begin 
  i = Isotoper.new
# TEMPORARILY, #generate_spectra will return a Molecule
  m = i.generate_spectra "PEPTIDE"
  dist = m.isotopic_distribution
  p i.graph_crap(dist)
=end
  # TODO normalization options
# OPT PARSER
  options = {:verbose => false, :charge_state => 1, :normalize => :first_peak, graph: false}; require 'optparse'
  parser = OptionParser.new do |opts|
    opts.banner  = "Usage #{File.basename(__FILE__)} sequence [options]"
    opts.separator "Output: [String] sequence.mgf [sequence.svg]"
    opts.on('-p', '--protons N', Integer, "Number of protons to add." ) do |p| #  Otherwise, I'll predict them based upon peptide sequence") do |p|
      options[:charge_state] = p ? p : 1
    end 
    Norm_opts = [:none, :first_peak, :sum]
    opts.on('-n', '--normalize [TYPE]', Norm_opts, "Set normalization style: #{Norm_opts.join(', ')}") do |n|
      options[:normalize] = n
    end
    opts.on('-g', '--graph', "Graph the output as well as outputting the mgf") do |g|
      options[:graph] = g 
    end
    opts.on_tail('-h', '--help', "Show this message") do 
      puts opts
      exit
    end
  end.parse! # OptionParser
# RUN
  i = Isotoper.new(options)
  ARGV.each do |pepseq|
    ms1_data = i.generate_spectra(pepseq)
    if options[:graph]
      i.graph(ms1_data)
    end
  end # ARGV.each
end

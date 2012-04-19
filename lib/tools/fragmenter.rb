require_relative 'fragmenter/masses'
require_relative 'charge_calculator'
module MS
  class Fragmenter
    include MS

    attr_accessor :list
    TableEntry = Struct.new(:ion, :seq, :mass, :charge, :composition_arr)
    Ion_Defaults = {:b => true, :y => true}
    Defaults = {:charge_states => true, :avg => false}
    def initialize(opts = {}, ion_opts = {})
      set_options(opts, ion_opts)
      self
    end
    def set_options(opts, ion_opts)
      #@opts = Default_fragments.merge(opts)
      opts = Defaults.merge(opts)
      ion_opts = Ion_Defaults.merge(ion_opts)
      @n_term_search_ion_types = []
      @c_term_search_ion_types = []
      @max_charge = 1 unless opts[:charge_states]
      #puts "options :charge_states = #{opts[:charge_states]}"
      ion_opts.each do |key, v|
        if v
          case key 
            when :b
              @n_term_search_ion_types.push(:b, :b_star, :b_not)
            when :a
              @n_term_search_ion_types.push(:a, :a_star, :a_not)
            when :c
              @n_term_search_ion_types << :c
            when :x
              @c_term_search_ion_types << :x
            when :y
              @c_term_search_ion_types.push(:y, :y_star, :y_not)
            when :z
              @c_term_search_ion_types << :z
          end
        end
      end
      @mass_list = opts[:avg] ? MS::AvgResidueMasses : MS::MonoResidueMasses
      #putsv "@mass_list: #{@mass_list}"
    end #set_options
      
    def calculate_fragments(sequence)
      arr = sequence.upcase.split('')
      out = [[],[]]
      (0..arr.size-2).each do |i|
        out[0] << arr[0..i].join
        out[1] <<  arr[(i+1)..-1].join
      end
      out
    end
# This fxn exists to provide the API consistent with John's request for the 689R class.
# Options may include a list of fragment classes as symbols (i.e. :b, :y)
    def fragment(pep_seq, options={}) # TODO handle an intensity option to handle normalization and scaling...?
      set_options(options) unless options.empty?
      generate_fragment_masses(pep_seq)
      @list.map(&:mass)
    end
    def generate_fragment_masses(sequence) # Returns the TableEntry object which should be easy to use for table generation
      @sequence = sequence
      @max_charge ||= MS::ChargeCalculator.charge_at_pH(MS::ChargeCalculator.identify_potential_charges(sequence), 2).ceil
      ### Calculate the base ion masses	
      n_terms, c_terms = calculate_fragments(sequence)
      n_terms.map! do |seq|
        mass = MS::NTerm
        seq.chars.map(&:to_sym).each do |residue|
          mass += @mass_list[residue]
        end
        [seq, mass]
      end
      c_terms.map! do |seq|
        mass = MS::CTerm
        seq.chars.map(&:to_sym).each do |residue|
          mass += @mass_list[residue]
        end
        [seq, mass]
      end
### Tablify and generate a comprehensive list of ions
      list = []
      send_to_list = lambda do |fragment_arr, iontypes_arr|
        fragment_arr.each do |n_terms|
          seq = n_terms.first
          mass = n_terms.last
          iontypes_arr.each do |iontype|
            (1..@max_charge).each do |charge|
              charge_legend = '+'*charge
              list << TableEntry.new("#{iontype}(#{seq.size})#{charge_legend}".to_sym, seq, charge_state(mass + IonTypeMassDelta[iontype], charge), charge)
            end # 1..max_charge
          end	# iontypes_arr
        end # fragment_arr
      end # lambda block
      send_to_list.call(n_terms, @n_term_search_ion_types)
      send_to_list.call(c_terms, @c_term_search_ion_types)		
      @list = list
    end	
    def to_mgf(seq = nil)
      if seq.nil?
        seq = @sequence
        list = @list
      else
        list = generate_fragment_masses(seq)
      end
      intensity = 1000 # An arbitrary intensity value
      output_arr = []
      output_arr << %Q{COM=Project: In-silico Fragmenter\nBEGIN IONS\nPEPMASS=#{precursor_mass(seq, @max_charge)}\nCHARGE=#{@max_charge}+\nTITLE=Label: Sequence is #{seq}}
      list.sort_by{|a| a.mass}.each do |table_entry|
        #	TableEntry = Struct.new(:ion, :seq, :mass, :charge)
        output_arr << "#{"%.5f" % table_entry.mass }\t#{intensity}"
      end
      output_arr << "END IONS"
      File.open("#{seq}.mgf", 'w') {|o| o.print output_arr.join("\n") }
      output_arr.join("\n")
    end
    def graph(list = nil)
      list ? list : list = @list
      require 'rserve/simpler'
      robj = Rserve::Simpler.new
      hash = {}
      hash["mass"] = list.map(&:mass)
      hash["intensity"] = list.map{ 1000.0} # Hacky standard intensity value
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
  end
end
          
       
######### Testing stuff
if $0 == __FILE__
  require 'optparse'
  options = {charge_states: true, avg: false}
  ion_options = {}
  parser = OptionParser.new do |opts|
    opts.banner = "Usage: #{File.basename(__FILE__)} sequence [options]"
    opts.separator "Output: [Array] (containing fragment ion masses)"

    opts.on('--ion_type a,b,c,x,y,z', Array, "Select ion types (default is b,y)") do |t|
      arr = t.map{|a| a.downcase.to_sym}
      hash = {}
      arr.each {|a| hash[a] = true}
      ion_options[:ion_types] = hash
    end
    opts.on('--[no-]charge_states', "Turn on or off the charge state output") do |s|
      options[:charge_states] = s
    end
    opts.on('-a', '--avg', "Use average masses to calculate ions instead of monoisotopic masses") do |a|
      options[:avg] = a
    end
    if ARGV.size == 0
      puts opts
      exit
    end
    opts.on('--[no-]charge_states', "Turn off output of multiple charge states in list") do |s|
      options[:charge_states] = s
    end
    opts.on() do 

    end

    opts.on_tail('-h', '--help', "Show this message") do 
      puts opts
      exit
    end
  end.parse!  # OptionParser
  if ARGV.size >= 1
    f = Fragmenter.new(options, ion_options)
    f.fragment(ARGV.first)
    puts "I graphed these fragments and wrote them to #{f.graph} for you."
  end  
end # $0 == __FILE__


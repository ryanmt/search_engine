require 'mspire/mass'
class Molecule
  Neutron = 1.00866491600
  Cutoff = 0.001
  attr_accessor :charge_state, :name, :counts, :z, :monoisotopic
  alias :charge_state :z
  def initialize(hash, opts={})
    # This takes a hash of symbol letters keys representing the type of element and an integer value representing the count in that 
    @counts = hash
    @name = opts[:name] ? opts[:name] : nil
    @charge_state = opts[:charge_state] ? opts[:charge_state] : nil
    @z = opts[:z] ? opts[:z] : @charge_state
  end
# Formula here is like: C2BrH12O14
  def mass(formula=nil)
# Charge state changes would be reflected in giving additional MS::Mass::H_PLUS 
    Mspire::Mass.formula_to_exact_mass(formula) + @charge_state * Mspire::Mass::H_PLUS
  end
  def mass_from_counts
    strings = [0,4].map{|i| @counts.map{|k,v| "#{k.to_s.capitalize}#{v+i}" if v > 0}.join }
    strings.map {|s| Mspire::Mass.formula_to_exact_mass(s) }
  end
  def isotopic_distribution
    arrs = {}; freq_arrs = {}
    arr_size = mass_from_counts.last
    @counts.keys.each do |element|
      next if @counts[element] < 1
      arr = Array.new(arr_size) { 0 }
      arrs[element] = arr
    end
    require_relative 'isotopes'
    require_relative '../fft_R'
    transformer = FFT.new
    arrs.keys.each do |element|
      isotopes = ::Isotopes[element]
      isotopes.each do |arr|
        arrs[element][arr[2]] = arr[0]
      end
      freq_arrs[element] = transformer.fft( arrs[element]).map {|a| a**@counts[element]}
    end
    @monoisotopic = mass(counts_to_formula)
    output_size = @monoisotopic + @counts.values.inject(0){|s,v| v > s ? v : s} 
    output = Array.new(arr_size) {1}
    output = output.zip(*freq_arrs.values).map {|list| list.compact.reduce(:*)}
    filtered = transformer.ifft(output).each_with_index.map {|v,i| [v,i]}.select {|v| v.first > Cutoff }
    masses_and_raw_intensities = filtered.zip(filtered.size.times.map {|i| @monoisotopic + i*Neutron}).map do |arr|
      [arr.last,arr[0].first]
    end
  end
  def counts_to_formula
    outstring = ""
    @counts.each do |k,v|
      outstring << k.to_s.capitalize + v.to_s unless v < 1
    end
    outstring
  end
end

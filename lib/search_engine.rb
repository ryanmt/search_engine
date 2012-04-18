require 'data_structs'
require 'spectrum'
require 'search_engine/searcher'


def putsv(string)
  puts(string) if $VERBOSE
end

module MS
  module SearchEngine
    # Default Options Hash: Tolerances in ppm for now
    Options = {decoy: true, 
      precursor_mass_tolerance: 10, 
      ms2_mass_tolerance: 300, 
      write_pepdb: false,
      write_spectradb: false,
      xcorr_range: 75,
      missed_cleavages: 2 } 
  end
end

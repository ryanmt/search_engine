require 'mspire/bin'
require 'tools/ppm'
module MS
  class Bin < Mspire::Bin
    def self.create_bins_by_ppm(start_spot, end_spot, ppm)
      bins = [Bin.new(start_spot, PPM.shift_up(start_spot, ppm), true)]
      while bins.last.end <= end_spot
        bins << Bin.new(bins.last.end, PPM.shift_up(bins.last.end, ppm), true)
      end
      lbin = bins.last
# Right now, I just force the size of the last spot to the end_spot
      bins[-1] = Bin.new(lbin.begin, end_spot, false)
# Return the bins
      bins
    end #self.create_bins_by_ppm
  end # class
end #module

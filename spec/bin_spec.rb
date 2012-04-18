require 'spec_helper'

describe Bin do 
  it 'can handle creation of bins with variable widths (ppm widths)' do 
    bins = Bin.create_bins_by_ppm(200, 2000, 100)
    [bins.first, bins[-2]].map do |bin| 
      ppms = PPM.calculate(bin.end, bin.begin)
      ppms.should be_within(1e-2).of(100)
    end
  end
end 

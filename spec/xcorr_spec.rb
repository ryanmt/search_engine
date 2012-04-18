require 'spec_helper'

describe XCorr do 
  it 'shifts the sums' do 
    a = [0,1,2,3,4,5]
    putsv "Shift_summer 1"
    XCorr.bin_shifter(a, 1).should == [1,2,3,4,5,0]
    putsv "Shift_summer 2"
    XCorr.bin_shifter(a, 2).should == [2,3,4,5,0,1]
  end
  it 'shifts to the negative' do 
    a = [0,1,2,3,4,5]
    putsv "Shift_summer -1"
    XCorr.bin_shifter(a, -1).should == [5,0,1,2,3,4]
    putsv "Shift_summer -3"
    XCorr.bin_shifter(a, -3).should == [3,4,5,0,1,2]
  end
  it 'handles a full range shift' do 
    range = (200..2000).to_a
    XCorr.bin_shifter(range, 75)[0,50].should == (275...325).to_a
    XCorr.bin_shifter(range, -75)[0,50].should == (1926...1976).to_a
  end
  it 'generates an xcorr value' do 
    scans = (1..5).map do |i| 
      1500.times.map { rand > 0.8 ? rand : rand(100000) }
    end
    resp = XCorr.compare_to(scans, scans.first)
    p resp.map {|l| l/resp.first}
    resp.should == [5,0,0,0,0]
  end
end

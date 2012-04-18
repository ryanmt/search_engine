describe PPM do 
  it 'calculates' do 
    ans = PPM.calculate 445.120025, 445.11
    ans.should  be_within(1e-3).of(22.522)
  end
  it 'tests PPM.within?' do 
    ans = PPM.within? 445.120025, 445.11, 25
    ans.should == true
  end
  it 'gives a Range window' do 
    ans = PPM.mass_to_window(445.120025, 25)
    ans.is_a?(Range).should == true
    ans.should == Range.new(445.108896999375, 445.131153000625)
  end
  it 'gives the shift up' do 
    ans = PPM.shift_up(445.120025, 20)
    ans.should be_within(1e-4).of(445.1289)
  end
  it 'gives the shift down' do 
    ans = PPM.shift_down(445.120025, 20)
    ans.should be_within(1e-4).of(445.1111)
  end
end

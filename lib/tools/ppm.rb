class PPM
  def self.mass_to_window(mass, ppm)
    tolerance = ppm/1e6*mass.to_f
    Range.new(mass-tolerance, mass+tolerance)
  end

  def self.within?(mass1, mass2, ppm)
    mass_to_window(mass1, ppm).cover?(mass2)
  end

  def self.calculate(mass1, mass2)
    mass1, mass2 = [mass1, mass2].map(&:to_f)
    ppm = (mass1-mass2)/mass1*1e6
  end

  def self.shift_up(mass, ppm)
    mass_to_window(mass,ppm).end
  end

  def self.shift_down(mass, ppm)
    mass_to_window(mass,ppm).begin
  end
end

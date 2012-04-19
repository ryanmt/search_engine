module MS
  Proton = 1.00782503207
	def precursor_mass(seq, charge)
		mass = NTerm + CTerm
		seq.chars.map(&:to_sym).each do |residue|
			mass += MS::MonoResidueMasses[residue]
		end
		charge_state(mass, charge)
	end
		
	def charge_state(mass, charge)
		if charge > 0
			(mass + charge) / charge.to_f
		else
			(mass - charge) / charge.to_f
		end
	end

	IonTypeMassDelta = {
		a:  (- 29.00273),
		a_star: (-(29.00273+17.02654)),
		a_not:	(-(17.02654 + 29.00273+18.01056)),
		b: (-1.00782),
		b_star: ( - 1.00782 - 17.02654),
		b_not: (-17.02654 - 1.00782 - 18.01056),
		c: (16.01872),
		x: (27.99491 - 1.00782),
		y: (1.00782),
		y_star: (1.00782  - 17.02654),
		y_not: (1.00782 - 18.01056),
		z: (- 16.01872)
	}
	
	NTerm =	1.00782

	CTerm =	27.99491 - 10.9742

	MonoResidueMasses = {
		:A => 71.037114, 
		:R => 156.101111,
		:N => 114.042927,
		:D => 115.026943,
		:C => 103.009185,
		:E => 129.042593,
		:Q => 128.058578, 
		:G => 57.021464,
		:H => 137.058912,
		:I => 113.084064,
		:L => 113.084064,
		:K => 128.094963,
		:M => 131.040485,
		:F => 147.068414,
		:P => 97.052764,
		:S => 87.032028,
		:T => 101.047679,
		:U => 150.95363,
		:W => 186.079313,
		:Y => 163.06332,
		:V => 99.068414,
    :* => 118.805716,
    :B => 172.048405,
    :X => 118.805716,
    :Z => 128.550585
	}
  AvgResidueMasses = {
    :* => 118.88603, 
    :A => 71.0779, 
    :B => 172.1405, 
    :C => 103.1429, 
    :D => 115.0874, 
    :E => 129.11398, 
    :F => 147.17386, 
    :G => 57.05132, 
    :H => 137.13928, 
    :I => 113.15764, 
    :K => 128.17228, 
    :L => 113.15764, 
    :M => 131.19606, 
    :N => 114.10264, 
    :O => 211.28076, 
    :P => 97.11518, 
    :Q => 128.12922, 
    :R => 156.18568, 
    :S => 87.0773, 
    :T => 101.10388, 
    :U => 150.0379, 
    :V => 99.13106, 
    :W => 186.2099, 
    :X => 118.88603, 
    :Y => 163.17326, 
    :Z => 128.6231
  }
end
################
=begin
Formula: H3N1

Monoisotopic mass :     17.02654

Formula: C1H1O1

Monoisotopic mass :     29.00273

Formula: H2O1

Monoisotopic mass :     18.01056

Formula: H1

Monoisotopic mass :      1.00782

Formula: H2N1

Monoisotopic mass :     16.01872

Formula: C1O1

Monoisotopic mass :     27.99491
=end

module Ms
	def precursor_mass(seq, charge)
		mass = NTerm + CTerm
		seq.chars.map(&:to_sym).each do |residue|
			mass += Ms::MonoResidueMasses[residue]
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
################
# This is straight from my pI calculator, and adds the fxn of calculating a maximum charge state for the total peptide, given the sequence.
#
#	Usage:  charge_at_pH(identify_potential_charges(peptide_sequence), pH_desired)
	PepCharges = Struct.new(:seq, :n_term, :c_term, :y_num, :c_num, :k_num, :h_num, :r_num, :d_num, :e_num, :pi)
	def identify_potential_charges(str)
		string = str.upcase
		first = string[0]; last = string[-1]
		puts string if first.nil? or last.nil?
		begin
			out = PepCharges.new(string, PkTable[first.to_sym][0], PkTable[last.to_sym][1], 0, 0, 0 ,0 ,0 ,0, 0)
		rescue NoMethodError
			abort string
		end
		string.chars.each do |letter|
			case letter
				when "Y" 
					out.y_num += 1 
				when "C"
					out.c_num += 1
				when "K"
					out.k_num += 1 
				when "H"
					out.h_num += 1
				when "R"
					out.r_num += 1
				when "D"
					out.d_num += 1
				when "E"
					out.e_num += 1
			end
		end
		out
	end # Returns the PepCharges structure

	PkTable = {
		:K => [2.18,8.95,10.53], 
		:E => [2.19,9.67,4.25], 
		:D => [1.88,9.60,3.65], 
		:H => [1.82,9.17,6.00],
		:R => [2.17,9.04,12.48],
		:Q => [2.17,9.13,nil],
		:N => [2.02,8.80,nil],
		:C => [1.96,10.28,8.18],
		:T => [2.11,9.62,nil],
		:S => [2.21,9.15,nil],
		:W => [2.38,9.39,nil],
		:Y => [2.20,9.11,10.07],
		:F => [1.83,9.13,nil],
		:M => [2.28,9.21,nil],
		:I => [2.36,9.68,nil],
		:L => [2.36,9.60,nil],
		:V => [2.32,9.62,nil],
		:P => [1.99,10.96,nil],
		:A => [2.34,9.69,nil],
		:G => [2.34,9.60,nil],
# These are the fringe cases... B and Z... Jerks, these are harder to calculate pIs
		:B => [1.95,9.20,3.65],
		:Z => [2.18,9.40,4.25],
		:X => [2.20,9.40,nil],
		:U => [1.96,10.28,5.20] # Unfortunately, I've only found the pKr for this... so I've used Cysteine's values.
	}

	def charge_at_pH(pep_charges, pH)
		charge = 0
		charge += -1/(1+10**(pep_charges.c_term-pH))
		charge += -pep_charges.d_num/(1+10**(PkTable[:D][2]-pH))
		charge += -pep_charges.e_num/(1+10**(PkTable[:E][2]-pH))
		charge += -pep_charges.c_num/(1+10**(PkTable[:C][2]-pH))
		charge += -pep_charges.y_num/(1+10**(PkTable[:Y][2]-pH))
		charge += 1/(1+10**(pH - pep_charges.n_term))
		charge += pep_charges.h_num/(1+10**(pH-PkTable[:H][2]))
		charge += pep_charges.k_num/(1+10**(pH-PkTable[:K][2]))
		charge += pep_charges.r_num/(1+10**(pH-PkTable[:R][2]))
		charge
	end
############################################################
end
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

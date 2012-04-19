module MS
  class ChargeCalculator
# This is straight from my pI calculator, and adds the fxn of calculating a maximum charge state for the total peptide, given the sequence.
#
#	Usage:  charge_at_pH(identify_potential_charges(peptide_sequence), pH_desired)
    PepCharges = Struct.new(:seq, :n_term, :c_term, :y_num, :c_num, :k_num, :h_num, :r_num, :d_num, :e_num, :pi)
    def self.identify_potential_charges(str)
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

    def self.charge_at_pH(pep_charges, pH)
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
    end #charge_at_pH
  end # class ChargeCalculator
end


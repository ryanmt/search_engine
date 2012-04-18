module MS
  module DataStructs
    class Peptide < String
      attr_accessor :proteins, :spectrum
    end

    class Match
      attr_accessor :xcorr, :ppm, :experimental_spectrum, :theoretical_spectrum
    end

    class SpectrumObjects
      attr_accessor :mzs, :intensities, :precursor_z, :precursor_mass
      def precursor_neutral_mass
        # neutralize the charge to calculate the base mass for matching
      end

      def mzs_and_intensities
        @mzs.zip(@intensities)
      end

      def normalize
        # Necessary?
      end
    end
    class TheoreticalSpectrum < SpectrumObjects
      attr_accessor :peptide
      def initialize(mzs, intensities = Array.new(mzs.size, 100.0), peptide)
        @mzs, @intensities, @peptide = mzs, intensities, peptide
      end
    end
  end
end

module MS
require 'tools/fragmenter/masses'
  module DataStructs
    class Peptide < String
      attr_accessor :proteins, :spectrum
      def precursor_neutral_mass
        MS::Isotoper.precursor_mass(self.to_s)
      end
    end

    class Match
      attr_accessor :xcorr, :ppm, :experimental_spectrum, :theoretical_spectrum, :qvalue
      def initialize(xcorr, spectrum, ppm)
        @xcorr, @theoretical_spectrum, @ppm = xcorr, spectrum, ppm
      end
    end

    class SpectrumObjects
      attr_accessor :mzs, :intensities, :precursor_z, :precursor_mass
      def precursor_neutral_mass 
        (@precursor_mass - MS::Proton)* @precursor_z
      end

      def mzs_and_intensities
        @mzs.zip(@intensities)
      end

      def normalize
        # Necessary?
      end
    end
    class TheoreticalSpectrum < SpectrumObjects
      attr_accessor :peptide, :bins
      def initialize(mzs, intensities = Array.new(mzs.size, 100.0), peptide)
        @mzs, @intensities, @peptide = mzs, intensities, peptide
        @precursor_z = 0
        @precursor_mass = @peptide.precursor_neutral_mass
        @bins = []
      end
    end #TheoreticalSpectrum
    class ExperimentalSpectrum < SpectrumObjects
      attr_accessor :bins, :match
      def initialize(mzs, intensities, precursor_mass, precursor_charge)
        @mzs, @intensities, @precursor_mass, @precursor_z = mzs, intensities, precursor_mass, precursor_charge
        @bins = []
      end
    end #ExperimentalSpectrum
  end
end

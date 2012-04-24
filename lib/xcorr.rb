
class XCorr
  def self.compare_to(theoretical_scans, experimental_scan) 
    num = MS::SearchEngine::Options[:xcorr_range]
    range = ((num * -1)..(-1)).to_a + (1..num).to_a
    correction_arrays = range.map do |i|
      bin_shifter(experimental_scan, i)
    end #.transpose.inject(:+).map {|a| a/150}
    y_prime = experimental_scan.zip(correction_arrays.transpose.inject(:+) ).map {|a| a.first - a.last/150 }
    theoretical_scans.map do |scan| 
      dot_product(scan, y_prime)
    end
  end

  # Comparative benchmarking @ http://blade.nagaokaut.ac.jp/cgi-bin/scat.rb/ruby/ruby-talk/277157
  # Here then is the winner (@ 1e6 sized arrays of floats less than 1): 
  def self.dot_product(arr1, arr2)
    sum = 0
    for i in 0...arr1.size
      sum += arr1[i] * arr2[i]
    end
    sum
  end
  def self.bin_shifter(arr, bin_shift)
    rollover = arr.size
    arr.each_index.map do |i| 
#      #DEBUGGING
#        putsv "i+bin_shift >= rollover: #{i+bin_shift} >= #{rollover}"
#        putsv "rollover-(i+bin_shift): #{rollover - (i + bin_shift)}"
#        putsv "(i+bin_shift): #{(i + bin_shift)}"
#      #Output
      i+bin_shift >= rollover ? arr[(rollover-(i+bin_shift)).abs] : arr[i+bin_shift]
    end
  end
end #class XCorr


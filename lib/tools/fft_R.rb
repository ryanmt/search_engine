require 'rserve/simpler'

def arrays_to_complex_array(real, img)
  raise StandardError unless real.size == img.size
  real.each_with_index.map {|n, i| Complex(n, img[i]) }
end
def complex_array_to_arrays(complex_arr)
  img = []
  real = complex_arr.map {|n| img.push n.imag; n.real}
  [real, img]
end
class FFT
  def initialize(connection=nil)
    @r = connection || Rserve::Simpler.new 
  end

  def fft(arr, options = {})
    defaults = {normalize: true}
    opts = defaults.merge(options)
    @r.command(a: arr) {"x <- fft(a)"}
    real = @r.converse "Re(x)"
    img = @r.converse "Im(x)"
    arrays_to_complex_array(real, img)
  end

  def ifft(complex_arr, options = {} )
    defaults = {normalize: true}
    opts = defaults.merge(options)
    norm = opts[:normalize] ? complex_arr.size : 1
    arr = complex_array_to_arrays(complex_arr)
# Join the two arrays in R
    @r.command(a: arr) do 
      %Q{ input_arr = NULL  
        for(i in seq(along=a[[1]])) {
          input_arr = c(input_arr, complex(real=a[[1]][i], imaginary=a[[2]][i]))
        }
      } 
    end
    @r.command {"x <- fft(input_arr, inverse=TRUE)/#{norm}"}
    real = @r.converse "Re(x)"
    #imag = @r.converse "Im(x)"
    real # [real, imag]
  end

end # FFT

if __FILE__ == $0
  f = FFT.new
  freq = f.fft([1,2,3,4])
  p freq
  norm = f.ifft(freq)
  p norm
end

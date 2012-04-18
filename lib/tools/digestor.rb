#!/usr/bin/env ruby
class Digestor
  DEFAULTS = {:missed_cleavages => 0, :minimum_length => 1}
  def initialize(enzyme)
    @regexp = Enzymes[enzyme.downcase]
  end
  def digest(string, opts = {}) # Returns an array of chomped sequences
    opts = DEFAULTS.merge(opts)
    holder = string.upcase.scan(@regexp).flatten.compact
    outs = holder
    if opts[:missed_cleavages] != 0
      outs << (1..(opts[:missed_cleavages]+1)).map do |cons|
        holder.each_cons(cons).map(&:join)
      end
    end
    outs = outs.flatten.map{|a| a if a.length >= opts[:minimum_length]}.compact
  end

  Enzymes = {'trypsin' => Regexp.new('(\D*?[KR])(?=[^P])|([^P].*?\b)'), 
    "trypsin/p" => Regexp.new('(\D*?[KR])|(.*?\b)'), 
    "arg-c" => Regexp.new('(\w+?[R])(?=[^P])|([^P].*?\b)'), 
    "asp-n" => Regexp.new('\w+?(?=[BD]|$)'),
    "chymotrypsin" => Regexp.new('(\D*?[FYWL])(?=[^P])|([^P].*?\b)'), 
    "cnbr+trypsin" => Regexp.new('(\D*?(?:(?:[KR](?=[^P]))|[M]))|(.*?\b)'),
    "cnbr" => Regexp.new('(\D*?M)|(.*?\b)')
  }

end

if __FILE__ == $0
  if ARGV.size < 3
    puts "Usage: #{File.basename(__FILE__)} <missed_cleavages> <enzyme> <seq> ... <seq(n)>"
    puts "Outputs: seq joined by comma, one per line"
    exit
  end
  cleavages = ARGV.shift
  enzyme = ARGV.shift
  d = Digestor.new(enzyme)
  ARGV.each do |seq|
    puts d.digest(seq, missed_cleavages: cleavages).join(", ")
  end
end

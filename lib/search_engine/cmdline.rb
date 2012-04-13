require 'optparse'
module MS
  module SearchEngine
    attr_accessor :options
    @options = {}
    class Cmdline
      def self.run(args)
        options = MS::SearchEngine.options
        parser = OptionParser.new do |opts|
          opts.on_tail('-h', '--help') do 
            puts opts
            return
          end
          opts.banner = "Usage: ruby search-spectra.rb input.fasta input.mzML"
        end
        parser.parse!(args)
        if args.size != 2
          puts parser.banner
          puts parser.summarize
          return
        end
      end #self.run
    end
  end
end

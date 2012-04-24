require 'optparse'
require_relative '../search_engine'
module MS
  module SearchEngine
    attr_accessor :options
    @options = {}
    class Cmdline
      def self.run(args)
        options = MS::SearchEngine::Options
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
        searcher = MS::SearchEngine::Searcher.new(args[0])
        searcher.search(args[1])
      end #self.run
    end
  end
end


if $0 == __FILE__
  $VERBOSE = 3
  p ARGV
  if ARGV.first == 'skip'
    #Loading from yaml takes 5 minutes, same as writing it out to YAML
    args = ["/home/ryanmt/Dropbox/HW/proteomics/bootcamp/search_engine/mystery_fasta.spectra.yml", "/home/ryanmt/Dropbox/HW/proteomics/bootcamp/search_engine/mystery_run.mzML"]
    putsv args.first
  else
    args = ["/home/ryanmt/Dropbox/HW/proteomics/bootcamp/search_engine/mystery_fasta.fasta", "/home/ryanmt/Dropbox/HW/proteomics/bootcamp/search_engine/mystery_run.mzML"]
    putsv args.first
  end
  MS::SearchEngine::Cmdline.run(args)
end

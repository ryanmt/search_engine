require 'fileutils'
require 'mspire/fasta'
module MS
  class Decoy
    Defaults = {concatenate: false, type: :reverse, :prefix => 'DECOY_'}
    def initialize(opts = {} )
      @opts = Defaults.merge(opts)
      self
    end
    def generate(input_file, output_file = nil)
      type_str = "_#{@opts[:type].to_s == 'randomize' ? 'shuffled' : 'reversed'}"
      outfile_addin = @opts[:concatenate] ? "_concatenated" + type_str : type_str
      output_file ||= (File.basename(input_file, '.fasta') + '_'  + outfile_addin + ".decoy.fasta").gsub(/_{2,}/, '_')
      entries = {}
      Mspire::Fasta.open(input_file) do |fasta| 
        fasta.each do |entry| 
          fix = entry.accession ? entry.accession : entry.entry_id
          entries[entry.header.sub(fix, @opts[:prefix] + fix)] = entry.sequence.send(@opts[:type])
        end
      end
      if @opts[:concatenate]
        FileUtils.cp(input_file, output_file)
      end  
      File.open(output_file, 'w') {|out| out.puts entries.map {|k,v| ">#{k}\n#{v}\n"} }
      output_file
    end
    class Cmdline
      class << self
        def run(argv)
          require 'optparse'
          options = {type: :reverse, concatenate: false}
          parser = OptionParser.new do |opts|
            opts.banner = "Usage: ruby dbdecoy.rb [options]"
            opts.on('--shuffle', "Shuffle the protein sequence instead of reversing it") do 
              options[:type] = :randomize
            end
            opts.on('--prefix STRING',String, "Change the accession modifying prefix ('DECOY_') to given string") do |str|
              options[:prefix] = str
            end
            opts.on('-c', '--concatenate', "Generates a concatenated decoy+real database instead of a separate, decoy only database.") do |c|
              options[:concatenate] = c
            end
            opts.on_tail('-h', '--help', 'Show this message') do 
              puts opts
              return
            end
          end
          parser.parse!(argv)
          if argv.size == 0
            puts parser.banner
            puts parser.summarize
            return
          end

          outf = Ryan::MS::Search::Decoy.new(options).generate(argv.first)
          [options, outf]
        end
      end
    end # Cmdline
  end # Decoy
end # MS

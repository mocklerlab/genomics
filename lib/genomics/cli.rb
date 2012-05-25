require 'thor'
require 'genomics'

module Genomics
  class CLI < Thor
    desc "identify", "Identifies the reciprocally best alignments between the supplied files."
    method_option :files, type: :array, required: true, aliases: '-f'
    method_option :output, type: :string, aliases: '-o'
    def identify
      command_options = {}
      command_options[:output_file] = options[:output] if options[:output]
      puts "#{Genomics::Alignment::RBB.identify(options[:files], command_options)} Reciprocal Best Alignments Identified"
    end
    
    desc "blastx", "Takes the supplied alignment file and generates a BLASTX GFF3 file from the results."
    method_option :input, type: :string, required: true, aliases: '-i'
    def blastx
      puts "BLASTX GFF3 successfully created." if Genomics::Alignment::BLASTX.transform(options[:input])
    end
    
    desc "est", "Takes the supplied alignment file and generates an EST GFF3 file from the results."
    method_option :input, type: :string, required: true, aliases: '-i'
    def est
      puts "EST GFF3 successfully created." if Genomics::Alignment::EST.transform(options[:input])
    end
  end
end
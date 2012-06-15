require 'thor'
require 'genomics'

module Genomics
  class Ests < Thor
    desc "test", "Identifies the reciprocally best alignments between the supplied files."
    method_option :files, type: :array, required: true, aliases: '-f'
    method_option :output, type: :string, aliases: '-o'
    def test
      command_options = {}
      command_options[:output_file] = options[:output] if options[:output]
      puts "#{Genomics::Operation::RBB.identify_orthologs(options[:files], command_options)} Reciprocal Best Alignments Identified"
    end
  end
end
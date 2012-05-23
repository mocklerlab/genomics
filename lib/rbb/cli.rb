require 'thor'
require 'rbb'

module Genomic
  module RBB
    class CLI < Thor
      desc "identify", "Identifies the reciprocally best alignments between the supplied files."
      method_option :files, type: :array, required: true, aliases: '-f'
      def identify
        puts "#{Genomics::RBB::Tasker.identify(*options[:files])} Reciprocal Best Alignments Identified"
      end
    end
  end
end
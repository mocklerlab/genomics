require 'thor'
require 'genomics'

require 'genomics/generators/test'

module Genomics
  class CLI < Thor
    desc "identify", "Identifies the reciprocally best alignments between the supplied files."
    method_option :files, type: :array, required: true, aliases: '-f'
    method_option :output, type: :string, aliases: '-o'
    def identify
      command_options = {}
      command_options[:output_file] = options[:output] if options[:output]
      puts "#{Genomics::Operation::RBB.identify_orthologs(options[:files], command_options)} Reciprocal Best Alignments Identified"
    end
  
    desc "rbb", "Identifies putative orthologs between proteomes, which are reciprocal best blasts of each other."
    method_option :protein_files, type: :array, required: true, aliases: '-p'
    method_option :database_files, type: :array, required: true, aliases: '-d'
    method_option :alignment_file_dir, type: :string, aliases: '-a', desc: 'The path to the directory where intermediate alignment files should be written.'
    method_option :output, type: :string, aliases: '-o'
    method_option :e_value, type: :numeric, aliases: '-e', desc: 'The evalue used as the cutoff for the alignment.'
    method_option :threads, type: :numeric, aliases: '-t', desc: 'The number of threads to use in each BLAST alignment.'
    def rbb
      proteomes = []
      options[:protein_files].each_with_index do |filepath, index|
        proteomes << { file: filepath, database: options[:database_files][index] }
      end
      command_options = { blast_options: {} }
      command_options[:alignment_file_dir] = options[:alignment_file_dir] if options[:alignment_file_dir]
      command_options[:blast_options][:threads] = options[:threads] if options[:threads]
      command_options[:blast_options][:e_value] = options[:e_value] if options[:e_value]
      
      puts "#{Genomics::Operation::RBB.perform(proteomes, command_options)} Reciprocal Best Alignments Identified"
    end
    
    desc "blastx", "Takes the supplied alignment file and generates a BLASTX GFF3 file from the results."
    method_option :input, type: :string, required: true, aliases: '-i'
    def blastx
      puts "BLASTX GFF3 successfully created." if Genomics::Alignment::BLASTX.transform(options[:input])
    end
    
    desc "est", "Takes the supplied EST FASTA file and the genome FASTA file and generates a GFF3 file from the alignment."
    method_option :est, type: :string, required: true, aliases: '-e'
    method_option :genome, type: :string, required: true, aliases: '-g'
    method_option :alignment_file_dir, type: :string, aliases: '-a', desc: 'The path to the directory where intermediate alignment files should be written.'
    # method_option :threads, type: :numeric, aliases: '-t'
    def est
      command_options = { alignment: {} }
      command_options[:alignment][:alignment_file_dir] = options[:alignment_file_dir] if options[:alignment_file_dir]
      puts "EST GFF3 successfully created." if Genomics::Operation::AlignESTs.perform(options[:est], options[:genome], command_options)
    end
  end
end
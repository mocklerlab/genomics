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
    
    ### Peptide BLASTX Commands ###
    
    desc "blastx", "Takes the query genome and the target proteome generating a GFF3 file describing the BLASTX alignment between the two."
    method_option :protome, type: :string, required: true, aliases: '-p'
    method_option :genome, type: :string, required: true, aliases: '-g'
    method_option :alignment, type: :string, aliases: '-a'
    method_option :output, type: :string, aliases: '-o'
    method_option :threads, type: :numeric, aliases: '-t', default: 1
    method_options id_prefix: :string, e_value: :numeric
    def blastx
      arguments = Hash[options.select{ |(key, value)| value }]
      arguments[:proteome_file] = arguments.delete(:proteome)
      arguments[:genome_file] = arguments.delete(:genome)
      arguments[:alignment_file] = arguments[:alignment] if arguments[:alignment]
      arguments[:output_file] = arguments[:output] if arguments[:output]
      
      # Execute the operation
      blastx = Genomics::Operation::BLASTX.new(arguments)
      puts "BLASTX GFF3 successfully created." if blastx.perform
    end
    
    desc "blastx", "Takes the query genome and the target proteome generating a BLASTX alignment between the two."
    method_option :protome, type: :string, required: true, aliases: '-p'
    method_option :genome, type: :string, required: true, aliases: '-g'
    method_option :output, type: :string, aliases: '-o'
    method_option :threads, type: :numeric, aliases: '-t', default: 1
    method_options e_value: :numeric
    def blastx_align
      arguments = Hash[options.select{ |(key, value)| value }]
      arguments[:task] = :align
      arguments[:proteome_file] = arguments.delete(:proteome)
      arguments[:genome_file] = arguments.delete(:genome)
      arguments[:alignment_file] = arguments[:output] if arguments[:output]
      
      # Execute the operation
      blastx = Genomics::Operation::BLASTX.new(arguments)
      puts "BLASTX alignment successfully created." if blastx.perform
    end
    
    desc "blastx_cluster", "Takes an alignment between a query genome and target proteome generating a GFF3 file."
    method_option :alignment, type: :string, required: true, aliases: '-a'
    method_option :output, type: :string, aliases: '-o'
    method_option :threads, type: :numeric, aliases: '-t'
    method_options id_prefix: :string
    def transcripts_cluster
      arguments = Hash[options.select{ |(key, value)| value }]
      arguments[:task] = :cluster
      arguments[:alignment_file] = arguments.delete(:alignment)
      arguments[:output_file] = arguments[:output] if arguments[:output]
      
      # Execute the operation
      blastx = Genomics::Operation::BLASTX.new(arguments)
      puts "BLASTX GFF3 successfully created." if blastx.perform
    end
    
    ### Transcript Alignment Commands ###
    
    desc "transcripts", "Takes the supplied transcript FASTA file and the genome FASTA file and generates a GFF3 file from the alignment."
    method_option :transcripts, type: :string, required: true, aliases: '-t'
    method_option :genome, type: :string, required: true, aliases: '-g'
    method_option :alignment, type: :string, aliases: '-a'
    method_option :output, type: :string, aliases: '-o'
    method_option :threads, type: :numeric, aliases: '-t', default: 1
    method_options id_prefix: :string, e_value: :numeric
    def transcripts
      arguments = Hash[options.select{ |(key, value)| value }]
      arguments[:transcripts_file] = arguments.delete(:transcripts)
      arguments[:genome_file] = arguments.delete(:genome)
      arguments[:alignment_file] = arguments[:alignment] if arguments[:alignment]
      arguments[:output_file] = arguments[:output] if arguments[:output]
      
      # Execute the operation
      aligner = Genomics::Operation::TranscriptAligner.new(arguments)
      puts "Transcript GFF3 successfully created." if aligner.perform
    end

    desc 'transcripts_align', 'Takes the supplied transcript FASTA file and the genome FASTA file and aligns them via BLAT.'
    method_option :transcripts, type: :string, required: true, aliases: '-t'
    method_option :genome, type: :string, required: true, aliases: '-g'
    method_option :output, type: :string, aliases: '-o'
    method_option :threads, type: :numeric, aliases: '-t'
    def transcripts_align
      arguments = Hash[options.select{ |(key, value)| value }]
      arguments[:task] = :align
      arguments[:transcripts_file] = arguments.delete(:transcripts)
      arguments[:genome_file] = arguments.delete(:genome)
      arguments[:alignment_file] = arguments[:output] if arguments[:output]
      
      aligner = Genomics::Operation::TranscriptAligner.new(arguments)
      puts "Alignment completed." if aligner.perform
    end
    
    desc 'transcripts_cluster', 'Takes the results of an alignment of transcripts against a genome and aggregates them into clusters ideally indicating intron/exon structure.'
    method_option :alignment, type: :string, required: true, aliases: '-a'
    method_option :output, type: :string, aliases: '-o'
    method_option :threads, type: :numeric, aliases: '-t'
    method_options id_prefix: :string, e_value: :numeric
    def transcripts_cluster
      arguments = Hash[options.select{ |(key, value)| value }]
      arguments[:task] = :cluster
      arguments[:alignment_file] = arguments.delete(:alignment)
      arguments[:output_file] = arguments[:output] if arguments[:output]
      
      aligner = Genomics::Operation::TranscriptAligner.new(arguments)
      puts "Clustering completed." if aligner.perform
    end
  end
end
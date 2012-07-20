require 'thor'
require 'genomics'

module Genomics
  class CLI < Thor
    
    ### Reciprocal Best Blast Commands ###
    
    desc "rbb", "Identifies putative orthologs between proteomes, which are reciprocal best blasts of each other."
    method_option :protein_files, type: :array, required: true, aliases: '-p', desc: 'A list of protein fasta files.'
    method_option :database_files, type: :array, required: true, aliases: '-d', desc: 'A list of corresponding BLASTDB files.'
    method_option :alignment_file_dir,  type: :string, 
                                        aliases: '-a', 
                                        desc: 'The path to the directory where intermediate alignment files should be written.'
    method_option :output, type: :string, aliases: '-o'
    method_option :e_value, type: :numeric, aliases: '-e', desc: 'The evalue used as the cutoff for the alignment.'
    method_option :threads, type: :numeric, aliases: '-t', desc: 'The number of threads to use in each BLAST alignment.'
    method_option :reciprocal, type: :boolean, default: true, desc: 'If set, only reciprocal best results will be retained.'
    method_option :detailed, type: :boolean, default: false, desc: 'If set, additional details about the alignment will be retained.'
    def rbb
      # Generate the arguments
      arguments = prepare_options(options)
      protein_files = arguments.delete(:protein_files)
      database_files = arguments.delete(:database_files)
      
      # Construct the proteomes data structure
      proteomes = []
      protein_files.each_with_index do |filepath, index|
        proteomes << { file: filepath, database: database_files[index] }
      end
      arguments[:proteomes] = proteomes
      
      rbb = Genomics::Operation::RBB.new(arguments)
      if results_count = rbb.perform
        puts "RBB successfully performed. (#{results_count} results identified)" 
      end
    end
    
    desc "rbb_align", "Performs the BLAST alignment between the proteomes generating input for RBB identification."
    method_option :protein_files, type: :array, required: true, aliases: '-p', desc: 'A list of protein fasta files.'
    method_option :database_files, type: :array, required: true, aliases: '-d', desc: 'A list of corresponding BLASTDB files.'
    method_option :alignment_file_dir,  type: :string, 
                                        required: true,
                                        aliases: '-a', 
                                        desc: 'The path to the directory where intermediate alignment files should be written.'
    method_option :e_value, type: :numeric, aliases: '-e', desc: 'The evalue used as the cutoff for the alignment.'
    method_option :threads, type: :numeric, aliases: '-t', desc: 'The number of threads to use in each BLAST alignment.'
    def rbb_align
      # Generate the arguments
      arguments = prepare_options(options.merge(task: :align))
      protein_files = arguments.delete(:protein_files)
      database_files = arguments.delete(:database_files)
      
      # Construct the proteomes data structure
      proteomes = []
      protein_files.each_with_index do |filepath, index|
        proteomes << { file: filepath, database: database_files[index] }
      end
      arguments[:proteomes] = proteomes
      
      # Execute the operation
      rbb = Genomics::Operation::RBB.new(arguments)
      puts "RBB alignment successfully performed." if rbb.perform
    end
    
    desc "rbb_ortholog", "Identifies putative orthologs from pairs of BLAST alignment files."
    method_option :alignment_files, type: :array, required: true, aliases: '-f', desc: 'A list of alignment files.'
    method_option :output, type: :string, aliases: '-o'
    method_option :reciprocal, type: :boolean, default: true, desc: 'If set, only reciprocal best results will be retained.'
    method_option :detailed, type: :boolean, default: false, desc: 'If set, additional details about the alignment will be retained.'
    def rbb_ortholog
      # Generate the arguments
      arguments = prepare_options(options)
      
      rbb = Genomics::Operation::RBB.new(arguments)
      if results_count = rbb.perform
        puts "RBB successfully performed. (#{results_count} results identified)" 
      end
    end
    
    ### Peptide BLASTX Commands ###
    
    desc "blastx", "Takes the query genome and the target proteome generating a GFF3 file describing the BLASTX alignment between the two."
    method_option :proteome, type: :string, required: true, aliases: '-p'
    method_option :genome, type: :string, required: true, aliases: '-g'
    method_option :alignment, type: :string, aliases: '-a'
    method_option :output, type: :string, aliases: '-o'
    method_option :threads, type: :numeric, aliases: '-t', default: 1
    method_options id_prefix: :string, e_value: :numeric
    def blastx
      arguments = prepare_options(options, { proteome: :proteome_file, 
                                             genome: :genome_file, 
                                             alignment: :alignment_file, 
                                             output_file: :output })
      
      # Execute the operation
      blastx = Genomics::Operation::BLASTX.new(arguments)
      puts "BLASTX GFF3 successfully created." if blastx.perform
    end
    
    desc "blastx_align", "Takes the query genome and the target proteome generating a BLASTX alignment between the two."
    method_option :protome, type: :string, required: true, aliases: '-p'
    method_option :genome, type: :string, required: true, aliases: '-g'
    method_option :output, type: :string, aliases: '-o'
    method_option :threads, type: :numeric, aliases: '-t', default: 1
    method_options e_value: :numeric
    def blastx_align
      arguments = prepare_options(options.merge(task: :align), { proteome: :proteome_file, 
                                                                 genome: :genome_file, 
                                                                 alignment: :alignment_file })
      # Execute the operation
      blastx = Genomics::Operation::BLASTX.new(arguments)
      puts "BLASTX alignment successfully created." if blastx.perform
    end
    
    desc "blastx_cluster", "Takes an alignment between a query genome and target proteome generating a GFF3 file."
    method_option :alignment, type: :string, required: true, aliases: '-a'
    method_option :output, type: :string, aliases: '-o'
    method_option :threads, type: :numeric, aliases: '-t'
    method_options id_prefix: :string
    def blastx_cluster
      arguments = prepare_options(options.merge(task: :cluster), alignment: :alignment_file, output: :output_file )
      
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
      arguments = Hash[options.to_a.select{ |(key, value)| !value.nil? }]
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
      arguments = Hash[options.to_a.select{ |(key, value)| !value.nil? }]
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
      arguments = Hash[options.to_a.select{ |(key, value)| !value.nil? }]
      arguments[:task] = :cluster
      arguments[:alignment_file] = arguments.delete(:alignment)
      arguments[:output_file] = arguments[:output] if arguments[:output]
      
      aligner = Genomics::Operation::TranscriptAligner.new(arguments)
      puts "Clustering completed." if aligner.perform
    end
    
    private
    
    # Takes a hash of options, removing keys with nil values and updates keys according to the map.
    #
    # * *Args*    :
    #   - +options+ -> The Hash of options
    #   - +name_map+ -> An optional Hash the specifies criteria for changing keys.
    # * *Returns* :
    #   - A Hash
    #
    def prepare_options(options, name_map = {})
      options_array = options.to_a.select { |(key, value)| !value.nil? }
      options_array.map! do |(key, value)|
        name_map[key.to_sym] ? [name_map[key.to_sym], value] : [key.to_sym, value]
      end
      
      Hash[options_array]
    end
  end
end
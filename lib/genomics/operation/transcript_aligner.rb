require 'pathname'

module Genomics
  module Operation
    # This class handles the operation of taking a set of transcripts e.g. (ESTs) and aligning them against a genomic scale sequence.
    #
    # * *Attributes* :
    #   - +task+ -> A Symbol specifying the particular step to be conducted.  Possible values include :align, :cluster, 
    #   :all (Default: :all)
    #   - +alignment_file+ -> A String specifying the path the alignment file.
    #   - +verbose+ -> A boolean that sets whether the current operation is printed to STDOUT (Default: true). 
    #   - +threads+ -> An integer that sets the maximum number of threads that can be used (Default: 1).
    #
    # * *Alignment Attributes* :
    #   - +transcripts_file+ -> A String specifying the path to the FASTA file with the transcript sequences.
    #   - +genome_file+ -> A String specifying the path to the FASTA file with the genomic sequence. 
    #   - +blat_path+ -> A String explicitly specifying the path to BLAT.
    #   - +blat_options+ -> A Hash of options to be passed directly to the Alignment::BLAT object used to conduct the alignment.
    #
    # * *Cluster Attributes* :
    #   - +output_file+ -> A String specifying the path of the generated GFF3 file.
    #   - +id_prefix+ -> A String that will be used as a prefix for the ID field in the resultant GFF3 file (Default: 'transcript_'). 
    #   - +e_value+ -> A float specifying the minimum EValue that a hit must have in order to be clustered (Default: 1e^-5).
    class TranscriptAligner
      
      attr_accessor :task, :verbose, :threads, :alignment_file,
                    :transcripts_file, :genome_file, :blat_path, :blat_options,
                    :output_file, :id_prefix, :e_value 
      
      def initialize(options = {})
        options = { task: :all, verbose: true, threads: 1, 
                    blat_options: {},
                    id_prefix: 'transcript_', e_value: 1e-5 }.merge(options)
        
        options.each do |name, value|
          send("#{name}=", value)
        end
      end
      
      # Carries out the process of generating a GFF3 file from alignments of the transcripts against the genome.
      #
      # * *Returns* :
      #   -
      #
      def perform
        result_file = nil
        
        # Create the alignment files
        result_file = generate_alignment if @task == :all || @task == :align
        
        # Identify the clusters
        result_file = identify_clusters if @task == :all || @task == :cluster
        
        result_file
      end
      
      private
      
      # Uses the transcripts_file and genome_file to generate the alignment file.
      #
      # * *Returns* :
      #   - A File object handling the results of the alignment.
      #
      def generate_alignment
        raise ArgumentError, 'Missing genome FASTA file.' unless @genome_file
        raise ArgumentError, 'Missing transcripts FASTA file.' unless @transcripts_file
        
        # Prepare the BLAT alignment
        blat = Alignment::BLAT.new(@blat_options.merge({ out_format: :tab, database: @genome_file }))
        
        # Optionally set a permanent file to write the results to
        @alignment_file ||= "#{@transcripts_file}.alignment"
        blat.output_file = @alignment_file
        
        puts "Running BLAT alignment..." if @verbose
        
        # Run
        result_file = blat.run(@transcripts_file)
        result_file.path
      end
      
      # Parses through the alignment file and generates a gff3 based on the alignment clusters.
      #
      # * *Returns* :
      #   - A String specifying the path of the resultant file.
      #
      def identify_clusters
        raise ArgumentError, 'Missing BLAT alignment file.' unless @alignment_file
        
        # Pull out all of the hits clustered by the query into an array
        puts "Reading alignment file..." if @verbose
        query_hits = []
        IO::BLASTFormat.open(@alignment_file) do |f|
          f.each_query do |query, hits| 
            query_hits << hits.select { |hit| hit.e_value <= @e_value }
          end
        end
        
        # Create a status bar to monitor thread process
        pbar = ProgressBar.new("Converting Hits", query_hits.size, STDOUT) if @verbose
        
        # Thread the creation of the entries from the hits
        entries = Utilities::Threader.thread(query_hits, threads: @threads) do |thread_query_hits|
          thread_query_hits.map do |hits|
            pbar.inc if @verbose
          
            # Cluster the hits
            clusters = Alignment::Aligner.cluster_hits(hits, cluster_on: :subject)
    
            # Convert the clusters to entries
            clusters.map { |clustered_hits| create_feature(clustered_hits) }
          end
        end.flatten

        # Write the file
        puts "Writing entries to file..." if @verbose
        @output_file ||= "#{@alignment_file}.gff3"
        IO::GFFFormat.open(@output_file, mode: 'w') do |f|
          f.puts_header
          f.puts(entries, progress_bar: true, id_prefix: @id_prefix)
        end
        
        @output_file
      end
    
      # Takes an array of hits and creates a GFF entry out of them.
      #
      # * *Args*    :
      #   - +hit_clusters+ -> An array of BLAST::Hits to source the entry.
      # * *Returns* :
      #   - A GFF::Feature representing the collection of hits.
      #
      def create_feature(hits)
        # Initialize the entry
        query, subject = hits.first.query, hits.first.subject
        feature = IO::GFF::Feature.new(seqid: subject, source: 'BLATN', type: :nucleotide_match, attributes: { 'Name' => query })

        # Determine the orientation based on the query start/end positions
        feature.strand = hits.first.subject_start < hits.first.subject_end ? '+' : '-'

        # Add each of the hits to the entry
        hits.each do |hit|
          feature.regions.create(start: hit.subject_start, 
                               end: hit.subject_end, 
                               score: hit.bit_score, 
                               attributes: { 'EValue' => hit.e_value, 'Target' => "#{query} #{hit.query_start} #{hit.query_end}" })
        end
      
        feature
      end
    end
  end
end
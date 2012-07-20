module Genomics
  module Operation
    # This class handles the operation of a BLASTP alignment of a proteome against a reference database and constructs a GFF3 file
    # containing the alignment results.
    #
    # * *Attributes* :
    #   - +task+ -> A Symbol specifying the particular step to be conducted.  Possible values include :align, :cluster, 
    #   :all (Default: :all)
    #   - +alignment_file+ -> A String specifying the path the alignment file.
    #   - +verbose+ -> A boolean that sets whether the current operation is printed to STDOUT (Default: true). 
    #   - +threads+ -> An integer that sets the maximum number of threads that can be used (Default: 1).
    #
    # * *Alignment Attributes* :
    #   - +genome_file+ -> A String specifying the path to the FASTA file with the genome contig sequences.
    #   - +proteome_file+ -> A String specifying the path to the FASTA file with the proteome sequences. 
    #   - +blast_path+ -> A String explicitly specifying the path to BLASTX.
    #   - +blast_options+ -> A Hash of options to be passed directly to the Alignment::BLAST object used to conduct the alignment.
    #   - +e_value+ -> A float specifying the minimum EValue that a hit must have in order to be retained (Default: 1e^-5).
    #
    # * *Cluster Attributes* :
    #   - +output_file+ -> A String specifying the path of the generated GFF3 file.
    #   - +id_prefix+ -> A String that will be used as a prefix for the ID field in the resultant GFF3 file (Default: 'blastx_'). 
    class BLASTP
      
      attr_accessor :task, :verbose, :threads, :alignment_file,
                    :genome_file, :proteome_file, :blast_path, :blast_options,
                    :output_file, :id_prefix, :e_value 
      
      def initialize(options = {})
        options = { task: :all, verbose: true, threads: 1, 
                    blast_path: :blastx, blast_options: {}, e_value: 1e-5,
                    id_prefix: 'blastx_' }.merge(options)
        
        options.each do |name, value|
          send("#{name}=", value)
        end
      end
      
      # Performs the operation based on the attributes.
      #
      # * *Returns* :
      #   - A String with the path to the resultant file.
      #
      def perform
        result_file = nil
        
        # Create the alignment files
        # result_file = generate_alignment if @task == :all || @task == :align
        
        # Identify the clusters
        result_file = identify_clusters if @task == :all || @task == :cluster
        
        result_file
      end
      
      private
      
      # Parses through the alignment file and generates a gff3 based on the alignment clusters.
      #
      # * *Returns* :
      #   - A String specifying the path of the resultant file.
      #
      def identify_clusters
        raise ArgumentError, 'Missing BLASTP alignment file.' unless @alignment_file
        
        # Go through all of the hits by query generating 
        puts "Reading alignment file..." if @verbose
        query_hits = []
        IO::BLASTFormat.open(@alignment_file) do |f|
          f.each_query do |query, hits| 
            query_hits << hits
          end
        end
        
        # Create a status bar to monitor thread process
        pbar = ProgressBar.new("Converting Hits", query_hits.size, STDOUT) if @verbose
        
        # Thread the creation of the features from the hits
        features = Utilities::Threader.thread(query_hits, threads: @threads) do |thread_query_hits|
          thread_query_hits.map do |hits|
            pbar.inc if @verbose
            
            create_feature(hits)
          end
        end.flatten
        
        # Write the file
        puts "Writing entries to file..." if @verbose
        @output_file ||= "#{@alignment_file}.gff3"
        IO::GFFFormat.open(@output_file, mode: 'w') do |f|
          f.puts_header
          f.puts(features, progress_bar: true, id_prefix: @id_prefix)
        end
        
        @output_file
      end
      
      # Takes an array of hits and creates a GFF entry out of them.
      #
      # * *Args*    :
      #   - +hits+ -> An array of BLAST::Hits to source the entry.
      # * *Returns* :
      #   - A GFF::Feature representing the collection of hits.
      #
      def create_feature(hits)
        # Initialize the entry
        query, subject = hits.first.query, hits.first.subject

        # Determine the orientation based on the query start/end positions
        feature.strand = hits.first.query_start < hits.first.query_end ? '+' : '-'
        
        # Convert the hit to a feature
        
        # Add each of the hits to the entry
        hits.map do |hit|
          feature = IO::GFF::Feature.new(seqid: query, source: 'BLASTP', type: :match, attributes: { 'Name' => subject })
          feature.regions.create(start: hit.query_start, 
                                 end: hit.query_end, 
                                 score: hit.bit_score, 
                                 attributes: { 'EValue' => hit.e_value, 'Target' => "#{subject} #{hit.subject_start} #{hit.subject_end}" })
                               
          feature                               
        end
        
      
        feature
      end
    end
  end
end

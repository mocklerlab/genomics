require 'pathname'

module Genomics
  module Operation
    # This class handles the operation logic of conducting a reciprocal best BLAST analysis.  The alignments can be generated and the 
    # results can be analyzed to identify the putative orthologs.
    #
    # * *Attributes* :
    #   - +task+ -> A Symbol specifying the particular step to be conducted.  Possible values include :align, :ortholog, 
    #   :all (Default: :all)
    #   - +alignment_files+ -> A String specifying the path the alignment files.
    #   - +verbose+ -> A boolean that sets whether the current operation is printed to STDOUT (Default: true). 
    #
    # * *Alignment Attributes* :
    #   - +proteomes+ -> An Array of hashes specifying the files to be used.
    #   - +blast_path+ -> A String giving the path to blastp. (Default: blastp)
    #   - +blast_options+ -> A Hash of options to be passed directly to the Alignment::BLAST object used to conduct the alignment.
    #   - +alignment_file_dir+ -> A String specifying the directory where the alignment files should be 
    #       written.  (Default: system temporary directory)
    #   - +e_value+ -> A float specifying the minimum EValue that a hit must have in order to be retained (Default: 1e^-5).
    #
    # * *Orthologs Attributes* :
    #   - +output_file+ -> A String specifying the path of the generated file. (Default: orthologs.tab)
    #   - +detailed+ -> A boolean indicating whether or not include extra details of the alignment between the two. (Default: false)
    #   - +reciprocal+ -> A boolean indicating whether or not to only include reciprocal best alignments. (Default: true)
    #
    class RBB
      
      attr_accessor :task, :verbose, :threads,
                    :proteomes, :blast_path, :blast_options, :alignment_file_dir, :e_value,
                    :output_file, :detailed, :reciprocal
      attr_reader :alignment_files
      
      def initialize(options = {})
        options = { task: :all, verbose: true, threads: 1,
                    blast_path: :blastp, blast_options: { out_format: :tab }, e_value: 1e-5,
                    output_file: 'orthologs.tab', detailed: false, reciprocal: true }.merge(options)
        
        options.each do |name, value|
          send("#{name}=", value)
        end
      end
      
      # Sets the alignment files, which are an internal array of File objects, from the strings or file objects provided.
      #
      # * *Args*    :
      #   - +new_alignment_files+ -> An Array of objects that respond to path.
      # * *Returns* :
      #   - An Array of File objects
      #
      def alignment_files=(new_alignment_files)
        @alignment_files = new_alignment_files.map { |file| File.open(file) }
      end
      
      # This runs the entire process of generating orthologs from two proteomes.
      #
      # * *Returns* :
      #   - Returns an Integer indicating the number of orthologs identified.
      #
      def perform
        result = nil
        
        # Create the alignment files
        result = generate_alignments if @task == :all || @task == :align
        
        # Identify the orthologs
        result = identify_orthologs if @task == :all || @task == :ortholog
        
        # Clean up resources
        @alignment_files.map(&:close)
        
        result
      end
         
      private   
         
      # This takes a pair of FASTA formated proteomes and mutually aligns them.
      #
      # * *Returns* :
      #   - An Array of Files containing the alignment results 
      #
      def generate_alignments
        raise ArgumentError, 'Missing proteome files.' unless @proteomes && @proteomes.any?
        
        # Create the database/query pairs
        alignment_combinations = @proteomes.product(@proteomes)
        alignment_combinations = alignment_combinations.select { |combination| combination[0] != combination[1] }
        alignment_combinations.map! { |combination| [combination[0][:file], combination[1][:database]] }

        # Iteratively process the pairs
        @alignment_files = alignment_combinations.map do |combination|
          # Create and run the alignment
          blast = Alignment::BLAST.new(@blast_path, @blast_options.merge(threads: @threads, e_value: @e_value))
          blast.database = combination[1]
          
          # Optionally set a permanent file to write the results to
          if @alignment_file_dir
            filename = "#{Pathname.new(combination[0]).basename}_vs_#{Pathname.new(combination[1]).basename}"
            blast.output_file = File.join(@alignment_file_dir, filename) 
          end
          
          puts "Running BLASTP alignment..." if @verbose
          
          # Get the results and clean up the file resoursces
          result_file = blast.run(combination[0])
          result_file
        end
      end   
      
      # Reads the two files supplied and identifies the reciprocal best hits, writing them to an output file.
      #
      # * *Returns* :
      #   - An Integer with the number of alignments written to file
      #
      def identify_orthologs
        raise ArgumentError, 'Missing alignment files.' unless @alignment_files && @alignment_files.any?
        
        # Get the best alignments for each query from each file
        best_query_hits, best_target_hits = @alignment_files.map { |alignment_file| parse_best_alignments(alignment_file) }

        # Determine which best alignments are pairwise reciprocal
        pbar = ProgressBar.new("Calc. Reciprocal", best_query_hits.size, STDOUT) if @verbose
        
        reciprocal = []
        best_query_hits.each do |query, hits|
          pbar.inc if @verbose
          
          # Go through each hit and see if it is reciprocal
          hits.each do |hit|
            reverse_hits = best_target_hits[hit.subject] || []
            reciprocal << [query, hit.subject] if reverse_hits.find { |reverse_hit| reverse_hit.subject == query }
          end
        end
        
        # Create a hash to look up which hits are reciprocal
        reciprocal_hash = {}
        reciprocal.each do |pair|
          reciprocal_hash[pair[0]] ||= []
          reciprocal_hash[pair[0]] << pair[1]
        end
        
        # Iterate through the BLAST file generating results
        pbar = ProgressBar.new("Generating Results", best_query_hits.size, STDOUT) if @verbose
        
        results = []
        IO::BLASTFormat.open(alignment_files.first) do |f|
          f.each_query do |query, hits|
            pbar.inc if @verbose
            reciprocal_hits = reciprocal_hash[query] || []
            
            # Remove duplicate hits, retaining only the best
            unique_hits = {}
            hits.each do |hit|
              unique_hits[hit.subject] = hit if unique_hits[hit.subject].nil? || unique_hits[hit.subject].bit_score < hit.bit_score
            end
            
            # Generate the results
            unique_hits.values.each do |hit|
              is_reciprocal = reciprocal_hits.include?(hit.subject)
              next if @reciprocal && !is_reciprocal
              
              # Generate the data to include in the results
              hit_data = [query, hit.subject]
              hit_data << is_reciprocal unless @reciprocal
              hit_data += [hit.e_value, hit.bit_score] if @detailed
              
              results << hit_data
            end
          end
        end
        
        # Write the results to file
        File.open(@output_file, 'w') do |f|
          f.puts results.sort.map { |result| result.join("\t") }
        end

        results.count
      end
            
      # Parses the alignment file specified and returns a hash with the queries as keys mapped to the best hits.
      #
      # * *Args* :
      #   - +alignment_file+ -> A BLASTFormat IO object that will be used for parsing.
      # * *Returns* :
      #   - A Hash of query_ids to hits
      #
      def parse_best_alignments(alignment_file)
        # Reduce the hits to a hash or reciprocally best matches
        alignments = {}
        IO::BLASTFormat.open(alignment_file) do |f|
          f.each do |hit|
            case
            when alignments[hit.query].nil? || alignments[hit.query].first.bit_score < hit.bit_score
              alignments[hit.query] = [hit]
            when alignments[hit.query].first.bit_score == hit.bit_score
              alignments[hit.query] << hit
            end
          end
        end
        
        alignments
      end
    end
  end
end
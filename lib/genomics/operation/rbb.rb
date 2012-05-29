require 'pathname'

module Genomics
  module Operation
    # This class handles the operation logic of conducting a reciprocal best BLAST analysis.  The alignments can be generated and the 
    # results can be compared to identify the putative orthologs.
    class RBB
      class << self
        
        # This runs the entire process of generating orthologs from two proteomes.
        #
        # * *Args*    :
        #   - +proteomes+ -> An Array of hashes specifying the files to be used.
        # * *Returns* :
        #   - Returns the number of orthologs identified.
        #
        def perform(proteomes, options = {})
          # Create the alignment files
          alignment_files = generate_alignments(proteomes, options)

          # Pass the alignment files onto the ortholog identification
          identify_orthologs(alignment_files)
        end
        
        # This takes a pair of FASTA formated proteomes and mutually aligns them.
        #
        # * *Args*    :
        #   - +proteomes+ -> An Array of hashes specifying the files to be used.
        # * *Options* :
        #   - +blast_path+ -> A String giving the path to blastp.
        #   - +blast_options+ -> A Hash of options to be passed directly to the Alignment::BLAST object used to conduct the alignment.
        #   - +alignment_file_dir+ -> A String specifying the directory where the alignment files should be written.  (Default: temporary directory)
        # * *Returns* :
        #   - An Array of file paths containing the alignment results 
        #
        def generate_alignments(proteomes, options = {})
          options = { blast_path: :blastp, blast_options: {} }.merge(options)
          options[:blast_options] = { e_value: 0.00001 }.merge(options[:blast_options])
          options[:blast_options][:out_format] = :tab
          
          # Create the database/query pairs
          alignment_combinations = proteomes.product(proteomes)
          alignment_combinations = alignment_combinations.select { |combination| combination[0] != combination[1] }
          alignment_combinations.map! { |combination| [combination[0][:file], combination[1][:database]] }

          # Break up the processing into threads
          # threads = alignment_combinations.map do |combination|
          #   Thread.new do
          alignment_combinations.map do |combination|
            # Create and run the alignment
            blast = Alignment::BLAST.new(options[:blast_path], options[:blast_options])
            blast.database = combination[1]
            
            # Optionally set a permanent file to write the results to
            if options[:alignment_file_dir]
              filename = "#{Pathname.new(combination[0]).basename}_vs_#{Pathname.new(combination[1]).basename}"
              blast.output_file = File.join(options[:alignment_file_dir], filename) 
            end
            
            # Get the results and clean up the file resoursces
            result_file = blast.run(combination[0])
            result_file.path
          end
          #   end
          # end
          
          # Collect the threads and return
          # threads.map { |thread| thread.value.path }
        end
        
        # Reads the two files supplied and identifies the reciprocal best hits, writing them to an output file.
        #
        def identify_orthologs(alignment_files, options = {})
          options = { output_file: 'orthologs.tab' }.merge(options)
          
          # Get the best alignments for each query from each file
          alignments1, alignments2 = alignment_files.map { |alignment_file| parse_alignments(alignment_file) }
          
          # Determine which best alignments are pairwise reciprocal
          pbar = ProgressBar.new("Calc. Reciprocal", alignments1.size, STDOUT)

          reciprocal = []
          alignments1.each do |query, hits|
            pbar.inc

            hits.each do |hit|
              reverse_hits = alignments2[hit.subject] || []
              reciprocal << [query, hit.subject] if reverse_hits.find { |hit2| hit2.subject == query }
            end
          end

          # Write the results to file
          File.open(options[:output_file], 'w') do |f|
            reciprocal.sort.each { |pair| f.puts pair.join("\t") }
          end

          reciprocal.count
        end
    
        private
    
        # Parses the alignment file specified and returns a hash with the queries as keys mapped to the row.
        # 
        def parse_alignments(alignment_file)
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
end
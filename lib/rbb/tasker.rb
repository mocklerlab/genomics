module Genomics
  module RBB
    class Tasker
      extend DNA::Utility::FileParser
      
      class << self
        # Reads the two files supplied and identifies the reciprocal best hits, writing them to an output file.
        #
        def identify(alignment_file1, alignment_file2, options = {})
          options = { output_file: 'orthologs.tab' }.merge(options)
          
          # Get the best alignments for each query from each file
          alignments1 = parse_alignments(alignment_file1)
          alignments2 = parse_alignments(alignment_file2)
          
          # Determine which best alignments are pairwise reciprocal
          pbar = ProgressBar.new("Calc. Reciprocal", alignments1.size, STDOUT)
          
          reciprocal = []
          alignments1.each do |query, alignment_hash|
            pbar.inc
            
            alignment_hash[:subjects].each do |subject|
              reciprocal << [query, subject] if alignments2[subject][:subjects].include?(query)
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
          alignments = {}
          
          # Create the progress bar
          `grep -v '#' #{alignment_file} | wc -l` =~ /\d+/
          pbar = ProgressBar.new("Parsing Alignments", $~[0].to_i, STDOUT)
          
          parse_file(alignment_file) do |row|
            pbar.inc
            
            # Retain the alignment pair(s) with the highest bit score
            bit_score = row[11].to_f
            case
            when alignments[row[0]].nil? || alignments[row[0]][:bit_score] < bit_score
              alignments[row[0]] = { subjects: [row[1]], bit_score: bit_score }
            when alignments[row[0]][:bit_score] == bit_score
              alignments[row[0]][:subjects] << row[1]
            end
          end
          
          alignments
        end
      end
    end
  end
end
module Genomics
  module Alignment
    class RBB
      class << self
        # Reads the two files supplied and identifies the reciprocal best hits, writing them to an output file.
        #
        def identify(alignment_files, options = {})
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
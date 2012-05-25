module Genomics
  module Alignment
    # This class handles conducting an alignment of ESTs onto a genome.
    class EST < Aligner
      class << self
        # Reads in the results from the provided file and transforms them into the specified format.
        #
        # * *Args*    :
        #   - +filename+ -> A string determining the file to be read and transformed into another filetype.
        # * *Returns* :
        #   - A boolean describing whether or not the transformation was successful.
        #
        def transform(filename, options = {})
          options = { format: :gff3, output_file: "#{filename}.gff3" }.merge(options)
          
          # Parse the file to get all of the alignment hits
          aggregated_hits = []
          IO::BLASTFormat.open(filename) do |f|
            aggregated_hits = f.entries(aggregate_hits: true)
          end
          
          # Perform the clustering.
          clustered_hits = cluster_aggregated_hits(aggregated_hits)
          
          # Generate the entries
          pbar = ProgressBar.new("Generating Entries", clustered_hits.size, STDOUT)
          
          entries = []
          clustered_hits.each do |hits|
            pbar.inc
            
            # Initialize the entry
            query, subject = hits.first.query, hits.first.subject
            entry = IO::GFF::Entry.new(seqid: subject, source: 'BLATN', type: :EST_match, attributes: { 'Name' => query })
            
            # Determine the orientation based on the query start/end positions
            entry.strand = hits.first.subject_start < hits.first.subject_end ? '+' : '-'

            # Add each of the hits to the entry
            hits.each do |hit|
              entry.regions.create(start: hit.subject_start, 
                                   end: hit.subject_end, 
                                   score: hit.bit_score, 
                                   attributes: { 'EValue' => hit.e_value, 'Target' => "#{query} #{hit.query_start} #{hit.query_end}" })
            end
            
            entries << entry
          end

          # Write the file
          puts "Writing entries to file..."
          IO::GFFFormat.open(options[:output_file], 'w') do |f|
            f.puts_header
            f.puts(entries)
          end
          
          true
        end
      end
    end
  end
end
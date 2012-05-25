module Genomics
  module IO
    # This class acts as an IO allowing tab delimited BLAST files to be read and written.
    # TODO: Expand with more details
    class BLASTFormat < FlatFileFormat
      
      # Iterates through each of the Hits in the file yield each successively to the block. 
      #
      def each
        super do |row|
          # Pull out the attribute values from the row
          attributes = { query:                row[0],
                         subject:              row[1],
                         query_start:          row[6].to_i,
                         query_end:            row[7].to_i,
                         subject_start:        row[8].to_i,
                         subject_end:          row[9].to_i,
                         e_value:              row[10], 
                         bit_score:            row[11].to_f,
                         percentage_identity:  row[2].to_f,
                         alignment_length:     row[3].to_i,
                         mismatches:           row[4].to_i,
                         gap_openings:         row[5].to_i }

          # Create the new Hit object and yield it to the block for processing
          yield BLAST::Hit.new(attributes)
        end
      end
      
      # Iterates through the file by the queries rather than by individual hits.  The block is successively yielded the string
      # name of the query and a list of the hits for that query.
      #
      def each_query
        # Set the state variables
        current_rows = []
        current_query = nil
        
        # Read through rows until the query changes.
        each do |hit|
          if hit.query == current_query
            current_rows << hit
          else
            yield current_query, current_rows if current_query
            current_query, current_rows = hit.query, [hit]
          end
        end
      end
      
      # Retrns a collection of all of the entries in the file.
      #
      # * *Options*    :
      #   - +aggregate_hits+ -> A boolean specifying whether the alignment hits in the file should be returned as a multi-dimmensional
      #                         datastructure.  If true, a hash of hashes is returned with the first keys being the query ids, and the
      #                         the second set of keys being the subject ids.
      # * *Returns* :
      #   - A collection of the individual BLAST::Hit objects parsed in the file.
      #
      def entries(options = {})
        options = { aggregate_hits: false }.merge(options)
        
        # Get the hits
        hits = []
        each { |hit| hits << hit }
        
        return hits unless options[:aggregate_hits]

        # Group the hits based on the specifics of what was matched.
        aggregated_hits = {}
        hits.each do |hit|
          aggregated_hits[hit.query] ||= {}
          aggregated_hits[hit.query][hit.subject] ||= []
          aggregated_hits[hit.query][hit.subject] << hit
        end
        
        aggregated_hits
      end
      
      # Writes the entries to the IO stream.
      #
      # * *Args*    :
      #   - +objects+ -> The Genomic::IO:GFF:Entry objects to be written successively to the stream.
      #
      def puts(*entries)
        entries = entries.first if entries.length == 1 && entries.first.is_a?(Array)
        
        @last_id ||= 0
        entries.sort.each do |entry|
          entry.attributes["ID"] ||= (@last_id += 1)
          @io.puts entry.to_gff
        end
      end
      
      # Write the pragma for a valid gff3 header.
      #
      def puts_header
        @io.puts '##gff-version 3'
      end
    end
  end
end
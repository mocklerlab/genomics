module Genomics
  module IO
    # This class acts as an IO allowing tab delimited BLAST files to be read and written.
    # TODO: Expand with more details
    class BLASTFormat < FlatFileFormat
      
      # Iterates through each of the Hits in the file yield each successively to the block. 
      #
      def each
        unless @format == :xml
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
                           identities:           (row[2].to_f * row[3].to_i / 100).round,
                           alignment_length:     row[3].to_i,
                           gap_openings:         row[5].to_i 
                         }

            # Create the new Hit object and yield it to the block for processing
            yield BLAST::Hit.new(attributes)
          end
        else
          # Load the entire file into memory and parse the xml structure yielding the individual hits
          root_element = Ox.load_file(@io.path).root
          root_element.locate('BlastOutput_iterations/Iteration').each do |iteration_node|
            # Pull out the query
            query = iteration_node.locate('Iteration_query-def/?').first
            
            # Start pulling out individual hits
            iteration_node.locate('Iteration_hits/Hit').each do |hit_node|
              subject = hit_node.locate('Hit_def/?').first
              
              # Create the individual hits from the HSPs
              hit_node.locate('Hit_hsps/Hsp').each do |hsp|
                # Pick out the individual HSP details
                attributes = { query: query, 
                               subject: subject,
                               bit_score:         hsp.locate('Hsp_bit-score/?').first.to_f,
                               e_value:           hsp.locate('Hsp_evalue/?').first,
                               query_start:       hsp.locate('Hsp_query-from/?').first.to_i,
                               query_end:         hsp.locate('Hsp_query-to/?').first.to_i,
                               subject_start:     hsp.locate('Hsp_hit-from/?').first.to_i,
                               subject_end:       hsp.locate('Hsp_hit-to/?').first.to_i,
                               query_frame:       hsp.locate('Hsp_query-frame/?').first.to_i,
                               subject_frame:     hsp.locate('Hsp_hit-frame/?').first.to_i,
                               identities:        hsp.locate('Hsp_identity/?').first.to_i,
                               positives:         hsp.locate('Hsp_positive/?').first.to_i,
                               gap_openings:      hsp.locate('Hsp_gaps/?').first.to_i,
                               alignment_length:  hsp.locate('Hsp_align-len/?').first.to_i,
                               query_sequence:    hsp.locate('Hsp_qseq/?').first,
                               subject_sequence:  hsp.locate('Hsp_hseq/?').first,
                               midline:           hsp.locate('Hsp_midline/?').first,
                             }
                             
                # Create the new Hit object and yield it to the block for processing
                yield BLAST::Hit.new(attributes)
              end
            end
            
          end
        end
      end
      
      # Iterates through the file by the queries rather than by individual hits.  The block is successively yielded the string
      # name of the query and a list of the hits for that query.
      #
      # * *Arguments*   :
      #   - +block+ -> A block that is passed the name of a query (String) and an Array of the Hits for the query.
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
      #   - +aggregate+ -> A boolean specifying whether the alignment hits in the file should be returned as a multi-dimmensional
      #                         datastructure.  If true, a hash of hashes is returned with the first keys being the query ids, and the
      #                         the second set of keys being the subject ids.
      #   - +sort+ -> A boolean specifying whether or not to sort the hits (Default false).
      #   - +transpose+ -> A boolean specifying whether or not to switch the query and subject values on the hit (Default false).
      # * *Returns* :
      #   - A collection of the individual BLAST::Hit objects parsed in the file or a multi-dimensional Hash.
      #
      def entries(options = {})
        options = { aggregate: false, sort: false, transpose: false }.merge(options)
        
        # Get the hits
        hits = []
        each { |hit| hits << hit }
        
        # Transpose the hits if selected
        hits.map!(&:transpose!) if options[:transpose]
        
        # Sort the hits if selected
        hits.sort! if options[:sort]

        return hits unless options[:aggregate]

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
module Genomics
  module IO
    # This class acts as an IO allowing tab delimited BLAST files to be read and written.
    # TODO: Expand with more details
    class BLASTFormat < FlatFileFormat
      
      # Iterates through each of the Hits in the file yielding each successively to the block. 
      #
      # * *Options*
      #   - +sort+ -> A boolean specifying whether or not to sort the hits for each query (Default false).
      #   - +subject_regex+ -> A RegEx used to reduce the subject name of the hit to a substring of the value in the file.  This is
      #   useful for parsing out extraneous information.
      #
      def each(options = {})
        options = { subject_regex: nil }.merge(options)
        
        unless @format == :xml
          super do |row|
            # Parse out the name of the subject if instructed
            subject = row[1]
            subject = $~.to_s if options[:subject_regex] && subject.match(options[:subject_regex])
            
            # Pull out the attribute values from the row
            attributes = { query:                row[0],
                           subject:              subject,
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

              # Parse out the name of the subject if instructed
              subject = $~.to_s if options[:subject_regex] && subject.match(options[:subject_regex])

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
      #   - +block+ -> A block th at is passed the name of a query (String) and an Array of the Hits for the query.
      # * *Options*
      #   - +sort+ -> A boolean specifying whether or not to sort the hits for each query (Default false).
      #   - +subject_regex+ -> A RegEx used to reduce the subject name of the hit to a substring of the value in the file.  This is
      #   useful for parsing out extraneous information.
      #
      def each_query(options = {})
        options = { sort: false, subject_regex: nil }.merge(options)
        
        # Set the state variables
        current_rows = []
        current_query = nil

        # Read through rows until the query changes.
        each(subject_regex: options[:subject_regex]) do |hit|
          if hit.query == current_query
            current_rows << hit
          else
            if current_query
              current_rows.sort! if options[:sort]
              yield current_query, current_rows 
            end

            current_query, current_rows = hit.query, [hit]
          end
        end
        
        # Yield the last query if any were found
        yield current_query, current_rows if current_query
      end
      
      # Iterates through the file by clusters of hits.  Hits are clustered together if they are proximal on the target.  This
      # determination is made while paying respect to orientation and to the increasing positions of individual hits and how they
      # align as a group. The block is successively yielded the string name of the query, the name of the subject, and the sorted 
      # Array of hits.
      #
      # * *Arguments*   :
      #   - +block+ -> A block th at is passed the name of a query (String) and an Array of the Hits for the query.
      # * *Options*
      #   - +cluster_on+ -> A symbol determining whether to cluster hits based on positions on the :query or :subject sequences.
      # * *Options*
      #   - +subject_regex+ -> A RegEx used to reduce the subject name of the hit to a substring of the value in the file.  This is
      #   useful for parsing out extraneous information.
      #
      def each_cluster(options = {})
        options = { cluster_on: :subject, subject_regex: nil }.merge(options)
        
        each_query(subject_regex: options[:subject_regex]) do |query, hits|
          # Cluster the hits
          clusters = Alignment::Aligner.cluster_hits(hits, cluster_on: options[:cluster_on])

          # Yeild the results
          clusters.each do |cluster|
            yield query, cluster.first.subject, cluster
          end
        end
      end
      
      # Returns an array of all of the hits in the file.
      #
      # * *Args*    :
      #   - +sort+ -> A boolean specifying whether or not to sort the hits (Default false).
      #   - +transpose+ -> A boolean specifying whether or not to switch the query and subject values on the hit (Default false).
      #   - +subject_regex+ -> A RegEx used to reduce the subject name of the hit to a substring of the value in the file.  This is
      #   useful for parsing out extraneous information.
      # * *Returns* :
      #   - An Array
      #
      def hits(options = {})
        options = { sort: false, transpose: false, subject_regex: nil }.merge(options)
        
        # Get the hits
        hits = []
        each(options[:subject_regex]) { |hit| hits << hit }
        
        # Transpose the hits if selected
        hits.map!(&:transpose!) if options[:transpose]
        
        # Sort the hits if selected
        hits.sort! if options[:sort]
        
        hits
      end
      
      # Returns a hash of queries paired with their matching hits.
      #
      # * *Args*    :
      #   - +sort+ -> A boolean specifying whether or not to sort the hits (Default false).
      #   - +transpose+ -> A boolean specifying whether or not to switch the query and subject values on the hit (Default false).
      #   - +subject_regex+ -> A RegEx used to reduce the subject name of the hit to a substring of the value in the file.  This is
      #   useful for parsing out extraneous information.
      # * *Returns* :
      #   - An Array
      #
      def query_hits(options = {})
        options = { sort: false, transpose: false, subject_regex: nil }.merge(options)
        
        # Get the hits
        hits_hash = {}
        each_query(sort: options[:sort], subject_regex: options[:subject_regex]) do |query, hits| 
          # hits.map!(&:transpose!) if options[:transpose]
          hits_hash[query] = hits
        end

        # TODO: Transpose.  This needs to actually interchange subject and query on the hash level
        
        hits_hash
      end
      
      # Returns a multi-dimmensional hash of queries matched to subjects and their respective clusters of hits.
      #
      # * *Args*    :
      #   - +cluster_on+ -> A symbol determining whether to cluster hits based on positions on the :query or :subject sequences.
      #   - +transpose+ -> A boolean specifying whether or not to switch the query and subject values on the hit (Default false).
      #   - +subject_regex+ -> A RegEx used to reduce the subject name of the hit to a substring of the value in the file.  This is
      #   useful for parsing out extraneous information.
      # * *Returns* :
      #   - An Array
      #
      def clustered_hits(options = {})
        options = { cluster_on: :subject, transpose: false, sort: false, subject_regex: nil }.merge(options)
        
        # Get the hits
        hits_hash = {}
        each_cluster(cluster_on: options[:cluster_on], subject_regex: options[:subject_regex]) do |query, subject, hits| 
          hits.map!(&:transpose!) if options[:transpose]
          hits_hash[query] ||= {}
          hits_hash[query][subject] ||= []
          hits_hash[query][subject] << hits
        end
        
        if options[:sort]
          sorted_hits_hash = {}
          # For each query go through and sort the subjects and their respective clusters
          hits_hash.each do |query, subject_hash|
            sorted_hits_hash[query] ||= {}
            
            sorted_subjects = subject_hash.sort_by do |subject, clusters|
              -clusters.map { |hits| hits.map(&:bit_score).max }.max
            end

            sorted_subjects.each do |subject, clusters|
              sorted_hits_hash[query][subject] ||= []
              
              # Sort the clusters, and sort the hits within the cluster
              clusters.sort_by! { |hits| -hits.sort!.max.bit_score }
              sorted_hits_hash[query][subject] = clusters
            end
          end
          
          hits_hash = sorted_hits_hash
        end
        # TODO: Transpose.  This needs to actually interchange subject and query on the hash level
        
        hits_hash
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
      # #TODO Deprecated
      # def entries(options = {})
      #   options = { aggregate: false, sort: false, transpose: false }.merge(options)
      #   
      #   # Get the hits
      #   hits = []
      #   each { |hit| hits << hit }
      #   
      #   # Transpose the hits if selected
      #   hits.map!(&:transpose!) if options[:transpose]
      #   
      #   # Sort the hits if selected
      #   hits.sort! if options[:sort]
      # 
      #   return hits unless options[:aggregate]
      # 
      #   # Group the hits based on the specifics of what was matched.
      #   aggregated_hits = {}
      #   hits.each do |hit|
      #     aggregated_hits[hit.query] ||= {}
      #     aggregated_hits[hit.query][hit.subject] ||= []
      #     aggregated_hits[hit.query][hit.subject] << hit
      #   end
      #   
      #   aggregated_hits
      # end
      
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
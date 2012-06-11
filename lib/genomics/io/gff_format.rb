module Genomics
  module IO
    # This class acts as an IO allowing GFF3 files to be read and written.
    class GFFFormat < FlatFileFormat
      
      # Iterates through each of the Features in the file yielding each successively to the block. 
      #
      def each(options = {})
        # Initialize the previously parsed data
        current_feature_hierarchy = []
        
        # Iterate through the rows
        super do |row|
          parsed_attributes = GFF::Feature.parse_attributes(row[8])
          id = parsed_attributes[:ID]
          
          # Handle the row depending on the state of iteration
          last_feature = nil
          while current_feature_hierarchy.any?
            # Get the feature at the bottom and check to see if matches the present row
            current_feature = current_feature_hierarchy.pop
            next unless current_feature.id

            if current_feature.id == (parsed_attributes[:Parent].is_a?(Array) ? parsed_attributes[:Parent].first : parsed_attributes[:Parent])
              last_feature = current_feature.features.create(parse_attributes(row))
              
              # Add both onto the hierarchy
              current_feature_hierarchy << current_feature << last_feature
              break
            elsif current_feature.id == id
              # Re-append the feature to the hierarchy
              current_feature_hierarchy << last_feature = current_feature
              break
            elsif parsed_attributes[:Derives_from] && current_feature.id == parsed_attributes[:Derives_from]
              # The feature is a special derivative instance
              last_feature = current_feature.derivatives.create(parse_attributes(row))
              
              # The derivative is an offshoot, only add the feature back onto the hierarchy
              current_feature_hierarchy << current_feature
              break
            elsif current_feature_hierarchy.empty?
              # None of the nodes were matched, so the current feature has been exhausted
              yield current_feature
            end
          end
          
          if current_feature_hierarchy.empty?
            ### Instantiate the next Feature ###
          
            # Pull out the attribute and create the object
            current_feature_hierarchy << last_feature = GFF::Feature.new(parse_attributes(row))          
          end
          
          # Add the region to the feature created by this row
          region_attributes = { start: row[3], end: row[4], score: row[5], phase: row[7], attributes: parsed_attributes }
          last_feature.regions.create(region_attributes)
          last_feature = nil
        end
        
        # Yield the last feature
        yield current_feature_hierarchy.first if current_feature_hierarchy.any?
      end
      
      # Writes the entries to the IO stream.  In order to maintain a valid gff3 file, it ensures that unique IDs are assigned to
      # each entry if they do not already exist.  For convenience, the entries are automatically sorted into the standard GFF3
      # order by increasing seqid and then start position along it.
      #
      # * *Args*    :
      #   - +objects+ -> The Genomic::IO:GFF:Entry objects to be written successively to the stream.
      #
      def puts(entries, options = {})
        options = { progress_bar: false, id_prefix: '' }.merge(options)
        entries = [entries] unless entries.is_a?(Array)
        
        # Create the progress bar
        pbar = ProgressBar.new("Writing #{Pathname.new(@io.path).basename}", entries.size, STDOUT) if options[:progress_bar]
        
        @last_id ||= 0
        entries.sort.each do |entry|
          pbar.inc if options[:progress_bar]
          
          entry.attributes["ID"] ||= "#{options[:id_prefix]}#{@last_id += 1}"
          @io.puts entry.to_gff
        end
      end
      
      # Write the pragma for a valid gff3 header.
      #
      def puts_header
        @io.puts '##gff-version 3'
      end
      
      private
      
      # Returns an attributes hash derived from the row.
      #
      def parse_attributes(row)
        { seqid:       row[0],
          source:      row[1],
          type:        row[2],
          strand:      row[6],
          attributes:  row[8] }
      end
      
    end
  end
end
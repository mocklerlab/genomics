module Genomics
  module IO
    # This class acts as an IO allowing GFF3 files to be read and written.
    class GFFFormat < FlatFileFormat
      GFF3_ID_REGEX = /ID=([A-Za-z0-9:\-_.]+);?/
      GFF3_NAME_REGEX = /Name=([A-Za-z0-9:\-_.]+);?/
      GFF3_PARENT_REGEX = /Parent=([A-Za-z0-9:\-_.]+);?/
      
      # Iterates through each of the Features in the file yielding each successively to the block. 
      #
      def each(options = {})
        # Initialize the previously parsed data
        current_feature = nil
        
        # Iterate through the rows
        super do |row|
          # Get the ID
          id = if row[8].match(GFF3_ID_REGEX)
            $~[1]
          else
            warn "Missing valid GFF3 ID: #{row.join("\t")}"
            nil
          end
          
          # Get the region specific attributes
          region_attributes = { start: row[3], end: row[4], score: row[5], phase: row[7], attributes: row[8] }
          
          # Check to see if the previous row has the same ID
          if current_feature && current_feature.id == id
            # Add the row to the feature
            # current_feature.regions.create(region_attributes)
          else
            # The feature has been completely read, so pass it to the block
            yield current_feature if current_feature
            
            ### Instantiate the next Feature ###
            
            # Pull out the attribute and create the object
            attributes = { seqid:       row[0],
                           source:      row[1],
                           type:        row[2],
                           strand:      row[6],
                           attributes:  row[8]
                         }
            current_feature = GFF::Feature.new(attributes)
          end
          
          # Add the region
          current_feature.regions.create(region_attributes)
        end
        
        # Yield the last feature
        yield current_feature if current_feature
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
    end
  end
end
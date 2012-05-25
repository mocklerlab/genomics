module Genomics
  module IO
    # This class acts as an IO allowing GFF3 files to be read and written.
    class GFFFormat < FlatFileFormat
      GFF3_ID_REGEX = /ID=([A-Za-z0-9:\-_.]+);?/
      GFF3_NAME_REGEX = /Name=([A-Za-z0-9:\-_.]+);?/
      GFF3_PARENT_REGEX = /Parent=([A-Za-z0-9:\-_.]+);?/
      
      # Writes the entries to the IO stream.  In order to maintain a valid gff3 file, it ensures that unique IDs are assigned to
      # each entry if they do not already exist.  For convenience, the entries are automatically sorted into the standard GFF3
      # order by increasing seqid and then start position along it.
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
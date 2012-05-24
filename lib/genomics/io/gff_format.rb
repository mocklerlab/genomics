module Genomics
  module IO
    # This class acts as an IO allowing GFF3 files to be read and written.
    class GFFFormat
      GFF3_ID_REGEX = /ID=([A-Za-z0-9:\-_.]+);?/
      GFF3_NAME_REGEX = /Name=([A-Za-z0-9:\-_.]+);?/
      GFF3_PARENT_REGEX = /Parent=([A-Za-z0-9:\-_.]+);?/
      
      class << self
        def open(filename, *args)
          if block_given?
            # Open the file using the block analogous to opening a standard file.
            File.open(filename, *args) do |f|
              yield self.new(f)
            end
          else
            # Return a new instance of the IO object, which must be explicitly closed
            self.new(File.open(filename, *args))
          end
        end
      end
      
      # This creates a Genomic::IO::GFFFormat object that wraps around a standard IO stream object.
      # This endows the IO stream with additional methods for accessing the records inside.
      #
      # * *Args*    :
      #   - +io+ -> An the IO object that is being read/written to. 
      #
      def initialize(io)
        @io = io
      end
      
      # Closes the IO object.
      #
      def close
        @io.close
      end
      
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
      
      # Rewinds the IO object if possible
      #
      def rewind
        @io.rewind
      end
      
      # Removes the file if possible.
      #
      def unlink
        @io.unlink
      end
    end
  end
end
module Genomics
  module IO
    # This class wraps a File IO and provides accessor methods for interacting with it.  It is intended to be subclassed
    # by classes that implement specific functionality for dealing with a particular file type.
    class FlatFileFormat
      
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
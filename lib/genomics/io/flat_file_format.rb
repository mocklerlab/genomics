require 'pathname'

module Genomics
  module IO
    # This class wraps a File IO and provides accessor methods for interacting with it.  It is intended to be subclassed
    # by classes that implement specific functionality for dealing with a particular file type.
    class FlatFileFormat
      include Enumerable
      
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
      
      # Executes the block for every line in the IO stream. Extra sanitation is applied to the contents of the file
      # to avoid encoding errors.
      #
      # * *Options*  :
      #   - +delimiter+ -> A string delimiter used to split the line into row values (Default: "\t").
      #   - +progress_bar+ -> A boolean specifying if the progress bar should be displayed (Default: true).
      #
      # * *Returns*  :
      #   - An Array to the block split on the optional delimiter provided.
      #
      def each(options = {})
        options = { delimiter: "\t", progress_bar: true }.merge(options)
        
        # Create the progress bar
        if options[:progress_bar]
          `wc -l #{@io.path}` =~ /\d+/
          pbar = ProgressBar.new("Parsing #{Pathname.new(@io.path).basename}", $~[0].to_i, STDOUT)
        end
        
        @io.each do |line|
          pbar.inc if options[:progress_bar]
          
          # Skip lines that are comments.
          # TODO: This eventually may want to be handled differently since some comments are actually pragmas
          begin
            next if line =~ /^#|>/
          rescue
            # Try to handle for bad character conversion
            line = Iconv.new('UTF-8//IGNORE', 'UTF-8').iconv(line)
            next if line =~ /^#|>/
          end
  
          yield line.split(options[:delimiter]).map(&:strip)
        end
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
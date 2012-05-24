module Genomics
  # This class implements a method for parsing through files and handing character conversion errors.
  class FileParser
    
    class << self
      # Parses the file line by line sending each line to the supplied block for processing.
      #
      def parse_file(file, delimiter = "\t", &action)
        raise "The file #{file} cannot be found" unless File.exists?(file)
  
        File.open(file) do |f|
          f.each_line do |line|
            # Skip lines that are comments.
            begin
              next if line =~ /^#|>/
            rescue
              # Try to handle for bad character conversion
              line = Iconv.new('UTF-8//IGNORE', 'UTF-8').iconv(line)
              next if line =~ /^#|>/
            end
    
            action.call(line.split(delimiter).map(&:strip))
          end
        end
      end
    end
  end
end

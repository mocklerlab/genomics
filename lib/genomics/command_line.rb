require 'open3'

module Genomics
  # This class provides methods to interface with programs on the command line.
  class CommandLine
    
    class << self
      
      # Allows commands to be run on the command line with the output from the command line being returned
      # via an IO stream.
      #
      # * *Args*    :
      #   - +command+ -> A string specifying the command line program to be run.
      #   - +block+ -> An optional block can be provided which will be passed the IO stream containing the STDOUT
      #   output of the command being run.
      # * *Returns* :
      #   - IO or boolean
      #
      def run(command)
        if block_given?
          Open3.popen3(command) do |stdin, stdout, stderror, wait_thread|
            # Raise the errors
            errors = stderror.read
            raise errors if errors != ""

            # Otherwise yield the output
            yield stdout
          end
          
          true
        else
          ::IO.popen(command)
        end
      end
    end
  end
end
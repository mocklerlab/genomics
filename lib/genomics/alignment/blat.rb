require 'tempfile'

module Genomics
  module Alignment
    # This object handles creating BLAT queries and sending them to the command line to be run.
    #
    # * *Options* :
    #
    class BLAT
      # A hash matching the optional attributes to the corresponding flags for command line execution
      OPTION_MAPPING = { header:                  ->(value) { "-noHead" unless value },
                         query_type:              ->(value) { "-q=#{value}"},
                         target_type:             ->(value) { "-t=#{value}"},
                         extend_through_large_n:  ->(value) { "-extendThroughN" if value }, 
                         out_format:              ->(value) {
                           numberCode = case value
                           when :tab then 'blast9'
                           else value
                           end
                           
                           "-out=#{numberCode}"
                         },
                         database_mask:           ->(value) { "-mask=#{value}" },
                         query_mask:              ->(value) { "-qMask=#{value}" },
                         minimum_match:           ->(value) { "-minMatch=#{value}" },
                         minimum_score:           ->(value) { "-minScore=#{value}" },
                         minimum_identity:        ->(value) { "-minIdentity=#{value}" },
                         maximum_gap:             ->(value) { "-maxGap=#{value}" },
                         tile_size:               ->(value) { "-tileSize=#{value}" },
                         step_size:               ->(value) { "-stepSize=#{value}" },
                         maximum_intron_size:     ->(value) { "-maxIntron=#{value}" }
                        }
      
      attr_accessor :database, :query_type, :target_type, :out_format, :e_value, :header, :extend_through_large_n, :database_mask, 
                    :query_mask, :minimum_match, :minimum_score, :minimum_identity, :maximum_gap, :tile_size, :step_size, 
                    :maximum_intron_size
      # :threads
      
      def initialize(options = {})
        options.each do |name, value|
          send("#{name}=", value)
        end
      end

      # Returns the previously created File object or creates a new Tempfile to store the output results.
      #
      # * *Returns* :
      #   - A File or Tempfile.
      #
      def output_file
        @output_file ||= Tempfile.new("blat_out")
      end
      
      # Uses the supplied File or path string to create the output file object, where the result from the alignment will be
      # written. 
      #
      # * *Args*    :
      #   - +new_ouput_file+ -> A File or String specifying the location of the file.
      # * *Returns* :
      #   - The newly created File object.
      #
      def output_file=(new_output_file)
        @output_file = File.new(new_output_file.respond_to?(:path) ? new_output_file.path : new_output_file.to_s, 'w')
      end
      
      # Runs the BLAT instance against the query provided returning a file to the results of the alignment. 
      #
      # * *Args*    :
      #   - +query+ -> A String or File.  The former can be an explicit query sequence to be used or the path to a file.
      # * *Returns* :
      #   - A Tempfile object containing the result of the BLAT alignment.
      #
      def run(query, options = {})
        options = { verbose: false }.merge(options)
        
        # Format the query for execution
        query_path = case
        when File.exists?(query) then query
        else
          # Write the query to a file.
          tempfile = Tempfile.new("blast_query")
          tempfile.puts ">Query"
          tempfile.puts query
          tempfile.close
          
          tempfile.path
        end
        
        # Execute the command
        CommandLine.run(prepare_command(query_path)) do |f|
          puts f.read
        end
        
        # Grab the output file
        output_file.rewind
        output_file
      end
      
      class << self
        
        # Returns true if BLAT is installed on the system and is accessible using the current path.
        #
        # * *Returns* :
        #   - A boolean
        #
        def installed?
          begin
            CommandLine.run('blat') { |f| f.read }
            true
          rescue
            false
          end
        end
      end
      
      private
      
      # Returns a string suitable for use in running the alignment on the commad line.
      #
      # * *Args*    :
      #   - +query_path+ -> A String with the path to the query file.
      # * *Returns* :
      #   -
      #
      def prepare_command(query_path)
        options = []
        OPTION_MAPPING.each do |option, value_lambda|
          options << value_lambda.call(self.send(option)) if self.send(option)
        end
        
        command = "blat #{database} #{query_path} #{output_file.path} #{options.compact.join(" ")}"
      end
    end
  end
end
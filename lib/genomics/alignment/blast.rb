require 'tempfile'

module Genomics
  module Alignment
    # This object handles creating BLAST queries and sending them to the command line to be run.
    class BLAST
      # The possible choices of scoring matricies supported by the standard BLAST distribution
      SCORING_MATRICES = %w{BLOSUM62 BLOSUM80 BLOSUM50 BLOSUM45 BLOSUM90 PAM30 PAM70 PAM250}
      
      # A hash matching the optional attributes to the corresponding flags for command line execution
      OPTION_MAPPING = { database:        ->(value) { "-db #{value}"},
                         output_file:     ->(value) { "-out #{value.path}"},
                         e_value:         ->(value) { "-evalue #{value.to_f}" }, 
                         out_format:      ->(value) {
                           numberCode = case value
                           when :xml then 5
                           when :tab then 7
                           when :csv then 10
                           else value
                           end
                           
                           "-outfmt #{numberCode}"
                         },
                         scoring_matrix:  ->(value) { "-matrix #{value}" },
                         allow_gaps:      ->(value) { '-ungapped -comp_based_stats F' if value },
                         filter_query:    ->(value) { '-dust no -seg no' unless value }
                        }
      
      attr_reader :blast_command
      attr_accessor :database, :e_value, :out_format, :scoring_matrix, :word_size, :allow_gaps, :filter_query
      
      # * *Args*    :
      #   - +blast_command+ -> A string with the path to the BLAST program to be run or a symbol represeting a standard BLAST
      #   algorithm assumed to be in the PATH.
      # * *Returns* :
      #   -
      # * *Raises* :
      #   - ++ ->
      #
      def initialize(blast_command, options = {})
        @blast_command = blast_command
        
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
        @output_file ||= Tempfile.new("#{@blast_command}_out")
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
        @output_file = File.new(new_output_file.respond_to?(:path) ? new_output_file.path : new_output_file.to_s)
      end
      
      # Runs the BLAST instance against the query provided returning a file to the results of the alignment. 
      #
      # * *Args*    :
      #   - +query+ -> A String or File.  The former can be an explicit query sequence to be used or the path to a file.
      # * *Returns* :
      #   - A Tempfile object containing the result of the BLAST alignment.
      # * *Raises* :
      #   - ++ ->
      #
      def run(query, options = {})
        options = { verbose: false }.merge(options)
        
        # Format the query for execution
        query_path = case
        when File.exists?(query) then query
        else
          # Write the query to a file.
          tempfile = Tempfile.new("#{@blast_command}_query")
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
          unless self.send(option).nil?
            options << value_lambda.call(self.send(option))
          end
        end
        
        command = "#{@blast_command} -query #{query_path} #{options.compact.join(" ")}"
      end
    end
  end
end
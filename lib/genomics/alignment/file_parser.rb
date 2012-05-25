require 'pathname'

module Genomics
  module Alignment
    # This module implements a method for parsing through files and handing character conversions errors.
    class FileParser < Genomics::FileParser
      
      class << self
        # Parses the supplied file into a list of alignment objects.
        #
        # * *Args*    :
        #   - +file+ -> A string representing the name of the file to parse. 
        # * *Options*    :
        #   - +aggregate_hits+ -> If true, rather than returning a list of hits, the hits are aggregated into a multilayered hash
        #                         of hits based on they query and subject.
        def parse_file(file, options ={})
          options = { aggregate_hits: false }.merge(options)
          
          # Create the progress bar
          `grep -v '#' #{file} | wc -l` =~ /\d+/
          pbar = ProgressBar.new("Parsing #{Pathname.new(file).basename}", $~[0].to_i, STDOUT)
          
          # Defer to the parent class to implement the details of parsing individual lines.
          hits = []
          super(file, "\t") do |row|
            pbar.inc
            
            attributes = { query:                row[0],
                           subject:              row[1],
                           query_start:          row[6].to_i,
                           query_end:            row[7].to_i,
                           subject_start:        row[8].to_i,
                           subject_end:          row[9].to_i,
                           e_value:              EValue.new(row[10]), 
                           bit_score:            row[11].to_f,
                           percentage_identity:  row[2].to_f,
                           alignment_length:     row[3].to_i,
                           mismatches:           row[4].to_i,
                           gap_openings:         row[5].to_i }
            
            hits << Hit.new(attributes)
          end
          
          return hits unless options[:aggregate_hits]

          # Group the hits based on the specifics of what was matched.
          aggregated_hits = {}
          hits.each do |hit|
            aggregated_hits[hit.query] ||= {}
            aggregated_hits[hit.query][hit.subject] ||= []
            aggregated_hits[hit.query][hit.subject] << hit
          end
          
          aggregated_hits
        end
      end
    end
  end
end
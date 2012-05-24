module Genomics
  module Alignment
    # This object represents a fundamental unit of an alignment.
    class Hit
      ATTRIBUTE_NAMES = [ :query, :subject, :query_start, :query_end, :subject_start, :subject_end, :e_value, :bit_score,
                          :percentage_identity, :alignment_length, :mismatches, :gap_openings]
      
      attr_accessor *ATTRIBUTE_NAMES
      
      alias :query_stop :query_end       
      alias :query_stop= :query_end=       
      
      def initialize(attributes = {})
        attributes.each do |name, value|
          send("#{name}=", value)
        end
      end
      
      def inspect
        attribute_strings = ATTRIBUTE_NAMES.map { |attribute| "#{attribute}: \"#{send(attribute)}\"" }
        "#<Genomics::Alignment::Hit #{attribute_strings.join(", ")} >" 
      end
      
      # Returns the length of the alignment as matched on the query sequence.  Strandedness issues are taken into account
      # when calculating this value.
      #
      # * *Returns* :
      #   - An integer describing the length of the alignment on the query sequence.
      #
      def query_length
        if query_end > query_start
          query_end - query_start + 1
        else
          query_start - query_end + 1
        end
      end
    end
  end
end
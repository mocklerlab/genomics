module Genomics
  module IO
    module BLAST
      # This object represents a fundamental unit of an alignment.
      class Hit
        ATTRIBUTE_NAMES = [ :query, :subject, :query_start, :query_end, :subject_start, :subject_end, :e_value, :bit_score,
                            :alignment_length, :gap_openings, :query_frame, :subject_frame,
                            :identities, :positives, :query_sequence, :subject_sequence, :midline]
      
        attr_accessor *ATTRIBUTE_NAMES
      
        alias :query_stop :query_end       
        alias :query_stop= :query_end=       
      
        def initialize(attributes = {})
          attributes.each do |name, value|
            send("#{name}=", value)
          end
        end
        
        def <=>(other)
          bit_score <=> other.bit_score
        end
        
        # Allows the e_value to be set for the hit, converting from both string and numerics to an EValue object.
        # 
        # * *Argument* :
        #   - +e_value+ -> A string or numerical value used to instantiate the e_value.
        # * *Returns* :
        #   - An integer describing the length of the alignment on the query sequence.
        #
        def e_value=(new_e_value)
          @e_value = new_e_value.is_a?(EValue) ? new_e_value : EValue.new(new_e_value)
        end
      
        def inspect
          attribute_strings = ATTRIBUTE_NAMES.map { |attribute| "#{attribute}: \"#{send(attribute)}\"" }
          "#<Genomics::Alignment::Hit #{attribute_strings.join(", ")}>" 
        end
        
        # Returns a boolean indicating the orientation of the hit relative to the strandedness of the subject or query.
        # The hit is on the forward strand if the positions are increasing.
        #
        # * *Options*    :
        #   - +:on+ -> A symbol indicating which of the query or the subject to use.  (Default: :subject)
        # * *Returns* :
        #   - A boolean indicating if the Hit is on the forward strand.
        #
        def forward_strand?(options = {})
          options = { on: :subject }.merge(options)
          
          case options[:on]
          when :subject
            subject_end > subject_start
          when :query
            query_end > query_start
          end 
        end
        
        # Returns the length of the alignment as matched on the specified sequence.  Strandedness issues are taken into account
        # when calculating this value.
        #
        # * *Options*    :
        #   - +:on+ -> A symbol indicating which of the query or the subject to use.  (Default: :subject)
        # * *Returns* :
        #   - An integer describing the length of the alignment on the query sequence.
        #
        def length(options = {})
          options = { on: :subject }.merge(options)
          
          value = case options[:on]
          when :subject
            forward_strand?(on: options[:on]) ? subject_end - subject_start : subject_start - subject_end
          when :query
            forward_strand?(on: :query) ? query_end - query_start : query_start - query_end
          end
          
          value + 1
        end
        
        # Returns the number of mismatches.
        #
        # * *Returns* :
        #   - Integer
        #
        def mismatches
          alignment_length - identities - gap_openigns
        end
        
        # Returns the hit object with the query and subject related fields interchanged.
        #
        # * *Returns* :
        #   - The Hit object.
        #
        def transpose!
          @query, @subject = subject, query
          @query_start, @query_end, @subject_start, @subject_end = subject_start, subject_end, query_start, query_end
          
          self
        end
      end
    end
  end
end
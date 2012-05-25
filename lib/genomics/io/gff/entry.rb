require 'stringio'

module Genomics
  module IO
    module GFF

      # Genomics::IO::GFFEntry represents a single entry in a GFF file.  This entry can consist of more than one
      # line, but all of the line are grouped together via so structure as reflected in the ID(s).
      class Entry
        VALID_STRANDS = %w{+ - .}
    
        attr_accessor :seqid, :source, :type, :score, :strand, :phase, :attributes
        attr_reader :regions
  
        # Creates a new object from the supplied attributes.
        #
        # * *Attributes*    :
        #   - +seqid+ -> The id of the landmark used to establish the coordinate system for the entry.
        #   - +source+ -> The algorithm or operating procedure used to generate the entry.
        #   - +type+ -> The term from the sequence ontology that defines the nature of the record.
        # * *Raises* :
        #   - +ArgumentError+ -> If any of the *seqid*, *source*, or *type* attributes are missing.
        #
        def initialize(__attributes__ = {})
          # Ensure that the required attribues are present.
          [:seqid, :source, :type].each do |required_attribute|
            raise ArgumentError, "Missing attribute #{required_attribute}" unless __attributes__[required_attribute] 
          end

          # Set default attribtues
          __attributes__ = { score: '.', strand: '.', phase: '.' }.merge(__attributes__)

          # Pull out the start and end attributes if present
          start = __attributes__.delete(:start)
          stop = __attributes__.delete(:end) || __attributes__.delete(:stop)
      
          # Create the regions object that holds the specifics of the sequence positioning.
          @regions = EntryRegions.new(self)

          @regions.create(start: start, end: stop) if start && stop
    
          __attributes__.each do |name, value|
            send("#{name}=", value)
          end
        end
    
        def <=>(other)
          seqid_sort = @seqid <=> other.seqid
          return seqid_sort if seqid_sort != 0

          start <=> other.start
        end
    
        # Returns the start position.
        #
        # * *Returns* :
        #   - An integer representing the start of the entry in the forward direction on the landmark sequence.
        def start
          regions.min.start
        end
  
        # Returns the stop position. 
        #
        # * *Returns* :
        #   - An integer representing the end of the entry in the forward direction on the landmark sequence.
        def end
          regions.min.end
        end
    
        # Returns the score associated with the entry.  If the entry contains regions with different scores, this returns
        # the smallest value.
        #
        def score
          @score || regions.map(&:score).min
        end
    
        # Sets the orientation of the entry relative to the landmark sequence.
        #
        # * *Args*    :
        #   - +strand+ -> A string definining the "strand" value of the entry. Possible values are: '+', '-', '.'
        # * *Returns* :
        #   -
        # * *Raises* :
        #   - +ArgumentError+ -> If the new strand is not of a valid type.
        #
        def strand=(new_strand)
          raise ArgumentError, "The strand provided is not valid." unless VALID_STRANDS.include?(new_strand)
          @strand = new_strand
        end
    
        # Returns true if the entry is oriented on the same strand as the landmark sequence.
        #
        def forward_strand?
          @strand == '+'
        end
    
        # Returns true if the entry is oriented on the reverse complement of the landmark sequence.
        #
        def reverse_strand?
          @strand == '-'
        end
    
        # Returns the entry as a string suitable for use in a GFF3 file.
        #
        # * *Returns* :
        #   - A string representing the entry in GFF3 format.
        #
        def to_gff
          gff_string = StringIO.new
      
          case type
          when :match, :EST_match
            # The entry is due to an alignment.  All of the regions will be printed as part of one discontiguous entry.
            common_values = "#{seqid}\t#{source}\t#{type}"
            @regions.sort.each do |region|
              segment_attributes = sort_attributes(region.attributes.merge(attributes))
              segment_attributes.map! { |attribute, value| "#{attribute}=#{value}"}
              gff_string.puts "#{common_values}\t#{region.start}\t#{region.stop}\t#{region.score}\t#{strand}\t#{phase}\t#{segment_attributes.join(";")}"
            end
          end
      
          gff_string.string
        end
    
        private
    
        # Sorts the attributes hash into the correct order for printing.
        # TODO: Slow and inefficient
        def sort_attributes(attributes_hash)
          attributes_order = %w{ID Name * Target Gap}
      
          attributes_hash.sort do |(a_key, a_value), (b_key, b_value)|
            # Reduce each of the attribute pairs to their correct index in the order
            a_key_val = attributes_order.include?(a_key) ? a_key : '*'
            b_key_val = attributes_order.include?(b_key) ? b_key : '*'
        
            index_sort = attributes_order.index(a_key_val) <=> attributes_order.index(b_key_val)
        
            index_sort == 0 ? a_key_val <=> b_key_val : index_sort
          end
        end
  
      end

      class EntryRegions
        include Enumerable
    
        def initialize(entry)
          @entry = entry
          @regions = []
        end
    
        def each
          @regions.each { |region| yield region }
        end
    
        # Creates a Genomics:IO:GFF:EntryRegion object using the supplied attributes and adds it to the collection.
        # This ensures that the orientation of the start and end positions is consistent with the strand of the entry, if defined.
        #
        # * *Attributes*    :
        #   - +start+ -> An integer specifying the position where the region begins.
        #   - +end+ -> An integer specifying the position where the region stops.
        # * *Returns* :
        #   - The newly created Genomics::IO:GFF:EntryRegion object.
        #
        def create(attributes = {})
          # Ensure that positions are always increasing.
          positions = [attributes.delete(:start), attributes.delete(:end)]
          positions.sort! if @entry.strand != '.'
          start, stop = positions
      
          # Add the region
          @regions << EntryRegion.new(attributes.merge(start: start, end: stop))
      
          @regions.last
        end
  
      end
  
      # Represents a contiguous part of the entry, which is taken to be on the same strand as the wrapping entry.
      #
      class EntryRegion
        attr_accessor :start, :end, :score, :phase, :attributes
    
        alias :stop :end
        alias :stop= :end=
    
        # * *Attributes*    :
        #   - +start+ -> The position on the landmark sequence where the region starts.
        #   - +end+ -> The positio on the landmark sequence where the region ends.
        #   - +score+ -> A numeric value representing the score of the sequence, which is typically uses for regions generated via alignments.
        #   - +phase+ -> 
        def initialize(__attributes__ = {})
          __attributes__.each do |name, value|
            send("#{name}=", value)
          end
        end
    
        def <=>(other)
          start <=> other.start
        end
      end
    end
  end
end
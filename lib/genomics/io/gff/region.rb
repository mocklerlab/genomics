module Genomics
  module IO
    module GFF
      # Genomics::IO::GFF::Regions extends enumerable to provide custom accessor functions for manipulating the regions of a
      # GFF::Feature.
      class Regions
        include Enumerable
    
        # TODO: Use method_undefined instead defering to the unerlying array.
        alias :length :count
        alias :size :count
    
        def initialize(feature)
          @feature = feature
          @regions = []
        end
      
        # Iterates through each of the Regions increaing in position.
        #
        # * *Options*   :
        #   - +strand_order+ -> A boolean that indicates whether the order should account for the strand the feature is on (Default: false). 
        #
        def each(options = {})
          options = { strand_order: false }.merge(options)
          
          # Sort the regions
          sorted_regions = @regions.sort
          sorted_regions.reverse! if options[:strand_order] && @feature.reverse_strand?

          sorted_regions.each { |region| yield region }
        end
        
        # Returns the last object under the default ordering.
        #
        # * *Returns* :
        #   - A Region
        #
        def last
          @regions.sort.last
        end
    
        # Creates a Genomics:IO:GFF:Region object using the supplied attributes and adds it to the collection.
        # This ensures that the orientation of the start and end positions is consistent with the strand of the feature, if defined.
        #
        # * *Attributes*    :
        #   - +start+ -> An integer specifying the position where the region begins.
        #   - +end+ -> An integer specifying the position where the region stops.
        # * *Returns* :
        #   - The newly created Genomics::IO:GFF:EntryRegion object.
        #
        def create(attributes = {})
          # Add the region
          @regions << Region.new(attributes)
      
          @regions.last
        end
        
      end
  
      # Represents a contiguous part of the Feature, which is taken to be on the same strand as the wrapping feature.
      #
      class Region
        attr_accessor :start, :end, :score, :phase, :attributes
    
        alias :stop :end
        alias :stop= :end=
    
        # * *Attributes*    :
        #   - +start+ -> The position on the landmark sequence where the region starts.
        #   - +end+ -> The position on the landmark sequence where the region ends.
        #   - +score+ -> A numeric value representing the score of the sequence, which is typically uses for regions generated via alignments.
        #   - +phase+ -> A intger indicate the number of bases that need to be removed from the beginning of this region to reach the next codon.
        def initialize(__attributes__ = {})
          __attributes__.each do |name, value|
            send("#{name}=", value)
          end
        end
        
        # Sets the attributes.
        #
        # * *Args*    :
        #   - +new_attributes+ -> A Hash or string in the valid GFF3 format.
        # * *Returns* :
        #   - The new attributes Hash
        def attributes=(new_attributes)
          @attributes = if new_attributes.is_a?(Hash)
            cloned_attributes = new_attributes.clone
            cloned_attributes.delete(:ID)
            cloned_attributes.delete(:Name)
            cloned_attributes
          else
            Feature.parse_attributes(new_attributes, except: [:ID, :Name])
          end
        end

        # Sets the end position of the region converting the argument to an integer.
        #
        # * *Returns* :
        #   - An integer
        #
        def end=(new_end)
          @end = new_end ? new_end.to_i : nil
        end
        
        # Returns the length of the region.
        #
        # * *Returns* :
        #   - An integer
        #
        def length
          @end - @start + 1
        end
        
        # Sets the phase of the region converting the argument to an integer.
        #
        # * *Returns* :
        #   - An integer
        #
        def phase=(new_phase)
          @phase = new_phase && new_phase != '.' ? new_phase.to_i : nil
        end
        
        # Sets the score of the region converting the argument to a float.
        #
        # * *Returns* :
        #   - An float
        #
        def score=(new_score)
          @score = new_score && new_score != '.' ? new_score.to_f : nil
        end
        
        # Sets the start position of the region converting the argument to an integer.
        #
        # * *Returns* :
        #   - An integer
        #
        def start=(new_start)
          @start = new_start ? new_start.to_i : nil
        end
        
        def <=>(other)
          start <=> other.start
        end
      end
    end
  end
end
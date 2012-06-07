require 'stringio'

module Genomics
  module IO
    module GFF
      # Genomics::IO::GFF::Featrue represents a single feature in a GFF file.  This entry can consist of more than one
      # line, but all of the line are grouped together via so structure as reflected in the ID(s).
      class Feature
        VALID_STRANDS = %w{+ - . ?}
    
        attr_accessor :seqid, :source, :type, :strand, :attributes
        attr_reader :regions
  
        # Creates a new object from the supplied attributes.
        #
        # * *Attributes*    :
        #   - +seqid+ -> The id of the landmark used to establish the coordinate system for the entry.
        #   - +source+ -> The algorithm or operating procedure used to generate the entry.
        #   - +type+ -> The term from the sequence ontology that defines the nature of the record.
        #   - +attributes+ ->
        # * *Raises* :
        #   - +ArgumentError+ -> If any of the *seqid*, *source*, or *type* attributes are missing.
        #
        def initialize(__attributes__ = {})
          # Ensure that the required attribues are present.
          # [:seqid, :source, :type].each do |required_attribute|
          #   raise ArgumentError, "Missing attribute #{required_attribute}" unless __attributes__[required_attribute] 
          # end

          # Set default attribtues
          __attributes__ = { strand: '.' }.merge(__attributes__)

          # Pull out the start and end attributes if present
          start = __attributes__.delete(:start)
          stop = __attributes__.delete(:end) || __attributes__.delete(:stop)
      
          # Create the regions object that holds the specifics of the sequence positioning.
          @regions = Regions.new(self)

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
    
        # Sets the attributes.
        #
        # * *Args*    :
        #   - +new_attributes+ -> A Hash or string in the valid GFF3 format.
        # * *Returns* :
        #   - The new attributes Hash
        def attributes=(new_attributes)
          @attributes = if new_attributes.is_a?(Hash)
            new_attribtues
          else
            Feature.parse_attributes(new_attributes, only: [:ID, :Name])
          end
        end
    
        # Returns the ID of the feature.
        #
        # * *Returns* :
        #   - A string
        #
        def id
          @attributes[:ID]
        end
        
        # Returns the Name of the feature.
        #
        # * *Returns* :
        #   - A string
        #
        def name
          @attributes[:Name]
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
        alias :stop :end
    
        # Returns the score associated with the entry.  If the entry contains regions with different scores, this returns
        # the largest value.
        #
        def score
          regions.map(&:score).max
        end
    
        # Sets the orientation of the entry relative to the landmark sequence.
        #
        # * *Args*    :
        #   - +strand+ -> A string definining the "strand" value of the entry. Possible values are: '+', '-', '.'
        # * *Returns* :
        #   - A string
        # * *Raises* :
        #   - +ArgumentError+ -> If the new strand is not of a valid type.
        #
        def strand=(new_strand)
          raise ArgumentError, "The strand provided is not valid." unless VALID_STRANDS.include?(new_strand)
          @strand = new_strand
        end
    
        # Sets the type ensuring that it is a valid symbol.
        #
        # * *Args*    :
        #   - +new_type+ -> A String os symbol specifying the type attribute
        # * *Returns* :
        #   - A symbol representing the new attribute value.
        #
        def type=(new_type)
          @type = new_type.to_sym
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
    
        class << self
          
          # Takes a GFF attribute string and converts it into a hash, properly accounting for special characters and encodings.
          #
          # * *Args*    :
          #   - +attributes_string+ -> A String represinting a GFF3 attributes field.
          # * *Options* : 
          #   - +only+  : A symbol or array of symbols explicitly specifying which attributes should be retained. 
          #   - +except+  : A symbol or array of symbols explicitly specifying which attributes should not be included. 
          # * *Returns* : 
          #   - A Hash
          #
          def parse_attributes(attributes_string, options = {})
            options = { only: [], except: [] }.merge(options)
            options[:only] = [options[:only]] unless options[:only].is_a?(Array)
            options[:except] = [options[:except]] unless options[:except].is_a?(Array)
            
            attributes_hash = {}
            
            # Split based on ; and then key value pairs based on =
            attribute_pairs = attributes_string.split(";")
            attribute_pairs.each do |pair|
              attribute, value = pair.split("=")
              attribute = attribute.to_sym
              
              # Determine if this attribute should be parsed
              next if options[:only].any? && !options[:only].include?(attribute)
              next if options[:except].any? && options[:except].include?(attribute)
              
              # Convert the value to an array based on the semantic character ,
              values = value.split(',')
              
              # Undo special character encodings
              values.each do |value|
                value.gsub!(/%\d[A-Z0-9]/) do |match|
                  { '%3D' => '=', '%3B' => ';', '%2C' => ',', '%09' => "\t"}[match.to_s]
                end
              end
              
              # Reduce single valued arrays
              value = values.size == 1 ? values.first : values
              
              attributes_hash[attribute] = value
            end
            
            attributes_hash
          end
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
    end
  end
end

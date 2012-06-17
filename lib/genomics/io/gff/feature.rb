require 'stringio'

module Genomics
  module IO
    module GFF
      # Genomics::IO::GFF::Features extends enumerable to provide custom accessor functions for manipulating the child features of a
      # Feature.
      class Features
        include Enumerable
    
        def initialize(feature)
          @feature = feature
          @features = []
        end
        
        # Iterates through each of the Features increaing in position.
        #
        # * *Options*   :
        #   - +strand_order+ -> A boolean that indicates whether the order should account for the strand the feature is on (Default: false). 
        #
        def each(options = {})
          options = { strand_order: false }.merge(options)
          
          # Sort the regions
          sorted_features = @features.sort
          sorted_features.reverse! if options[:strand_order] && @feature.reverse_strand?
          
          sorted_features.each { |feature| yield feature }
        end
        
        # Takes the given query and finds the feature that has a matching ID or returns nil.
        #
        # * *Args*    :
        #   - +id+ -> The id of the feature to search for
        # * *Returns* :
        #   - A Feature
        # TODO: Use method missing to add arbitrary finders.
        def find_by_id(id)
          @features.find { |feature| feature.id == id }
        end
        
        # Returns the last object under the default ordering.
        #
        # * *Returns* :
        #   - A Feature
        #
        def last
          @features.sort.last
        end
    
        # Creates a Genomics:IO:GFF:FEature object using the supplied attributes and adds it to the collection.
        #
        # * *Attributes*    :
        #   - +start+ -> An integer specifying the position where the region begins.
        #   - +end+ -> An integer specifying the position where the region stops.
        # * *Returns* :
        #   - Genomics::IO:GFF:Feature
        #
        def create(attributes = {})
          # Automatically set some attributes if possible
          attributes = { seqid: @feature.seqid, source: @feature.source, strand: @feature.strand }.merge(attributes)
          
          @features << Feature.new(attributes)
          @features.last
        end
  
        def respond_to?(method_name, include_private = false)
          @features.respond_to?(method_sym) ? true : super
        end

        def method_missing(method_sym, *arguments, &block)
          if @features.respond_to?(method_sym)
            @features.send(method_sym, *arguments, &block)
          else
            super
          end
        end
        
      end
      
      # Genomics::IO::GFF::Feature represents a single feature in a GFF file.  This entry can consist of more than one
      # line, but all of the line are grouped together via so structure as reflected in the ID(s).
      class Feature
        VALID_STRANDS = %w{+ - . ?}
    
        attr_accessor :seqid, :source, :type, :strand, :attributes
        attr_reader :features, :regions, :derivatives
  
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
          # TODO: This was causing problems when features got serialized to a db, since null features were being instantiated
          # Ensure that the required attribues are present.
          # [:seqid, :source, :type].each do |required_attribute|
          #   raise ArgumentError, "Missing attribute #{required_attribute}" unless __attributes__[required_attribute] 
          # end
          
          @attributes = {}
          
          # Set default attribtues
          __attributes__ = { strand: '.' }.merge(__attributes__)

          # Pull out the start and end attributes if present
          start = __attributes__.delete(:start)
          stop = __attributes__.delete(:end) || __attributes__.delete(:stop)
      
          # Create the regions object that holds the specifics of the sequence positioning.
          @regions = Regions.new(self)
          @features = Features.new(self)
          @derivatives = Features.new(self)
          
          @regions.create(start: start, end: stop) if start && stop
    
          __attributes__.each do |name, value|
            send("#{name}=", value)
          end
        end
    
        def <=>(other)
          # Sort by the landmark sequence, this could possibly be numerical or a mix
          seqid_sort = if seqid =~ /\d+$/
            common_regex = /#{other.seqid.gsub(/\d+$/, '(\d+)')}/
            other_number = $~[0] if $~
            if other_number && seqid.match(common_regex)
              $~[1].to_i <=> other_number.to_i
            else
              @seqid <=> other.seqid
            end
          else
            @seqid <=> other.seqid
          end
          return seqid_sort if seqid_sort != 0
          
          # Sort by position on the landmark sequence
          start_sort = start <=> other.start
          return start_sort if start_sort != 0
          
          # Otherwise sort by the feature ordered in a way to provide meaning
          Feature.type_comparison(self, other)
        end
      
        # Returns a string describing the state of the object.
        #
        # * *Returns* :
        #   - The new attributes Hash
        #
        def inspect
          attributes = ["seqid=\"#{@seqid}\"", "source=\"#{@source}\"", "type=\"#{@type}\"", "strand=\"#{@strand}\""]
          attributes += @attributes.map { |(attribute, value)| "#{attribute}=\"#{value}\"" }
          
          "#<#{self.class} #{attributes.join(', ')}>"
        end
    
        # Sets the attributes.
        #
        # * *Args*    :
        #   - +new_attributes+ -> A Hash or string in the valid GFF3 format.
        # * *Returns* :
        #   - The new attributes Hash
        # TODO: Remove the parent attribute and construct an association out of it.
        def attributes=(new_attributes)
          acceptable_attributes = [:ID, :Name, :Note, :Alias, :Parent]
          @attributes = if new_attributes.is_a?(Hash)
            # Symbolify the keys in the hash
            symbolized_attributes = Hash[new_attributes.map { |(attribute, value)| [attribute.to_sym, value] }]
            
            # Pull out the acceptable attributes
            detected_attributes = {}
            acceptable_attributes.each do |attribute|
              detected_attributes[attribute] = symbolized_attributes[attribute] if symbolized_attributes[attribute]
            end
            
            detected_attributes
          else
            Feature.parse_attributes(new_attributes, only: acceptable_attributes)
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
        
        # Allows the ID attribute to be directly set.  This also cascades, updating the ids of child features.
        #
        # * *Returns* :
        #   - A string
        #
        def id=(new_id)
          @features.each do |feature| 
            feature.id = feature.id.gsub(/#{@attributes[:ID]}/, new_id) if feature.id =~ /#{@attributes[:ID]}/
            
            # Update child parent attributes as well
            if feature.attributes[:Parent].is_a?(Array)
              feature.attributes[:Parent].each do |parent|
                parent.gsub!(/#{@attributes[:ID]}/, new_id) if parent =~ /#{@attributes[:ID]}/
              end
            elsif (parent = feature.attributes[:Parent]) =~ /#{@attributes[:ID]}/
              feature.attributes[:Parent] = parent.gsub(/#{@attributes[:ID]}/, new_id) 
            end
          end
          
          @attributes[:ID] = new_id
        end
        
        # Returns the Name of the feature.
        #
        # * *Returns* :
        #   - A string
        #
        def name
          @attributes[:Name]
        end
        
        # Allows the Name attribute to be directly set.  This also cascades, updating the names of child features.
        #
        # * *Returns* :
        #   - A string
        #
        def name=(new_name)
          @features.each do |feature| 
            feature.name = feature.name.gsub(/@attribtues[:Name]/, new_name) if feature.name =~ /#{@attributes[:Name]}/
          end
          
          @attributes[:Name] = new_name
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
        # * *Options* :
        #   - +strand_order+ -> A boolean indicating whether the features/regions should ordered based on the natural direction of the
        #   strand that they lie on (Default: true).
        #   - +sequential_order+ -> A boolean signifying that multipart regions should be broken up and printed based on their positions
        #   and not as a single contiguous block (Default: true).
        # * *Returns* :
        #   - A string representing the entry in GFF3 format.
        #
        def to_gff(options = {})
          options = { strand_order: true, sequential_order: true }.merge(options)
          gff_string = StringIO.new
      
          # Go through the regions in order printing them
          @regions.each(strand_order: options[:strand_order]) do |region|
            gff_string.puts to_gff_for_region(region)
          end
          
          gff_string.string
        end
        
        # Takes a region and only prints the region and the subfeatures that fall in it.  Takes the same options as to_gff.
        #
        # * *Args*
        #   - +region+ -> The region to print gff subfeatures for
        # * *Returns* :
        #   - An array of hashes with feature and region keys
        #
        def to_gff_for_region(region, options = {})
          options = { strand_order: true, sequential_order: true }.merge(options)
          gff_string = StringIO.new
          
          # Print the region for the feature
          common_values = "#{seqid}\t#{source}\t#{type}"
          region_attributes = Feature.encode_attributes(region.attributes.merge(attributes))
          gff_string.puts "#{common_values}\t#{region.start}\t#{region.stop}\t#{region.score || '.'}\t#{strand}\t#{region.phase || '.'}\t#{region_attributes}"
          
          # Find all of the subfeatures in the region and print those features
          contained_features = @features.select { |feature| region.include?(feature) }
          
          # TODO: Implement the false variant for sequential sorting
          
          # Get all of the regions associated with all of the child features
          child_feature_regions = contained_features.map do |feature|
            feature.regions.map do |region|
              { feature: feature, region: region }
            end
          end.flatten

          # Sort the regions
          child_feature_regions.sort! do |a, b|
            region_sort = options[:strand_order] && reverse_strand? ? b[:region].stop <=> a[:region].stop : a[:region] <=> b[:region]
            region_sort != 0 ? region_sort : Feature.type_comparison(a[:feature], b[:feature])
          end

          # Print each of the regions
          child_feature_regions.each do |region_feature_hash|
            gff_string.puts region_feature_hash[:feature].to_gff_for_region(region_feature_hash[:region]) 
          end
          
          gff_string.string
        end
          
        class << self
          
          # Takes a hash of attribtues and encodes them properly for injection into a GFF3 file.
          #
          # * *Args*    :
          #   - +attributes+ -> A hash of attributes to be written as a string
          # * *Returns* :
          #   - A string
          #
          def encode_attributes(attributes)
            # Encode the GFF special characters
            encoded_attributes = {}
            attributes.each do |attribute, value|
              # The value could possibly be an array, so convert everything to an array
              values = !value.is_a?(Array) ? [value] : value
              values.map! do |value|
                value.to_s.gsub(/[=;,\t]/) do |match|
                  { '=' => '%3D', ';' => '%3B', ',' => '%2C', "\t" => '%09'}[match.to_s]
                end
              end.join(',')

              encoded_attributes[attribute] = values.join(',')
            end

            # Join key/value pairs and then concatenate
            sort_attributes(encoded_attributes).map do |(name, value)|
              "#{name}=#{value}"
            end.join(';')
          end
          
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
          
          # Returns an integer reflecting an ordering comparison between two features based on type.
          #
          # * *Args*    :
          #   - +a+ -> The first feature to compare.
          #   - +b+ -> The second feature to compare.
          # * *Returns* : 
          #   - -1, 0, or 1
          #
          def type_comparison(a, b)
            # Set the precendence for features
            type_order = [:exon,
                          :intron,
                          :CDS,
                          :five_prime_UTR,
                          :three_prime_UTR,
                          :transcription_start_site,
                          :transcription_end_site,
                          :start_codon,
                          :stop_codon ]
            (type_order.index(a.type) || type_order.length) <=> (type_order.index(b.type) || type_order.length)
          end
          
          private
          
          # Sorts the attributes hash into the correct order for printing.
          # TODO: Slow and inefficient
          def sort_attributes(attributes_hash)
            attributes_order = %w{ID Name * Target Gap}

            attributes_hash.sort do |(a_key, a_value), (b_key, b_value)|
              a_key = a_key.to_s
              b_key = b_key.to_s
              
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
end

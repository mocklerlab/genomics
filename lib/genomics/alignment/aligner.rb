module Genomics
  module Alignment
    class Aligner
      class << self
        
        private
        
        # Create the list of hits to be converted into single entries based on query, subject, strand, and separation i.e hits widely 
        # will be treated as distinct matches to the query.
        #
        # * *Args*    :
        #   - +aggregated_hits+ -> A multi-dimmensional hash with query and then subject keys mapping to hits.
        # * *Returns* :
        #   - An Array of arrays of which the later correspond to the clusters.
        # TODO: Pending deletion
        def cluster_aggregated_hits(aggregated_hits, options = {})
          # Create the progress bar
          total_hit_count = aggregated_hits.values.inject(0) { |sum, subject_hash| sum + subject_hash.values.length }
          pbar = ProgressBar.new("Clustering", total_hit_count, STDOUT)

          clustered_hits = []
          aggregated_hits.each do |query, subjects|
            subjects.each do |subject, hits|
              pbar.inc
              
              # Separate the hits by strand
              positive_hits, negative_hits = [], []
              hits.sort_by(&:query_start).each { |hit| hit.query_end > hit.query_start ? positive_hits << hit : negative_hits << hit }
              
              # Call the hit clusters based on separation on the query sequence.  Use the average length of the hit to set the scale.
              average_length = hits.inject(0) { |sum, hit| sum + hit.query_length } / hits.size.to_f
              clustered_hits += cluster(positive_hits, average_length * 10) if positive_hits.any?
              clustered_hits += cluster(negative_hits, average_length * 10) if negative_hits.any?
            end
          end
          
          clustered_hits
        end
        
        # Returns a set of clusters of hits that can be taken as a single discontiguous entry.
        #
        # * *Args*    :
        #   - +hits+ -> The Hit objects to be clustered together, which are assumed to be sorted in increasing position.
        #   - +cutoff+ -> A numeric value setting the maximum distance that two hits can be separate and still be called in the same cluster.
        # * *Returns* :
        #   - An array of arrays of clustered hits determined by the cutoff.
        #
        def cluster(hits, options = {})
          options = { cluster_on: :query, cutoff: 10000 }.merge(options)
          
          # Set the start stop methods
          start_method, stop_method = options[:cluster_on] == :query ? [:query_start, :query_end] : [:subject_start, :subject_end]
          
          clusters = []
          running_hits = []
          hits.each do |hit|
            if running_hits.empty? || hit.send(start_method) - running_hits.last.send(stop_method) < options[:cutoff]
              running_hits << hit
            else
              # Add the running hits as a cluster and reset
              clusters << running_hits
              running_hits = [hit]
            end
          end
          
          # Add any remaining running hits as a cluster
          clusters << running_hits if running_hits.any?
          
          clusters
        end
        
        def cluster_hits(hits, options = {})
          options = { cluster_on: :query }.merge(options)
          
          # Separate the hits based on their orientations
          positive_hits, negative_hits = [], []
          case options[:cluster_on]
          when :query
            hits.sort_by(&:query_start).each { |hit| hit.query_end > hit.query_start ? positive_hits << hit : negative_hits << hit }
          when :subject
            hits.sort_by(&:subject_start).each { |hit| hit.subject_end > hit.subject_start ? positive_hits << hit : negative_hits << hit }
          end
          
          # Call the hit clusters based on separation on the query sequence.  Use the average length of the hit to set the scale.
          average_length = hits.inject(0) { |sum, hit| sum + hit.query_length } / hits.size.to_f
          clustered_hits = []
          clustered_hits += cluster(positive_hits, cutoff: average_length * 10, cluster_on: :query) if positive_hits.any?
          clustered_hits += cluster(negative_hits, cutoff: average_length * 10, cluster_on: :query) if negative_hits.any?
          
          clustered_hits
        end
      end
    end
  end
end
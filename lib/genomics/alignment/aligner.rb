module Genomics
  module Alignment
    class Aligner
      class << self
        
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
          options = { cluster_on: :subject, cutoff: 10000 }.merge(options)
          
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
        
        
        # Returns a set of clusters of hits, where each cluster can be taken as a single discontiguous entity.  Clusters are determined
        # by hits that are in proximity of each other on the same strand on the same sequence.
        #
        # * *Args*    :
        #   - +hits+ -> The Hit objects to be clustered together, which are assumed to be sorted in increasing position.
        #   - +cutoff+ -> A numeric value setting the maximum distance that two hits can be separate and still be called in the same cluster.
        # * *Options*    :
        #   - +cluster_on+ -> A symbol specifying which sequence object to position based cluster on.  This can be ither the :subject, or 
        #   the :query. (Default :subject)
        # * *Returns* :
        #   - An array of arrays of clustered hits determined by the cutoff.
        #
        def cluster_hits(hits, options = {})
          options = { cluster_on: :subject }.merge(options)
          cluster_on = options[:cluster_on]
          
          # Separate the hits based on their orientations
          categorized_hits = { forward: {}, reverse: {} }
          sort_symbol = cluster_on == :subject ? :subject_start : :query_start
          hits.sort_by(&sort_symbol).each do |hit|
            strand = hit.forward_strand?(on: cluster_on) ? :forward : :reverse
            categorized_subject_hits = categorized_hits[strand][hit.subject] ||= []
            categorized_subject_hits << hit
          end
          
          # Get the average length, which will determine the cutoff
          average_length = hits.inject(0) { |sum, hit| sum + hit.length(on: cluster_on) } / hits.size.to_f

          # Call the hit clusters based on separation.
          clustered_hits = []
          categorized_hits.each do |strand, query_hits|
            query_hits.each do |query, hits|
              clustered_hits += cluster(hits, cutoff: average_length * 10, cluster_on: options[:cluster_on]) if hits.any?
            end
          end
          
          clustered_hits
        end
      end
    end
  end
end
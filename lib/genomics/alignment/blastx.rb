module Genomics
  module Alignment
    # This class handles conducting a blastx alignment of peptides onto a geome.
    class BLASTX < Aligner
      class << self
        # Reads in the results from the provided file and transforms them into the specified format.
        #
        def transform(filename, options = {})
          options = { format: :gff3, output_file: "#{filename}.gff3" }.merge(options)
          
          # Parse the file to get all of the alignment hits
          aggregated_hits = []
          IO::BLASTFormat.open(filename) do |f|
            aggregated_hits = f.entries(aggregate_hits: true)
          end
          
          # Create the list of hits to be converted into single entries based on query, subject, strand, and separation
          # i.e hits widely separated will be treated as distinct matches to the query.
          pbar = ProgressBar.new("Clustering", aggregated_hits.values.inject(0) { |sum, subject_hash| sum + subject_hash.values.length }, STDOUT)

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
          
          # Generate the entries
          pbar = ProgressBar.new("Generating Entries", clustered_hits.size, STDOUT)
          
          entries = []
          clustered_hits.each do |hits|
            pbar.inc
            
            # Initialize the entry
            query, subject = hits.first.query, hits.first.subject
            entry = IO::GFF::Entry.new(seqid: query, source: 'BLASTX', type: :match, attributes: { 'Name' => subject })
            
            # Determine the orientation based on the query start/end positions
            entry.strand = hits.first.query_start < hits.first.query_end ? '+' : '-'
            
            # Add each of the hits to the entry
            hits.each do |hit|
              entry.regions.create(start: hit.query_start, 
                                   end: hit.query_end, 
                                   score: hit.bit_score, 
                                   attributes: { 'Score' => hit.e_value, 'Target' => "#{subject} #{hit.subject_start} #{hit.subject_end}" })
            end
            
            entries << entry
          end
          
          # Write the file
          puts "Writing entries to file..."
          IO::GFFFormat.open(options[:output_file], 'w') do |f|
            f.puts_header
            f.puts(entries)
          end
          
          true
        end
        
        private
        
        # Returns a set of clusters of hits that can be taken as a single discontiguous entry.
        #
        # * *Args*    :
        #   - +hits+ -> The Hit objects to be clustered together, which are assumed to be sorted in increasing position.
        #   - +cutoff+ -> A numeric value setting the maximum distance that two hits can be separate and still be called in the same cluster.
        # * *Returns* :
        #   - An array of arrays of clustered hits determined by the cutoff.
        #
        def cluster(hits, cutoff)
          clusters = []
          running_hits = []
          hits.each do |hit|
            if running_hits.empty? || hit.query_start - running_hits.last.query_end < cutoff
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
      end
    end
  end
end
module Genomics
  module Alignment
    # This class handles conducting an alignment of ESTs onto a genome.
    class EST < Aligner
      class << self
        # Reads in the results from the provided file and transforms them into the specified format.
        #
        # * *Args*    :
        #   - +filename+ -> A string determining the file to be read and transformed into another filetype.
        # * *Returns* :
        #   - A boolean describing whether or not the transformation was successful.
        #
        def transform(filename, options = {})
          options = { format: :gff3, output_file: "#{filename}.gff3", threads: 1 }.merge(options)
          
          # Pull out all of the hits clustered by the query into an array
          query_hits = []
          IO::BLASTFormat.open(filename) do |f|
            f.each_query { |query, hits| query_hits << hits }
          end
          
          # Create a status bar to monitor thread process
          num_queries = query_hits.size
          pbar = ProgressBar.new("Converting Hits", num_queries, STDOUT)
          
          # Create the number of threads specified and split up the array to be processed
          threads = []
          queries_for_thread = (num_queries / options[:threads].to_f).ceil
          query_hits.each_slice(queries_for_thread) do |thread_query_hits|
            threads << Thread.new do
              thread_query_hits.map do |hits|
                pbar.inc
                
                # Cluster the hits
                clusters = cluster_hits(hits, cluster_on: :subject)

                # Convert the clusters to entries
                clusters.map { |clustered_hits| create_entry(clustered_hits) }
              end
            end
          end
          
          # Merge all of the entries together
          entries = threads.map(&:value).flatten
          
          # entries = []
          # IO::BLASTFormat.open(filename) do |f|
          #   # Turn each query into one or more entries
          #   f.each_query do |query, hits|
          #     # # Break the hits up by contig and then cluster
          #     # contig_hits_hash = {}
          #     # hits.each do |hit|
          #     #   # Grows huge!!!
          #     #   contig_hits = contig_hits_hash[hit.subject] ||= []
          #     #   contig_hits << hit
          #     # end
          #     
          #     clusters = cluster_hits(hits, cluster_on: :subject)
          #     # contig_hits_hash.each do |contig, contig_hits| 
          #     #   clusters += cluster_hits(contig_hits, cluster_on: :subject)
          #     # end
          # 
          #     # Convert the clusters to entries
          #     entries += clusters.map { |clustered_hits| create_entry(clustered_hits) }
          #   end
          # end

          # Write the file
          puts "Writing entries to file..."
          IO::GFFFormat.open(options[:output_file], 'w') do |f|
            f.puts_header
            f.puts(entries, progress_bar: true, id_prefix: 'EST')
          end
          
          true
        end
      
        private
      
        # Takes an array of hits and creates a GFF entry out of them
        #
        # * *Args*    :
        #   - +hit_clusters+ -> An array of BLAST::Hits to source the entry.
        # * *Returns* :
        #   - A GFF::Entry representing the collection of hits.
        #
        def create_entry(hits)
          # Initialize the entry
          query, subject = hits.first.query, hits.first.subject
          entry = IO::GFF::Entry.new(seqid: subject, source: 'BLATN', type: :EST_match, attributes: { 'Name' => query })

          # Determine the orientation based on the query start/end positions
          entry.strand = hits.first.subject_start < hits.first.subject_end ? '+' : '-'

          # Add each of the hits to the entry
          hits.each do |hit|
            entry.regions.create(start: hit.subject_start, 
                                 end: hit.subject_end, 
                                 score: hit.bit_score, 
                                 attributes: { 'EValue' => hit.e_value, 'Target' => "#{query} #{hit.query_start} #{hit.query_end}" })
          end
        
          entry
        end
      end
    end
  end
end
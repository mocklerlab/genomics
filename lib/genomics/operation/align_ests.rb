require 'pathname'

module Genomics
  module Operation
    # This class handles the operation logic of aligning ESTs to a genome.
    class AlignESTs
      class << self
        
        # This runs the entire process of generating orthologs from two proteomes.
        #
        # * *Args*    :
        #   - +est_file+ -> A String specifying the path of the FASTA file with the EST sequences.
        #   - +genome_file+ -> A String specifying the path of the FASTA file with the genome sequence.
        # * *Returns* :
        #   - Returns the number of orthologs identified.
        #
        def perform(est_file, genome_file, options = {})
          options = { alignment: {}, clustering: {} }.merge(options)
          
          # Create the alignment files
          alignment_file = generate_alignment(est_file, genome_file, options[:alignment])
        
          # Pass the alignment file to the cluster identification
          identify_est_clusters(alignment_file, options[:clustering])
        end
        
        # This a FASTA formatted set of ESTs and a FASTA formatted genome, aligning the ESTs to the genome.
        #
        # * *Args*    :
        #   - +est_file+ -> An String representing the path to the EST file to be used.
        # * *Options* :
        #   - +blat_path+ -> A String giving the path to blat.
        #   - +blat_options+ -> A Hash of options to be passed directly to the Alignment::BLAT object used to conduct the alignment.
        #   - +alignment_file_dir+ -> A String specifying the directory where the alignment files should be written.  (Default: 
        #   temporary directory)
        # * *Returns* :
        #   - An Array of file paths containing the alignment results 
        #
        def generate_alignment(est_file, genome_file, options = {})
          options = { blat_path: :blat, blat_options: {} }.merge(options)
          options[:blat_options] = { }.merge(options[:blat_options])
          options[:blat_options][:out_format] = :tab
          options[:blat_options][:database] = genome_file
          
          # Prepare the BLAT alignment
          blat = Alignment::BLAT.new(options[:blat_options])
          
          # Optionally set a permanent file to write the results to
          if options[:alignment_file_dir]
            filename = "#{Pathname.new(est_file).basename}_vs_#{Pathname.new(genome_file).basename}"
            blat.output_file = File.join(options[:alignment_file_dir], filename) 
          end
          
          # Run
          result_file = blat.run(est_file)
          result_file.path
        end
        
        # Reads the the alignment file provided and identifies groups of alignments that compose a single EST match againt the genome
        # accounting for the presence of introns.
        #
        # * *Args*    :
        #   - +filename+ -> A string specifying the alignment file to be parsed.
        # * *Options*  :
        #   - +threads+ -> An integer number of threads to use for the clustering (Default: 1).
        #   - +id_prefix+ -> A String that will be used as a prefix for the ID field in the resultant GFF3 file (Default: 'EST'). 
        #   - +e_value+ -> A float specifying the minimum EValue that a hit must have in order to be clustered (Default: 1e^-5).
        # * *Returns* :
        #   - A boolean describing whether or not the transformation was successful.
        #
        def identify_est_clusters(filename, options = {})
          options = { format: :gff3, output_file: "#{filename}.gff3", threads: 1, id_prefix: 'EST', e_value: 1e-5 }.merge(options)
          
          # Pull out all of the hits clustered by the query into an array
          query_hits = []
          IO::BLASTFormat.open(filename) do |f|
            f.each_query do |query, hits| 
              query_hits << hits.select { |hit| hit.e_value <= options[:e_value] }
            end
          end
          
          # Create a status bar to monitor thread process
          pbar = ProgressBar.new("Converting Hits", query_hits.size, STDOUT)
          
          # Thread the creation of the entries from the hits
          entries = Utilities::Threader.thread(query_hits, threads: options[:threads]) do |thread_query_hits|
            thread_query_hits.map do |hits|
              pbar.inc
            
              # Cluster the hits
              clusters = Alignment::Aligner.cluster_hits(hits, cluster_on: :subject)
      
              # Convert the clusters to entries
              clusters.map { |clustered_hits| create_entry(clustered_hits) }
            end
          end.flatten

          # Write the file
          puts "Writing entries to file..."
          IO::GFFFormat.open(options[:output_file], mode: 'w') do |f|
            f.puts_header
            f.puts(entries, progress_bar: true, id_prefix: options[:id_prefix])
          end
          
          true
        end    
        
        private
    
        # Takes an array of hits and creates a GFF entry out of them.
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
require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

module Genomics
  module Operation
    describe TranscriptAligner do
      let(:alignment_file) { File.join(SPEC_PATH, 'fixtures', 'operation', 'est_results.tab') }
  
      describe '#perform' do
        context 'for task :cluster' do
          let(:transcript_aligner) { TranscriptAligner.new(task: :cluster, alignment_file: alignment_file)  }

          it "should raise an error if missing an alignment file" do
            transcript_aligner.alignment_file = nil
            expect { transcript_aligner.perform }.to raise_error(/Missing BLAT alignment file/)
          end
          
          it "should raise an error for an invalid file" do
            transcript_aligner.alignment_file = 'invalid_file'
            expect { transcript_aligner.perform }.to raise_error(/No such file or directory/)
          end
          
          it "should wite the results to a gff3 file and return the path" do
            results_file = transcript_aligner.perform
            File.exists?(results_file).should be(true)
            File.unlink(results_file)
          end
          
          it "should assign IDs and names to the entries" do
            results_file = transcript_aligner.perform

            IO::GFFFormat.open(results_file) do |f|
              f.each do |feature|
                feature.name.should_not be_nil
                feature.id.should_not be_nil
              end
            end

            # File.unlink(results_file)
          end
        end
      end
  
      describe "#identify_est_clusters" do
        # TODO: Implement this after a GFF parser is written, otherwise it it too tedious.
        # it "should treat EST matches distantly separated as different entries" do
        #   BLASTX.identify_est_clusters(results_file)
        #   
        #   File.open("#{results_file}.gff3") do |f|
        #     f.each_line do |line|
        #       next if line =~ /^#/
        #       values = line.split("\t")
        #     
        #       if values[]
        #       values[8].should match(IO::GFFFormat::GFF3_ID_REGEX)
        #       values[8].should match(IO::GFFFormat::GFF3_NAME_REGEX)
        #     end
        #   end
        #   
        #   File.unlink("#{results_file}.gff3")
        # end
      end
    end
  end
end
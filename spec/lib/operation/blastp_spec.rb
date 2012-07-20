require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

module Genomics
  module Operation
    describe BLASTP do
      let(:alignment_file) { File.join(SPEC_PATH, 'fixtures', 'operation', 'blastp', 'blastp_alignment.tab') }
  
      describe '#perform' do
        context "for task :cluster" do
          let(:blastp) { BLASTP.new(task: :cluster, alignment_file: alignment_file)  }
          
          it "should raise an error if missing the alignment_file file" do
            blastp.alignment_file = nil
            expect { blastp.perform }.to raise_error(/Missing BLASTP alignment file/)
          end
          
          it "should raise an error for an invalid file" do
            blastp.alignment_file = 'invalid_file'
            expect { blastp.perform }.to raise_error(/No such file or directory/)
          end

          it "should wite the results to a gff3 file" do
            results_file = blastp.perform
            File.exists?(results_file).should be(true)
            File.unlink(results_file)
          end

          it "should assign IDs and names to the entries" do
            results_file = blastx_aligner.perform

            IO::GFFFormat.open(results_file) do |f|
              f.each do |feature|
                feature.id.should match(/^blastx_/)
                feature.name.should_not be_nil
              end
            end

            File.unlink(results_file)
          end
        end
      end
  
    end
  end
end
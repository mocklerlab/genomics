require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

module Genomics
  module Alignment
    describe BLAST do
      let(:database) { File.join(SPEC_PATH, 'fixtures', 'Athaliana.fa') }
      let(:blast) { BLAST.new('blastp', database: database) }
  
      context '#instance_methods' do
        describe '#run' do
          it "should run the BLAST alignment for the program specified" do
            results_file = blast.run('AAAAA', verbose: true)
            results_file.should be_a(Tempfile)
          end

          it "should write the alignment results to the returned file" do
            blast.out_format = :tab
            results_file = blast.run('AAAAA')
            results_file.read.should_not be(nil)
          end
        end
      end
    end
  end
end
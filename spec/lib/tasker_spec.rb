require File.expand_path(File.dirname(__FILE__) + '/../spec_helper')

module Genomics
  module RBB
    describe Tasker do
      let(:query_alignment_file) { File.join(SPEC_PATH, 'fixtures', 'query_alignment.tab') }
      let(:target_alignment_file) { File.join(SPEC_PATH, 'fixtures', 'target_alignment.tab') }
  
      describe "identify" do
        it "should raise an error if either file doesn't exist" do
          expect { Tasker.identify(query_alignment_file, 'invalid_file') }.to raise_error(/cannot be found/)
          expect { Tasker.identify('invalid_file', target_alignment_file) }.to raise_error(/cannot be found/)
          expect { Tasker.identify('invalid_file', 'invalid_file') }.to raise_error(/cannot be found/)
        end
        
        it "should return the number of reciprocal best blast results identified" do
          Tasker.identify(query_alignment_file, target_alignment_file).should eq(3)
        end
        
        it "should print the results to a file" do
          Tasker.identify(query_alignment_file, target_alignment_file)
          
          orthologs = []
          File.open('orthologs.tab') do |f|
            f.each_line do |line|
              query, target = line.split
              orthologs << [query, target]
            end
          end
          File.unlink('orthologs.tab')
          
          orthologs.should have(3).things
          orthologs.should include(['AT2G15240.1', 'Bradi4g08230.1'])
          orthologs.should include(['AT2G15240.1', 'Bradi4g08230.2'])
          orthologs.should include(['ATCG00670.1', 'Bradi1g05750.1'])
        end
      end
    end
  end
end
require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

module Genomics
  module Alignment
    describe RBB do
      let(:query_alignment_file) { File.join(SPEC_PATH, 'fixtures', 'alignment', 'RBB', 'Athaliana_vs_Spolyrrhiza.tab') }
      let(:target_alignment_file) { File.join(SPEC_PATH, 'fixtures', 'alignment', 'RBB', 'Spolyrrhiza_vs_Athaliana.tab') }
  
      describe "identify" do
        it "should raise an error if either file doesn't exist" do
          expect { RBB.identify([query_alignment_file, 'invalid_file']) }.to raise_error(/No such file or directory/)
          expect { RBB.identify(['invalid_file', target_alignment_file]) }.to raise_error(/No such file or directory/)
          expect { RBB.identify(['invalid_file', 'invalid_file']) }.to raise_error(/No such file or directory/)
        end
        
        it "should return the number of reciprocal best blast results identified" do
          RBB.identify([query_alignment_file, target_alignment_file]).should eq(1)
        end
        
        it "should print the results to a file" do
          RBB.identify([query_alignment_file, target_alignment_file])
          
          orthologs = []
          File.open('orthologs.tab') do |f|
            f.each_line do |line|
              query, target = line.split
              orthologs << [query, target]
            end
          end
          File.unlink('orthologs.tab')
          
          orthologs.should have(1).things
          orthologs.should include(['AT1G01010.1', 'Spipo0007S24960.1'])
        end
      end
    end
  end
end
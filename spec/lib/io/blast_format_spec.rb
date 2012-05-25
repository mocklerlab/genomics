require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

module Genomics
  module IO
    describe BLASTFormat do
      
      context "instance_methods" do
        let(:blast_file_path) { File.join(SPEC_PATH, 'fixtures', 'io', 'blast_format.tab') }
        let(:blast_format) { BLASTFormat.open(blast_file_path) }
        
        describe '#each' do
          it "should iterate through each of the hits" do
            blast_format.each do |hit|
              hit.should be_a(BLAST::Hit)
            end
          end
          
          it "should instantiate each hit with the correct values" do
            blast_format.each do |hit|
              if hit.query == 'scaffold00001' && hit.subject == 'AT1G10170.1'
                hit.percentage_identity.should eq(61.62)
                hit.alignment_length.should eq(1011)
                hit.mismatches.should eq(320)
                hit.gap_openings.should eq(15)
                hit.query_start.should eq(1780005)
                hit.query_end.should eq(1777093)
                hit.subject_start.should eq(207)
                hit.subject_end.should eq(1189)
                hit.e_value.to_f.should eq(0.0)
                hit.bit_score.should eq(1258)
              end
            end
          end
        end
      end
    end
  end
end
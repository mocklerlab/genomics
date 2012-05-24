require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

module Genomics
  module Alignment
    describe FileParser do
      let(:results_file) { File.join(SPEC_PATH, 'fixtures', 'alignment', 'blastx_results.tab') }
  
      describe "#parse" do
        it "should raise an error for an invalid file" do
          expect { FileParser.parse_file('invalid_file') }.to raise_error(/cannot be found/)
        end
        
        it "should return a list of hits" do
          FileParser.parse_file(results_file).each do |hit|
            hit.should be_a(Hit)
          end
        end
        
        it "should return hits with correct values" do
          hits = FileParser.parse_file(results_file)
          
          hit = hits.first
          hit.query.should eq('scaffold00001')
          hit.subject.should eq('AT1G10170.1')
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
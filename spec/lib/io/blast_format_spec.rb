require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

module Genomics
  module IO
    describe BLASTFormat do
      
      context "instance_methods" do
        let(:blast_file_path) { File.join(SPEC_PATH, 'fixtures', 'io', 'blast_format.tab') }
        let(:blast_xml_file_path) { File.join(SPEC_PATH, 'fixtures', 'io', 'blast_format.xml') }
        
        describe '#each' do
          context 'for a tab file' do
            let(:blast_format) { BLASTFormat.open(blast_file_path) }
            
            it "should iterate through each of the hits" do
              blast_format.each do |hit|
                hit.should be_a(BLAST::Hit)
              end
            end

            it "should instantiate each hit with the correct values" do
              blast_format.each do |hit|
                if hit.query == 'scaffold00001' && hit.subject == 'AT1G10170.1'
                  hit.identities.should eq(623)
                  hit.alignment_length.should eq(1011)
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
          
          context 'for an xml file' do
            let(:blast_format) { BLASTFormat.open(blast_xml_file_path, format: :xml) }
            
            it "should iterate through each of the hits" do
              blast_format.each do |hit|
                hit.should be_a(BLAST::Hit)
              end
            end

            it "should instantiate each hit with the correct values" do
              blast_format.each do |hit|
                if hit.query == 'AT1G01010.1' && hit.subject == 'Spipo0001S01170.1'
                  hit.identities.should eq(12)
                  hit.positives.should eq(18)
                  hit.alignment_length.should eq(39)
                  hit.gap_openings.should eq(0)
                  hit.query_start.should eq(56)
                  hit.query_end.should eq(94)
                  hit.subject_start.should eq(195)
                  hit.subject_end.should eq(233)
                  hit.e_value.to_f.should eq(3.19656)
                  hit.bit_score.should eq(22.7126)
                  hit.query_frame.should eq(0)
                  hit.subject_frame.should eq(0)
                  hit.query_sequence.should eq('YKSRDAMWYFFSRRENNKGNRQSRTTVSGKWKLTGESVE')
                  hit.subject_sequence.should eq('FKAGDAALRFFSWARRRPGFRHTTATYNAMLYMAGEARE')
                  hit.midline.should eq('+K+ DA   FFS      G R +  T +    + GE+ E')
                end
              end
            end
          end
        end
        
        describe '#each_query' do
          it "should iterate through each query, hits pair" do
            blast_format.each_query do |query, hits|
              query.should be_a(String)
              hits.should be_an(Array)
              hits.each do |hit|
                hits.should be_a(BLAST::Hit)
              end
            end
          end
        end
      end
    end
  end
end
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
          let(:blast_format) { BLASTFormat.open(blast_file_path) }
          
          it "should iterate through each query, hits pair" do
            debugger
            blast_format.each_query do |query, hits|
              query.should be_a(String)
              hits.should be_an(Array)
              hits.each do |hit|
                hit.should be_a(BLAST::Hit)
              end
            end
          end
        end
        
        describe '#each_cluster' do
          let(:blast_format) { BLASTFormat.open(blast_xml_file_path, format: :xml) }
          
          it "should iterate through each query, hits pair" do
            blast_format.each_cluster do |query, subject, hits|
              query.should be_a(String)
              subject.should be_a(String)
              hits.should be_an(Array)
              hits.each do |hit|
                hit.should be_a(BLAST::Hit)
              end
            end
          end
        end
        
        describe '#query_hits' do
          let(:blast_format) { BLASTFormat.open(blast_xml_file_path, format: :xml) }
          
          it "should return a hash of hits" do
            blast_format.query_hits.should be_a(Hash)
            blast_format.query_hits.each do |query, hits|
              query.should be_a(String)
              hits.first.should be_a(BLAST::Hit)
            end
          end
        end
        
        describe '#clustered_hits' do
          let(:blast_format) { BLASTFormat.open(blast_xml_file_path, format: :xml) }
          
          it "should return a hash of hits" do
            blast_format.clustered_hits.should be_a(Hash)
            blast_format.clustered_hits.each do |query, subject_hash|
              query.should be_a(String)
              subject_hash.each do |subject, clusters|
                subject.should be_a(String)
                clusters.should be_a(Array)
                clusters.first.should be_a(Array)
              end
            end
          end
        end
        
        describe '#hits' do
          let(:blast_format) { BLASTFormat.open(blast_xml_file_path, format: :xml) }
          
          it "should return an Array of Hits" do
            blast_format.hits.should be_an(Array)
            blast_format.hits.each do |entry|
              entry.should be_a(BLAST::Hit)
            end
          end
          
          context 'with the :sort option' do
            it "should return an Array of Hits sorted by bit score" do
              entries = blast_format.entries(sort: true)
              last_bit_score = entries.first.bit_score
              entries.each do |entry|
                entry.bit_score.should be <= last_bit_score
                last_bit_score = entry.bit_score
              end
            end
          end
          
          context 'with the :aggregate option' do
            it "should return a Hash of Hits" do
              blast_format.entries(aggregate: true).should be_an(Hash)
              blast_format.entries(aggregate: true).each do |query, subject_hash|
                query.should be_a(String)
                subject_hash.each do |subject, hits|
                  subject.should be_a(string)
                  hits.first.should be_a(BLAT::Hit)
                end
              end
            end
            
            context 'with the :sort option' do
              it "should return a Hash of sorted Hits" do
                entries = blast_format.entries(aggregate: true, sort: true)
                entries.each do |query, subjects|
                  last_bit_score = subjects.values.first.first.bit_score

                  subjects.each do |subject, hits|
                    # Ensure that the subjects are inserted based on their maximal bit scores
                    last_hit_bit_score = hits.first.bit_score
                    last_hit_bit_score.should be <= last_bit_score
                    last_bit_score = last_hit_bit_score

                    # Check that the hits are sorted
                    hits.each do |hit|
                      hit.bit_score.should be <= last_hit_bit_score
                      last_hit_bit_score = hit.bit_score
                    end
                  end
                end
              end
            end
          end
        end
      end
    end
  end
end
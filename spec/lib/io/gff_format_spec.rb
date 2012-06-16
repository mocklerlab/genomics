require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

module Genomics
  module IO
    describe GFFFormat do
      # let(:genome_file_path) { File.join(SPEC_PATH, 'fixtures', 'io', 'genome.gff3') }
      let(:genome_file_path) { File.join(SPEC_PATH, 'fixtures', 'io', 'augustus_predictions.gff3') }
      let(:est_matches_file_path) { File.join(SPEC_PATH, 'fixtures', 'io', 'est_matches.gff3') }
      
      context 'instance_methods' do
        describe '#each' do
          shared_examples_for 'any gff iterator' do
            it "should iterate through each of the GFF entries" do
              gff_format.each do |feature|
                feature.should be_a(GFF::Feature)
              end
            end
          end
          
          context 'for EST_matches' do
            let(:gff_format) { GFFFormat.open(est_matches_file_path) }
            
            it_should_behave_like 'any gff iterator'
            
            it "should only have entries without features" do
              gff_format.each do |feature|
                feature.should have(0).features
              end
            end

            it "should instantiate features with the correct attributes" do
              gff_format.each do |feature|
                if feature.name == 'FL8T4MU02JLNKH' 
                  feature.id.should eq('EST9')
                  feature.should have(3).regions
                end
              end
            end
          end
          
          context 'for genome' do
            let(:gff_format) { GFFFormat.open(genome_file_path) }
            
            it_should_behave_like 'any gff iterator'

            it "should have entries with features" do
              gff_format.each do |feature|
                debugger
                if feature.name == "AT1G01020"
                  feature.should have(2).features
                end
              end
            end
            
            it "should instantiate features with the correct attributes" do
              gff_format.each do |feature|
                if feature.name == 'AT1G01020' 
                  feature.attributes[:Note].should_not be_nil
                end
              end
            end
          end
        end
        
        describe '#puts' do
          let(:gff_format) { GFFFormat.open(genome_file_path) }
          
          it "should write entries in valid gff3 format" do
            GFFFormat.open('test.gff3', mode: 'w') do |f|
              gff_format.each do |feature|
                f.puts feature
              end
            end
          end
        end
      end
    end
    
    module GFF
      describe EntryRegions do
        let(:entry_regions) { EntryRegions.new(nil) }
        

        
        # describe "#create" do
        #   it "should return an EntryRegion object" do
        #     entry_regions.create(start: 100, end: 500).should be_a(EntryRegion)
        #   end
        #   
        #   it "should have the correct start and end positions" do
        #     entry_region = entry_regions.create(start: 100, end: 500)
        #     entry_region.start.should eq(100)
        #     entry_region.stop.should eq(500)
        #   end
        # end
      end
    end
  end
end
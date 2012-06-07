require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

module Genomics
  module IO
    describe GFFFormat do
      let(:gff_file_path) { File.join(SPEC_PATH, 'fixtures', 'io', 'gff_format.gff3') }
      
      context 'instance_methods' do
        let(:gff_format) { GFFFormat.open(gff_file_path) }
        
        describe '#each' do
          it "should iterate through each of the GFF entries" do
            gff_format.each do |feature|
              feature.should be_a(GFF::Feature)
            end
          end
          
          it "should instantiate features with the correct attributes" do
            gff_format.each do |feature|
              if feature.name == 'FL8T4MU02JLNKH' 
                debugger
                feature.id.should eq('EST9')
                feature.should have(2).regions
                feature.attributes[:EValue].should eq('4.40e-19')
                feature.attributes[:Target].should eq('FL8T4MU02JLNKH 412 458')
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
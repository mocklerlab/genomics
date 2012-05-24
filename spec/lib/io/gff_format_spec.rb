require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

module Genomics
  module IO
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
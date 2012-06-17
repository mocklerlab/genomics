require File.expand_path(File.dirname(__FILE__) + '/../../../spec_helper')

module Genomics
  module IO
    module GFF
      describe Feature do
        let(:feature) { FactoryGirl.build(:feature) }
        
        context '#instance_methods' do
          describe '#<=>' do
            it "should use numerical suffixes for sorting" do
              (FactoryGirl.build(:feature, seqid: 'scaffold_10') <=> FactoryGirl.build(:feature, seqid: 'scaffold_9')).should eq(1)
            end
          end
          
          describe '#id=' do
            it "should set the ID attribute in the attributes hash" do
              feature.id = "New ID"
              feature.attributes[:ID].should eq("New ID")
            end
          end
          
        end
      end
      
      describe Features do
        let(:feature) { FactoryGirl.build(:feature) }
        let(:features) { Features.new(feature) }
        
        context '#creation' do
          describe 'create' do
            it "should return a new Feature" do
              features.create(start: 10, stop: 1000).should be_a(Feature)
            end
            
            it "should set the parrent attribute of the child feature" do
              features.create(start: 10, stop: 1000).attributes[:Parent].should eq(feature.id)
            end
          end
        end
      end
    end
  end
end
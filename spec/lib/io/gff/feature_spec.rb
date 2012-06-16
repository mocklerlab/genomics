require File.expand_path(File.dirname(__FILE__) + '/../../../spec_helper')

module Genomics
  module IO
    module GFF
      describe Feature do
        let(:feature) { FactoryGirl.build(:feature) }
  
        context '#instance_methods' do
          describe '#id=' do
            it "should set the ID attribute in the attributes hash" do
              feature.id = "New ID"
              feature.attributes[:ID].should eq("New ID")
            end
          end
          
        end
      end
    end
  end
end
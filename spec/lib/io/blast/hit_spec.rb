require File.expand_path(File.dirname(__FILE__) + '/../../../spec_helper')

module Genomics
  module IO
    module BLAST
      describe Hit do
        let(:hit) { FactoryGirl.build(:hit) }
  
        context '#instance_methods' do
          describe '#transpose' do
            it "should return the hit object with the query and subject associated attributes switched" do
              original_hit = hit.clone
              hit.transpose!
              hit.query.should eq(original_hit.subject)
            end
          end
          
        end
      end
    end
  end
end
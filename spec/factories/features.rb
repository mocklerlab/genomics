FactoryGirl.define do
  factory :feature, class: Genomics::IO::GFF::Feature do
    sequence(:id) { |n| "Feature#{n}" }
  end
  
end
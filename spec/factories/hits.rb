FactoryGirl.define do
  factory :hit, class: Genomics::IO::BLAST::Hit do
    query               "FL8T4MU02I22OH"
    subject             "scaffold00001"
    percentage_identity 98.73
    alignment_length    474
    mismatches          5
    gap_openings        1
    query_start         1
    query_end           473
    subject_start       7166403
    subject_end         7165930
    e_value	            '2.9e-266'
    bit_score           914.0
  end
end
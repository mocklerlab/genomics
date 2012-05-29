Feature: CLI
  In order to more easily perform bioinformatic comutations
  As a CLI
  I want to be able to execute basic commands

  Scenario: Two Proteomes
    When I run `genomics rbb --protein_files=/Users/mcentee/Documents/genomics/spec/fixtures/Athaliana.fa /Users/mcentee/Documents/genomics/spec/fixtures/Spolyrrhiza.fa --database_files=/Users/mcentee/Documents/genomics/spec/fixtures/Athaliana.fa /Users/mcentee/Documents/genomics/spec/fixtures/Spolyrrhiza.fa`
    Then the output should contain "6 Reciprocal Best Alignments Identified"

  Scenario: Two Alignment Files
    When I run `genomics identify --files=/Users/mcentee/Documents/genomics/spec/fixtures/alignment/rbb/Athaliana_vs_Spolyrrhiza.tab /Users/mcentee/Documents/genomics/spec/fixtures/alignment/rbb/Spolyrrhiza_vs_Athaliana.tab`
    Then the output should contain "1 Reciprocal Best Alignments Identified"

  Scenario: BLASTX Alignment File
    When I run `genomics blastx -i=/Users/mcentee/Documents/genomics/spec/fixtures/alignment/blastx_results.tab`
    Then the output should contain "BLASTX GFF3 successfully created."
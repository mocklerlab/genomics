Feature: CLI
  In order to more easily perform bioinformatic comutations
  As a CLI
  I want to be able to execute basic commands

  Scenario: Two Alignment Files
    When I run `genomics identify --files=/Users/mcentee/Documents/dna_rbb/spec/fixtures/query_alignment.tab /Users/mcentee/Documents/dna_rbb/spec/fixtures/target_alignment.tab`
    Then the output should contain "3 Reciprocal Best Alignments Identified"

  Scenario: BLASTX Alignment File
    When I run `genomics blastx -i=/Users/mcentee/Documents/dna_rbb/spec/fixtures/alignment/blastx_results.tab`
    Then the output should contain "BLASTX GFF3 successfully created."
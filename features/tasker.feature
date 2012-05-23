Feature: Tasker
  In order to calculate the reciprocal best alignment hits
  As a CLI
  I want to be able to return results

  Scenario: Two Alignment Files
    When I run `rbb identify --files=/Users/mcentee/Documents/dna_rbb/spec/fixtures/query_alignment.tab /Users/mcentee/Documents/dna_rbb/spec/fixtures/target_alignment.tab`
    Then the output should contain "3 Reciprocal Best Alignments Identified"
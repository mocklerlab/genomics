require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

module Genomics
  module Sequence
    describe CodonTable do
      
      context "instance_methods" do
        let(:codon_table) { Genomics::Sequence::CodonTable[1] }
        
        describe '#[]' do
          it "should return the correct residue" do
            codons = {
              'ttt' => 'F', 'tct' => 'S', 'tat' => 'Y', 'tgt' => 'C',
              'tta' => 'L', 'tca' => 'S', 'taa' => '*', 'tga' => '*',
            }
            
            codons.each do |codon, residue|
              codon_table[codon].should eq(residue)
            end
          end
          
          it "should return the residue accounting for ambiguity codes" do
            ambiguous_codons = { 'ccn' => 'P', 'ctn' => 'L', 'aar' => 'K', 'agy' => 'S' }
            
            ambiguous_codons.each do |codon, residue|
              codon_table[codon].should eq(residue)
            end            
          end
          
          it "should return nil if the codon doesn't match a single residue" do
            unknown_codons = %w{aam ann rac crr}
            
            unknown_codons.each do |codon, residue|
              codon_table[codon].should eq(nil)
            end            
          end
        end
      end
    end
  end
end
module Genomics
  module Sequence
    # This class wraps bioruby's codon table and provides additional functionality such as ambiguous codon identification.
    class CodonTable < Bio::CodonTable
      
      # Translate a codon into a relevant amino acid.  This method is used for translating a DNA sequence into 
      # amino acid sequence.  If the codon contains an ambiguity code, it attempts to accound for that in the 
      # mapping.
      #
      # * *Args*    :
      #   - +codon+ -> A String representing the codon to be mapped to an amino acid.
      # * *Returns* :
      #   - A String
      #
      def [](codon)
        # Expand any ambiguity codes
        expanded_codons = if codon =~ /[rykwmsbdhvn]/
          codon.split('').inject(['']) do |expanded_bases, base| 
            expanded_bases.product(NA::IUPAC_CODES[base] || [base]).map(&:join)
          end
        else [codon]
        end

        # Get the residues for the codons
        residues = expanded_codons.map do |expanded_codon|
          @table[expanded_codon]
        end
        
        residues.uniq.length == 1 ? residues.first : nil
      end
            
    end
  end
end
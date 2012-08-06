module Genomics
  module Sequence
    # TODO: Respect case when storing sequences, since this can have meaning such as masking, etc.
    class NA < Bio::Sequence::NA
      # A mapping of the 11 standard IUPAC nucleic acid ambiguity codes.
      IUPAC_CODES = {
        'r' => %w{a g}, 'y' => %w{c t}, 'k' => %w{t g}, 'w' => %w{t a}, 'm' => %w{c a}, 's' => %w{c g},
        'b' => %w{c t g}, 'd' => %w{a t g}, 'h' => %w{a t c}, 'v' => %w{a c g},
        'n' => %w{a c g t}
      }
      
      # Translate into an amino acid sequence.
      #
      # * *Args*    :
      #   - +frame+ -> An Integer indicating the reading frame of the translation (1, 2, 3) 
      #   or negatives for reverse complement (Default 1).
      #   - +table+ -> An Integer specifying the index of the codon table (Default 1).
      #   - +unknown+ -> A String that will be used for codons that cannot be determined in translation (Default 'X').
      # * *Returns* :
      #   -
      #
      def translate(frame = 1, table = 1, unknown = 'X')
        # FIXME: Shamelessly copied for bioruby
        if table.is_a?(Bio::CodonTable)
          ct = table
        else
          ct = Genomics::Sequence::CodonTable[table]
        end
        
        naseq = self.dna
        case frame
        when 1, 2, 3
          from = frame - 1
        when 4, 5, 6
          from = frame - 4
          naseq.complement!
        when -1, -2, -3
          from = -1 - frame
          naseq.complement!
        else
          from = 0
        end
        
        nalen = naseq.length - from
        nalen -= nalen % 3
        aaseq = naseq[from, nalen].gsub(/.{3}/) {|codon| ct[codon] or unknown}
        return Bio::Sequence::AA.new(aaseq)
      end
    end
  end
end
module Genomics
  module Alignment
    # This class represets an e_value, which formats itself properly for printing and limits the floating point effects of extremely
    # small values.  It effectively sub classes a Numeric by defering all numerical computation to the float class.
    class EValue
      attr_accessor :coefficient, :exponent

      def initialize(coefficient, exponent = 0)
        value = coefficient.to_f * 10 ** exponent
        @coefficient, @exponent = ("%.2e" % value).split('e').map(&:to_f)
        @exponent = @exponent.to_i
      end
      
      def <=>(other)
        exponent_comparison = @exponent <=> other.exponent
        return exponent_comparison unless exponent_comparison == 0
        
        @coefficient <=> other.coefficient
      end
      
      def to_s
        "%.2e" % to_f
      end
      
      def to_f
        @coefficient.to_f * 10 ** @exponent
      end
      
      def method_missing(name, *args, &block)
        ret = to_f.send(name, *args, &block)
        ret.is_a?(Numeric) ? EValue.new(ret) : ret
      end
    end
  end
end
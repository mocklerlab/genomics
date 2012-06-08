module Genomics
  module Utilities
    # This provides a flexible framework to split a repetative operation on a list into threads and compile the results.
    class Threader
      class << self
        
        # Takes a list and divides it evenly between the the specified number of threads defering to the provided block
        # to implement the computation in the thread.
        #
        # * *Args*    :
        #   - +list+ -> A list to be split into sub-arrays and processed.
        # * *Options* :
        #   - +threads+ -> An integer specifying the number of threads to use.
        # * *Returns* :
        #   -
        #
        def thread(list, options = {})
          options = { threads: 2 }.merge(options)
          
          # Determine the maximum number of elements for each thread.
          elements_for_thread = (list.size / options[:threads].to_f).ceil
          
          # Create the threads and yeild to the block
          threads = []
          list.each_slice(elements_for_thread) do |sub_list|
            threads << Thread.new { yield sub_list }
          end
          
          # Collected all of the results from the threads and join
          threads.map(&:value)
        end
      end
    end
  end
end
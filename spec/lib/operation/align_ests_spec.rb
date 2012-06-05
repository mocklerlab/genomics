require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

module Genomics
  module Operation
    describe AlignESTs do
      let(:results_file) { File.join(SPEC_PATH, 'fixtures', 'operation', 'est_results.tab') }
  
      describe "#identify_est_clusters" do
        it "should raise an error for an invalid file" do
          expect { AlignESTs.identify_est_clusters('invalid_file') }.to raise_error(/No such file or directory/)
        end
        
        it "should wite the results to a gff3 file" do
          AlignESTs.identify_est_clusters(results_file)
          File.exists?("#{results_file}.gff3").should be(true)
          # File.unlink("#{results_file}.gff3")
        end
        
        it "should assign IDs and names to the entries" do
          AlignESTs.identify_est_clusters(results_file)
          
          File.open("#{results_file}.gff3") do |f|
            f.each_line do |line|
              next if line =~ /^#/
              values = line.split("\t")
        
              values[8].should match(IO::GFFFormat::GFF3_ID_REGEX)
              values[8].should match(IO::GFFFormat::GFF3_NAME_REGEX)
            end
          end
          
          # File.unlink("#{results_file}.gff3")
        end

        # TODO: Implement this after a GFF parser is written, otherwise it it too tedious.
        # it "should treat EST matches distantly separated as different entries" do
        #   BLASTX.identify_est_clusters(results_file)
        #   
        #   File.open("#{results_file}.gff3") do |f|
        #     f.each_line do |line|
        #       next if line =~ /^#/
        #       values = line.split("\t")
        #     
        #       if values[]
        #       values[8].should match(IO::GFFFormat::GFF3_ID_REGEX)
        #       values[8].should match(IO::GFFFormat::GFF3_NAME_REGEX)
        #     end
        #   end
        #   
        #   File.unlink("#{results_file}.gff3")
        # end
      end
    end
  end
end
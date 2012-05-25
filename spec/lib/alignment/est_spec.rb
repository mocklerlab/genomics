require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

module Genomics
  module Alignment
    describe EST do
      let(:results_file) { File.join(SPEC_PATH, 'fixtures', 'alignment', 'est_results.tab') }
  
      describe "#transform" do
        it "should raise an error for an invalid file" do
          expect { EST.transform('invalid_file') }.to raise_error(/No such file or directory/)
        end
        
        it "should wite the results to a gff3 file" do
          EST.transform(results_file)
          File.exists?("#{results_file}.gff3").should be(true)
          # File.unlink("#{results_file}.gff3")
        end
        
        # it "should assign IDs and names to the entries" do
        #   BLASTX.transform(results_file)
        #   
        #   File.open("#{results_file}.gff3") do |f|
        #     f.each_line do |line|
        #       next if line =~ /^#/
        #       values = line.split("\t")
        # 
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
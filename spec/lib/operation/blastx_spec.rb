require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

module Genomics
  module Operation
    describe BLASTX do
      let(:alignment_file) { File.join(SPEC_PATH, 'fixtures', 'operation', 'blastx_results.tab') }
      let(:genome_file) { File.join(SPEC_PATH, 'fixtures', 'operation', 'Spolyrrhiza.fa') }
      let(:proteome_file) { File.join(SPEC_PATH, 'fixtures', 'Athaliana.fa') }
  
      describe '#perform' do
        context 'for task :align' do
          let(:blastx_aligner) { BLASTX.new(task: :align, genome_file: genome_file, proteome_file: proteome_file)  }
          
          it "should raise an error if missing an genome file" do
            blastx_aligner.genome_file = nil
            expect { blastx_aligner.perform }.to raise_error(/Missing genome FASTA file/)
          end
          
          it "should raise an error if missing a proteome file" do
            blastx_aligner.proteome_file = nil
            expect { blastx_aligner.perform }.to raise_error(/Missing proteome FASTA file/)
          end
          
          it "generate a BLAST alignment file" do
            results_file = blastx_aligner.perform
            File.exists?(results_file).should be(true);
            File.unlink(results_file)
          end
        end
        
        context "for task :cluster" do
          let(:blastx_aligner) { BLASTX.new(task: :cluster, alignment_file: alignment_file)  }
          
          it "should raise an error if missing the alignment_file file" do
            blastx_aligner.alignment_file = nil
            expect { blastx_aligner.perform }.to raise_error(/Missing BLASTX alignment file/)
          end
          
          it "should raise an error for an invalid file" do
            blastx_aligner.alignment_file = 'invalid_file'
            expect { blastx_aligner.perform }.to raise_error(/No such file or directory/)
          end

          it "should wite the results to a gff3 file" do
            results_file = blastx_aligner.perform
            File.exists?(results_file).should be(true)
            File.unlink(results_file)
          end

          it "should assign IDs and names to the entries" do
            results_file = blastx_aligner.perform

            IO::GFFFormat.open(results_file) do |f|
              f.each do |feature|
                feature.id.should match(/^blastx_/)
                feature.name.should_not be_nil
              end
            end

            File.unlink(results_file)
          end
        end
      end
  
    end
  end
end
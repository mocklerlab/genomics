require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

module Genomics
  module Operation
    describe RBB do
      let(:query_alignment_file) { File.join(SPEC_PATH, 'fixtures', 'alignment', 'RBB', 'Athaliana_vs_Spolyrrhiza.tab') }
      let(:target_alignment_file) { File.join(SPEC_PATH, 'fixtures', 'alignment', 'RBB', 'Spolyrrhiza_vs_Athaliana.tab') }
      let(:athaliana_proteome) { File.join(SPEC_PATH, 'fixtures', 'Athaliana.fa') }
      let(:spolyrrhiza_proteome) { File.join(SPEC_PATH, 'fixtures', 'Spolyrrhiza.fa') }

      describe '#perform' do
        it "should return the number of orthologs" do
          alignment_files = RBB.perform([{ file: athaliana_proteome, database: athaliana_proteome }, { file: spolyrrhiza_proteome, database: spolyrrhiza_proteome }])
          File.unlink('orthologs.tab')
          alignment_files.should be_a(Integer)
          alignment_files.should eq(6)
        end
      end
  
      describe '#generate_alignments' do
        it "should return the alignment files" do
          alignment_files = RBB.generate_alignments([{ file: athaliana_proteome, database: athaliana_proteome }, { file: spolyrrhiza_proteome, database: spolyrrhiza_proteome }])
          alignment_files.each do |alignment_file|
            alignment_file.should be_a(String)
          end
        end
      end
  
      describe '#identify_orthologs' do
        it "should raise an error if either file doesn't exist" do
          expect { RBB.identify_orthologs([query_alignment_file, 'invalid_file']) }.to raise_error(/No such file or directory/)
          expect { RBB.identify_orthologs(['invalid_file', target_alignment_file]) }.to raise_error(/No such file or directory/)
          expect { RBB.identify_orthologs(['invalid_file', 'invalid_file']) }.to raise_error(/No such file or directory/)
        end
        
        it "should return the number of reciprocal best blast results identified" do
          RBB.identify_orthologs([query_alignment_file, target_alignment_file]).should eq(1)
        end
        
        it "should print the results to a file" do
          RBB.identify_orthologs([query_alignment_file, target_alignment_file])
          
          orthologs = []
          File.open('orthologs.tab') do |f|
            f.each_line do |line|
              query, target = line.split
              orthologs << [query, target]
            end
          end
          File.unlink('orthologs.tab')
          
          orthologs.should have(1).things
          orthologs.should include(['AT1G01010.1', 'Spipo0007S24960.1'])
        end
      end
    end
  end
end
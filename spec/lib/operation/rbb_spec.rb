require File.expand_path(File.dirname(__FILE__) + '/../../spec_helper')

module Genomics
  module Operation
    describe RBB do
      describe '#perform' do
        let(:athaliana_proteome) { File.join(SPEC_PATH, 'fixtures', 'Athaliana.fa') }
        let(:spolyrrhiza_proteome) { File.join(SPEC_PATH, 'fixtures', 'Spolyrrhiza.fa') }
        let(:proteome_files) { [{ file: athaliana_proteome, database: athaliana_proteome }, 
                                { file: spolyrrhiza_proteome, database: spolyrrhiza_proteome }] }
        let(:query_alignment_file) { File.join(SPEC_PATH, 'fixtures', 'alignment', 'RBB', 'Athaliana_vs_Spolyrrhiza.tab') }
        let(:target_alignment_file) { File.join(SPEC_PATH, 'fixtures', 'alignment', 'RBB', 'Spolyrrhiza_vs_Athaliana.tab') }
        
        context 'for task :align' do
          let(:rbb_aligner) { RBB.new(task: :align, proteomes: proteome_files ) }
                                                                
          it "should raise an error if missing protomes files" do
            rbb_aligner.proteomes = []
            expect { rbb_aligner.perform }.to raise_error(/Missing proteome files/)
          end
          
          it "should return the alignment files" do
            alignment_files = rbb_aligner.perform
            alignment_files.each do |alignment_file|
              alignment_file.should be_a(Tempfile)
            end
          end
        end
        
        context 'for task :ortholog' do
          let(:rbb_orthologer) { RBB.new(task: :ortholog, alignment_files: [query_alignment_file, target_alignment_file]) }
          
          it "should raise an error for missing alignment files" do
            expect do 
              rbb_orthologer.alignment_files = []
              rbb_orthologer.perform
            end.to raise_error(/Missing alignment files/)
          end
          
          it "should raise an for invalid files" do
            expect do 
              rbb_orthologer.alignment_files = [query_alignment_file, 'invalid_file']
              rbb_orthologer.perform
            end.to raise_error(/No such file or directory/)
            
            expect do 
              rbb_orthologer.alignment_files = ['invalid_file', target_alignment_file]
              rbb_orthologer.perform
            end.to raise_error(/No such file or directory/)
            
            expect do 
              rbb_orthologer.alignment_files = ['invalid_file', 'invalid_file']
              rbb_orthologer.perform
            end.to raise_error(/No such file or directory/)
          end
          
          it "should return the number of results identified" do
            rbb_orthologer.perform.should eq(1)
          end
          
          it "should print the results to a file" do
            rbb_orthologer.perform

            orthologs = []
            File.open('orthologs.tab') do |f|
              f.each_line do |line|
                orthologs << line.split
              end
            end
            File.unlink('orthologs.tab')

            orthologs.should have(1).things
            orthologs.should include(['AT1G01010.1', 'Spipo0007S24960.1'])
          end

          context 'with the :detailed option' do
            it "should top alignment details" do
              rbb_orthologer.detailed = true
              rbb_orthologer.perform

              begin
                File.open('orthologs.tab') do |f|
                  f.each_line do |line|
                    line.split.should have(4).things
                  end
                end
              ensure
                File.unlink('orthologs.tab')
              end
            end
          end
          
          context 'without the :reciprocal option' do
            it "should top alignment details" do
              rbb_orthologer.reciprocal = false
              rbb_orthologer.perform

              begin
                File.open('orthologs.tab') do |f|
                  f.each_line do |line|
                    line.split.should have(3).things
                  end
                end
              ensure
                File.unlink('orthologs.tab')
              end
            end
          end
        end
        
      end
    end
  end
end
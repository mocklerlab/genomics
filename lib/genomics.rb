require "genomics/version"

# Require external gems
require 'progressbar'
require 'ox'

##### Require internal files #####

require 'genomics/command_line'

### Require the alignment module ###
require 'genomics/alignment/blast'
require 'genomics/alignment/blat'

require 'genomics/alignment/aligner'

### Require the IO module ###
# Require the base IO wrapper first
require 'genomics/io/flat_file_format'

# Require various file formats
require 'genomics/io/blast_format'
require 'genomics/io/gff_format'

# Require BLAST data structures
require 'genomics/io/blast/e_value'
require 'genomics/io/blast/hit'

# Require GFF data sctructures
require 'genomics/io/gff/entry'
require 'genomics/io/gff/feature'
require 'genomics/io/gff/region'

### Require the Operations module ###
require 'genomics/operation/blastp'
require 'genomics/operation/blastx'
require 'genomics/operation/transcript_aligner'
require 'genomics/operation/rbb'

### Require the Utilities module ###
require 'genomics/utilities/threader'

module Genomics
end

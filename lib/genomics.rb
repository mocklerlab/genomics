require "genomics/version"

# Require external gems
require 'progressbar'

##### Require internal files #####

# Require the alignment module
require 'genomics/alignment/blastx'
require 'genomics/alignment/rbb'

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

module Genomics
end

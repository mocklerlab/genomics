require "genomics/version"

# Require external gems
require 'progressbar'

# Require internal files 
require 'genomics/file_parser'

# Require the alignment module
require 'genomics/alignment/blastx'
require 'genomics/alignment/hit'
require 'genomics/alignment/rbb'
require 'genomics/alignment/e_value'
require 'genomics/alignment/file_parser'

# Require the IO module
require 'genomics/io/blast_format'
require 'genomics/io/gff_format'
require 'genomics/io/flat_file_format'
require 'genomics/io/gff/entry'

module Genomics
end

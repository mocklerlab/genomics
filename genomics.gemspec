# -*- encoding: utf-8 -*-
$:.push File.expand_path("../lib", __FILE__)
require "genomics/version"

Gem::Specification.new do |s|
  s.name        = "genomics"
  s.version     = Genomics::VERSION
  s.authors     = ["Connor McEntee"]
  s.email       = ["sallustfire@gmail.com"]
  s.homepage    = ""
  s.summary     = %q{Takes the results of pairwise alignments and determines the reciprocally best alignments.}
  s.description = %q{Use it to parse the results of two alignment files, from which it calculates the reciprocally best alignments.}

  s.rubyforge_project = "genomics"

  s.files         = `git ls-files`.split("\n")
  s.test_files    = `git ls-files -- {test,spec,features}/*`.split("\n")
  s.executables   = `git ls-files -- bin/*`.split("\n").map{ |f| File.basename(f) }
  s.require_paths = ["lib"]
  
  s.add_dependency "thor"
  s.add_dependency "ruby-progressbar"

  s.add_development_dependency "ruby-debug19"
  s.add_development_dependency "rspec"
  s.add_development_dependency "cucumber"
  s.add_development_dependency "aruba"
end
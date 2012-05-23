# -*- encoding: utf-8 -*-
$:.push File.expand_path("../lib", __FILE__)
require "rbb/version"

Gem::Specification.new do |s|
  s.name        = "rbb"
  s.version     = Genomics::RBB::VERSION
  s.authors     = ["Connor McEntee"]
  s.email       = ["sallustfire@gmail.com"]
  s.homepage    = ""
  s.summary     = %q{TODO: Write a gem summary}
  s.description = %q{TODO: Write a gem description}

  s.rubyforge_project = "rbb"

  s.files         = `git ls-files`.split("\n")
  s.test_files    = `git ls-files -- {test,spec,features}/*`.split("\n")
  s.executables   = `git ls-files -- bin/*`.split("\n").map{ |f| File.basename(f) }
  s.require_paths = ["lib"]
  
  s.add_dependency "thor"
  s.add_dependency "dna", ">= 0.4.1"
  
  s.add_development_dependency "rspec"
  s.add_development_dependency "cucumber"
  s.add_development_dependency "aruba"
  # s.add_development_dependency "spork"
end

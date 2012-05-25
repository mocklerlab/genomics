require 'rubygems'
require 'bundler/setup'
require 'spork'
require 'spork/ext/ruby-debug'
require 'rspec'
require 'ruby-prof'

require 'genomics'


Spork.prefork do
  # Set an environment variable denoting the spec path
  SPEC_PATH = File.dirname(__FILE__)
end

Spork.each_run do
  # This code will be run each time you run your specs.
  $".grep(/genomics\/lib/).each {|e| $".delete(e) && require(e) }
end
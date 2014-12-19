# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'kepler/version'

Gem::Specification.new do |spec|
  spec.name          = "kepler"
  spec.version       = Kepler::VERSION
  spec.authors       = ["Jordan Stephens"]
  spec.email         = ["iam@jordanstephens.net"]
  spec.summary       = %q{A ruby gem for working with two-body Keplerian Orbits.}
  spec.description   = %q{A ruby gem for working with two-body Keplerian Orbits.}
  spec.homepage      = ""
  spec.license       = "MIT"

  spec.files         = `git ls-files -z`.split("\x0")
  spec.executables   = spec.files.grep(%r{^bin/}) { |f| File.basename(f) }
  spec.test_files    = spec.files.grep(%r{^(test|spec|features)/})
  spec.require_paths = ["lib"]

  spec.add_development_dependency "rake", "~> 10.0"
  spec.add_development_dependency "rspec", "~> 3.0"
  spec.add_development_dependency "pry-byebug", "~> 2.0"
end

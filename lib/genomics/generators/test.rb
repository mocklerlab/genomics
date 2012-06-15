require 'thor/group'

module Genomics
  module Generators
    class Test < Thor::Group
      argument :group, :type => :string
      argument :name, :type => :string
      include Thor::Actions
      
      def self.source_root
        File.dirname(__FILE__) + "/recipe"
      end
      
      def create_group
        empty_directory(group)
      end
      
      def copy_recipe
        template("recipe.txt", "#{group}/#{name}.txt")
      end
    end
  end
end
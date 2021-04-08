#!/usr/bin/env ruby

require 'optparse'
require 'optparse/time'
require 'ostruct'
require 'pp'

# This script expects a headerless input file with three columns:
# 1: The TF gene name/identifier
# 2: The target gene name/identifier
# 3: The p-value of the interaction between the TF and target gene

class OptionParse
	def self.parse(args)
		options = OpenStruct.new
		options.file_name=""
		options.pval = 1
		opts = OptionParser.new do |opts|
			opts.banner = "Usage: estimate_affinity.rb [options]"
			opts.separator ""
			opts.separator "Specific options:"
			# Mandatory argument.
			opts.on("-i", "--input FILE",
				"") do |file_name|
				options.file_name = file_name
			end
			opts.on("-p", "--pval FLOAT",
				"") do |pval|
				options.pval = pval.to_f
			end
		end
		
		opts.parse!(args)
		options
	end 
end 

def rank_array(arr)
	a_len = arr.length
	rankpos_arr = (0..a_len-1).to_a.sort_by{ |x| arr[x] }
	rank_arr = Array.new(a_len)
	rankpos_arr.each_with_index do |x,i|
		rank_arr[x] = i.to_f / (a_len-1).to_f
	end
	return rank_arr
end

options = OptionParse.parse(ARGV)

$stderr.puts 'options set:'
$stderr.puts '  input file: ' + options.file_name
$stderr.puts '  pval cutoff: ' + options.pval.to_s

targets = Hash.new
file = File.open(options.file_name)
file.each_with_index do |line,i|
    next if line.chomp! =~ /^$|#/
	dline = line.split("\t")
	key = dline[0] + "\t" + dline[1]
	if options.pval > dline[2].to_f
		if targets[key].nil?
			targets[key] = Array.new
		end
		targets[key] << - Math.log10(dline[2].to_f)
	end
end

file.close

$stderr.puts 'read file in'

site_sum = Array.new
site_max = Array.new
site_count = Array.new

$stderr.puts 'number of targets: ' + targets.length.to_s

targets.each do |key, value|
	site_sum << value.inject(:+)
	site_max << value.max
	site_count << value.length
end

$stderr.puts 'processed targets'

rank_site_sum =  rank_array(site_sum)
rank_site_max = rank_array(site_max)
rank_site_count = rank_array(site_count)

index=0
targets.each do |key, value|
	lineout = key + "\t"  + value.inject(:+).to_s + "\t" + rank_site_sum[index].to_s
	lineout = lineout + "\t" + value.max.to_s + "\t" + rank_site_max[index].to_s
	lineout = lineout + "\t" + value.length.to_s + "\t" + rank_site_count[index].to_s
	lineout = lineout + "\t" + (rank_site_sum[index] > rank_site_max[index] ? rank_site_sum[index].to_s : rank_site_max[index].to_s)
	puts lineout
	index = index + 1
end




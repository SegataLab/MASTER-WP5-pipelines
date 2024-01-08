
#count_reads.rb

if ARGV[0]==nil
puts ""
puts ".........................................................................."
puts "... Sorry! you forgot to add the folder name which contains .bam files ..."
puts "...          try again adding this information, please                 ..."
puts ".........................................................................."
puts ""
exit
elsif ARGV[0]=~/\//
`ls "#{ARGV[0]}"* > list.txt`
else 
`ls "#{ARGV[0]}"/* > list.txt`
end

l=''
l << `wc -l list.txt`
l=l.split("\s")[0]
m=0

sample=''
out2=File.new("report_different_match.txt","w")
out3=File.new("report_different_family.txt","w")

hfam={}
ff=File.open("phenotypes.txt").each_line do |line|
line.chomp!
cols=line.split("\t")
hfam[cols[0]]=cols[1]
end
ff.close

genes=[]
hgenes={}
gg=File.open("gene_list.txt").each_line do |line|
line.chomp!
if line =~ /SN:(\S+)/
genes << $1
end
end
gg.close

samples=[]
hsamples={}
ll=File.open("list.txt").each_line do |line|
line.chomp!
if line =~ /\/(.*)\.\w\w\w/
samples << $1
end
end
ll.close

genes.each_index{|a| hgenes[genes[a]]=a}
samples.each_index{|b| hsamples[samples[b]]=b}
matrix=Array.new(samples.length()){|i| Array.new(genes.length()) {|x| 0}}

aa=File.open("list.txt").each_line do |file|
file.chomp!
if file =~ /\/(.*).bam/ or file =~ /\/(.*).sam/ 
sample=$1
m+=1
puts "...processing sample #{sample} (#{m}/#{l})"
puts hsamples[sample]
`samtools sort -n "#{file}" > file1.bam`
`samtools view -h file1.bam > file2.sam`
read=''
        bb=File.open("file2.sam").each_line do |line|
 	line.chomp!
	cols=line.split("\t")
	if cols[0]!=read
	matrix[hsamples[sample]][hgenes[cols[2]]]+=1
        elsif cols[6]!="="
	out2.puts line
		if hfam[cols[2]]!=hfam[cols[6]]
		out3.puts "#{hfam[cols[2]]}\t#{hfam[cols[6]]}\t#{line}"
		matrix[hsamples[sample]][hgenes[cols[2]]]+=1		
		end
	end
   end
   bb.close
end
end
aa.close

matfile=File.new("matrix_Resfinder.txt","w")
mm=-1
genes.each {|i| matfile.print "#{i}\t"}
matfile.print "\n"
matrix.each {|i| i.each {|x| matfile.print "#{x}\t"}}
mm+=1
matfile.print "#{samples[mm]}\n"

puts "---------------------------------------------------------------"
puts ""
puts " Analysis performed correctly, if you have any issue please contact to jcobd@unileon.es"
puts ""


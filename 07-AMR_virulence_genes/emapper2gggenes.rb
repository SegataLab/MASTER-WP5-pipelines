
hcoord={}

aa=File.open("#{ARGV[0]}.emapper.gff").each_line do |line|
line.chomp!
col=line.split("\t")
if line =~ /ID=(\S+)\;em\_target/
#puts $1
    if  col[6]=="+"
    hcoord[$1]="#{col[3]}\t#{col[4]}\t1"
    else hcoord[$1]="#{col[3]}\t#{col[4]}\t0"
    end
end
end
aa.close

bb=File.open("#{ARGV[0]}.emapper.annotations").each_line do |line|
line.chomp!
col=line.split("\t")
if line =~/##/
elsif line =~/#/
#puts "contig\tstart\tend\tstrand\tanotation\tCOG_category\tDescription\tPreferred_name\tPFAMs"
puts "contig\tstart\tend\tstrand\t#{line}"
elsif col[0]=~/(.*)\_\d+/
#puts "#{$1}\t#{hcoord[col[0]]}\t#{col[6]}\t#{col[7]}\t#{col[8]}\t#{col[20]}"
puts "#{$1}\t#{hcoord[col[0]]}\t#{line}"
end
end
bb.close


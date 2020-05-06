#!/usr/bin/perl

while(<>){
	chomp;
	@arr = split /\t/;

	if($arr[2] ne "\*" && $arr[6] ne "\*"){
		print "$_\n";
	}
}

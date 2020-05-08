#!/usr/bin/perl

while(<>){
	chomp;
	@arr = split /\t/;

	if($arr[2] ne "\*"){
		print "$_\n";
	}
}

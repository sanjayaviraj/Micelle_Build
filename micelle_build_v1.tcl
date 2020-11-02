##########################################################
# Script to build a micelle 						     #
# Open the monomer pdb and source this script            #
# one residue id = one monomer 							 # 
# Sanjaya Viraj Bandara| sanjayavirajbandara@gmail.com   #
# 17/09/2020											 #
##########################################################

package require topotools
package require pbctools

set nMol 60 ; #type number of monomers
set SphericalRadius 40; #type sphereical radius

#multifly molecules 
proc multiplyMol {nMol} {
	set sel [atomselect top all]
	set nAtomsRef [molinfo top get numatoms]
	set nAtomsNew [expr {$nAtomsRef * $nMol}] ; # number of atoms in micell
	mol new atoms $nAtomsNew
	animate dup top

	set proplist {name type radius resname mass x y z}
	set atomdata {}

	for {set i 0} {$i < $nMol} {incr i} {
		#chage resname so can be easly identify the single monomers
		$sel set resname $i
		$sel moveby {50 0 0}
		set atomdataUpdate [$sel get $proplist]
		set atomdata [concat $atomdata $atomdataUpdate]
	}
	set selNewMol [atomselect top all]
	$selNewMol set $proplist $atomdata

	topo retypebonds 
	vmdcon -info "assigned [topo numbondtypes] bond types to [topo numbonds] bonds:"
	vmdcon -info "bondtypes: [topo bondtypenames]"
	
	topo guessangles
	vmdcon -info "assigned [topo numangletypes] angle types to [topo numangles] angles:"
	vmdcon -info "angletypes: [topo angletypenames]"
	
	mol bondsrecalc top
	mol reanalyze top
}

proc placeInSphere {nMol SphericalRadius} {
	for {set i 0} {$i < $nMol} {incr i} {
		set selresidue [atomselect 2 "resname $i"]
		#move the molecule to the center
		$selresidue moveby [vecscale -1 [measure center $selresidue]]

		#generate random numbers
		proc rand_range {min max} { 
			return [expr int(rand()*($max-$min+1)) + $min] 
		}
			
		#move on the sphere
		set theta [rand_range 1 360]
		set psi [rand_range 1 360]
		set x [expr $SphericalRadius*sin(double($theta)) * cos(double($psi))]
		set y [expr $SphericalRadius*sin (double($theta)) * sin (double($psi))]
		set z [expr $SphericalRadius * cos (double($theta)) ]

		#orientation
		$selresidue move [transvec "[list $x $y $z]"]
		$selresidue moveby [list $x $y $z]
	}
}
#call multiply process
multiplyMol $nMol
#write psf and pdb files
set selNew [atomselect 1 all]
$selNew writepdb P_$nMol.pdb
$selNew writepsf P_$nMol.psf

mol delete 0
mol delete 1
#load exported pdb and psf files
mol new P_$nMol.pdb autobonds yes waitfor all
mol addfile P_$nMol.psf type {psf} 
#call the place in the sphere process
placeInSphere $nMol $SphericalRadius

topo guessangles
topo guessdihedrals
topo guessimpropers

set selNew [atomselect 2 all]
$selNew writepdb P_$nMol.pdb
$selNew writepsf P_$nMol.psf

puts "DONE"
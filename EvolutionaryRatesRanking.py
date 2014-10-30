#!/usr/bin/python
from __future__ import division

import dendropy
import os
import glob

#########################################################################################################################
#This script was written by Nathan Whelan.  

# THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
# THE CONTRIBUTORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF 
# OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE 
# SOFTWARE.
##########################################################################################################################

###################################################################################################
#This script will create a list of alignments sorted by evolutionary rate
#as measured by tree length divided by taxa. 
#THIS SHOULD BE EXECUTED IN A FOLDER WITH SINGLE GENE TREES. As written works with RAxML_bipartition trees.
#If your alignments had different format than change variable directly below
treeFileFormat='*.fas.out' #This should be extension for gene trees in newick format
slow2Fast=True  ##Change to False if you want fastest first in the list descending to slow
path='' #PATH WITH TREES. As written, directory from which this script is executed
####################################################################################################

#Function that returns evol rate for a dendropy.Tree object
def evolRate(TREE):
	tree_length=TREE.length()
	numberOfTaxa=len(TREE.leaf_nodes())
	return(tree_length/numberOfTaxa)

##Goes through each tree and calculates its evolutionary rate and places it in dictionary ratesDict
ratesDict={} #Create empty dictionary.
for infile in glob.glob(os.path.join(path, treeFileFormat)): 
	myTree=dendropy.Tree(stream=open(infile),schema="newick")
	X=evolRate(myTree)
	ratesDict[infile]=X

##Creates a list from dictionary ratesDict that is sorted by rate. 
##If slow2Fast = True then slowest evolving gene will be first followed by increasingly fast
if slow2Fast == True:
	ratesList=sorted(ratesDict, key=ratesDict.__getitem__) 
else:
	ratesList=sorted(ratesDict, key=ratesDict.__getitem__, reverse=True)

#this will write a file with alignments in order of slowest evolRate to fastest
with open("evolRatesOrdered.txt","w") as evolRatesFile: #Old file will be overwritten if it exists 
	for i in ratesList:	
		y=i.split(".")		#This and the next line makes assumption the trees are RAxML_bipartitions.alignment.out
		evolRatesFile.write(y[1] + "." + y[2] + "\n") #If used RAxML then list will be in format of files used for individual RAxML gene trees

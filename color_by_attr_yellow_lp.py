# Colors and saves movies of a Chimera model from various angles

# Imports
import getopt, sys # Allows for command line arguments
import os
import chimera
from chimera import runCommand as rc

try:
	opts, args = getopt.getopt(sys.argv[1:], '', ['rep=', 'test=', 'qID=', 'rID=', 'lengthID=', 'cutoff=', 'colorType=', 'testType=', 'attr=', 'minVal=', 'maxVal=', 'frameCount='])
except getopt.error as message:
    raise chimera.NonChimeraError("%s: %s" % (__name__, message)) # Prints error if present

rep = ''
test = '' # doesn't actually exist in GUI3 script
qID = '' 
rID = ''
lengthID = '' # amino acid count
cutoff = ''
colorType = '' # rg, yb, or om
testType = '' # flux or corr
attr = '' # delta, pval, or stat
max_value = ''
min_value = ''
frame_count = ''

# Assigns option values to variables
for o,a  in opts:
    if o == '--rep':
        rep = a
    if o == '--test':
        test = a
    if o == '--qID':
        qID = a
    if o == '--rID':
        rID = a
    if o == '--lengthID':
        lengthID = a
    if o == '--cutoff':
        cutoff = a
    if o == '--colorType':
        colorType = a
    if o == '--testType':
        metric = a
    if o == '--attr':
        attr = a
    if o == '--minVal':
        min_value = a
    if o == '--maxVal':
        max_value = a
    if o == '--frameCount':
        frame_count = a

# Ensures that all variables have been assigned and ends script if any are missing
assert(len(rep) != 0)
assert(len(test) != 0)
assert(len(qID) != 0)
assert(len(rID) != 0)
assert(len(lengthID) != 0)
assert(len(cutoff) != 0)
assert(len(colorType) != 0)
assert(len(metric) != 0)
assert(len(attr) != 0)
assert(len(min_value) != 0)
assert(len(max_value) != 0)
assert(len(frame_count) != 0)

model = "1"
pdb_file = "%s.pdb" % (rID)
attr_file = "attribute_files/attr_%s_%s_%s_%s_%s.dat" % (qID, rID, attr, metric, cutoff)
central_residue = int(float(lengthID)/2)

if colorType == "rg":
 min_color = "red"
 max_color = "green"
if colorType == "yb":
 min_color = "yellow"
 max_color = "blue"
if colorType == "om":
 min_color = "orange"
 max_color = "magenta"
if colorType == "bw":
 min_color = "blue"
 max_color = "white"
if colorType == "rw":
 min_color = "red"
 max_color = "white"
if colorType == "wg":
 min_color = "white"
 max_color = "green"
if colorType == "wo":
 min_color = "white"
 max_color = "orange"

mid_color = ""
mid_value = ""
if attr == "delta" or attr == "deltaKL":
 mid_color = "white"
 mid_value = (float(min_value) + float(max_value))/2

rep_key = ""

mol = chimera.openModels.open(pdb_file)[0] # Opens molecule
rc("open md:%s.meta" % (rID))
# Run Preliminary Chimera Commands
rc("windowsize 1500 1080")
rc("lighting reflectivity 0")
rc("lighting brightness .85")
rc("lighting sharpness 10")
rc("center #%s:.A" % (model))
rc("scale 0.4")
rc("~ribbon")
#rc("~surface")
rc("~display")
rc("background solid black")

if rep == "surface":
 rc("surface #%s" % (model))
 rep_key = "s"
if rep == "ribbon":
 rc("ribbon #%s" % (model))
 rep_key = "r"
if rep == "ribbonsurface":
 rc("ribbon #%s" % (model))
 rc("surface #%s" % (model))
 rc("surftransparency 65 #%s" % (model))
 rep_key = "r,s"
rc("defattr %s raiseTool false" % (attr_file))

col_string = "%s %s %s %s %s %s novalue yellow" % (max_value, max_color, mid_value, mid_color, min_value, min_color)

rc("colorkey  .08,.1  .10,.800  %s" % (col_string))

rc("rangecolor %s,%s %s #%s" % (attr, rep_key, col_string, model)) 
rc("rangecolor %s,a %s #%s" % (attr, col_string, model))
# open again with ligand
rc("open md:%s.meta" % (rID))

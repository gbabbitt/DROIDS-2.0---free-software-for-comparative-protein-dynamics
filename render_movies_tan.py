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
if attr == "delta":
 mid_color = "white"
 mid_value = (float(min_value) + float(max_value))/2

out_folder = "Videos" # generalize
#out_folder = "%s/%s_%s_%s_%s_%s" % (out_folder_pre, rID, qID, rep, attr, metric)
rep_key = ""

#if not os.path.exists(out_folder_pre):
#    os.makedirs(out_folder_pre)
if not os.path.exists(out_folder):
    os.makedirs(out_folder)

#########################################################################
####  make 6 fixed point of view movies and 2 rolling movies   ##########
#########################################################################

# open and setup pdb file
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
rc("~surface")
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
 rc("surftransparency 50 #%s" % (model))
 rep_key = "r,s"
rc("defattr %s raiseTool false" % (attr_file))
col_string = "%s %s %s %s %s %s novalue tan" % (max_value, max_color, mid_value, mid_color, min_value, min_color)
rc("colorkey  .08,.1  .10,.800  %s" % (col_string))
rc("rangecolor %s,%s %s #%s" % (attr, rep_key, col_string, model)) 
rc("rangecolor %s,a %s #%s" % (attr, col_string, model))

###########################################
# set position and initial scale for movie Z1 fixed
rc("scale 0.2")
# Record movie
i = 0
out_file = "%s/%s_%s_%s_%s_%s_viewZ1_%i" % (out_folder, rID, qID, rep, attr, metric, i)
rc("movie record")
rc("coordset #%s 1, holdSteady #%s:%s@CA load true" %(model, model, central_residue))
rc("scale 1.02 100")
rc("wait %s" % (frame_count))
rc("movie stop")
rc("movie status")
rc("movie encode %s.mp4 wait true" % (out_file))
rc("wait")
rc("scale 0.7")

##########################################
# set position and initial scale for movie Z2 fixed
rc("scale 0.2")
rc("turn y 180")
# Record movie 2
i = 1
out_file = "%s/%s_%s_%s_%s_%s_viewZ2_%i" % (out_folder, rID, qID, rep, attr, metric, i)
rc("movie record")
rc("coordset #%s 1, holdSteady #%s:%s@CA load true" %(model, model, central_residue))
rc("scale 1.02 100")
rc("wait %s" % (frame_count))
rc("movie stop")
rc("movie status")
rc("movie encode %s.mp4 wait true" % (out_file))
rc("wait")
rc("scale 0.7")

##############################################
# set position and initial scale for movie X2 fixed
rc("scale 0.2")
rc("turn y 90")
# Record movie 3
i = 2
out_file = "%s/%s_%s_%s_%s_%s_viewX2_%i" % (out_folder, rID, qID, rep, attr, metric, i)
rc("movie record")
rc("coordset #%s 1, holdSteady #%s:%s@CA load true" %(model, model, central_residue))
rc("scale 1.02 100")
rc("wait %s" % (frame_count))
rc("movie stop")
rc("movie status")
rc("movie encode %s.mp4 wait true" % (out_file))
rc("wait")
rc("scale 0.7")

##############################################
# set position and initial scale for movie X1 fixed
rc("scale 0.2")
rc("turn y 180")
# Record movie 4
i = 3
out_file = "%s/%s_%s_%s_%s_%s_viewX1_%i" % (out_folder, rID, qID, rep, attr, metric, i)
rc("movie record")
rc("coordset #%s 1, holdSteady #%s:%s@CA load true" %(model, model, central_residue))
rc("scale 1.02 100")
rc("wait %s" % (frame_count))
rc("movie stop")
rc("movie status")
rc("movie encode %s.mp4 wait true" % (out_file))
rc("wait")
rc("scale 0.7")

##############################################
# set position and initial scale for movie Y1 fixed
rc("scale 0.2")
rc("turn x 90")
# Record movie 5
i = 4
out_file = "%s/%s_%s_%s_%s_%s_viewY1_%i" % (out_folder, rID, qID, rep, attr, metric, i)
rc("movie record")
rc("coordset #%s 1, holdSteady #%s:%s@CA load true" %(model, model, central_residue))
rc("scale 1.02 100")
rc("wait %s" % (frame_count))
rc("movie stop")
rc("movie status")
rc("movie encode %s.mp4 wait true" % (out_file))
rc("wait")
rc("scale 0.7")

##############################################
# set position and initial scale for movie Y2 fixed
rc("scale 0.2")
rc("turn y 180")
# Record movie 6
i = 5
out_file = "%s/%s_%s_%s_%s_%s_viewY2_%i" % (out_folder, rID, qID, rep, attr, metric, i)
rc("movie record")
rc("coordset #%s 1, holdSteady #%s:%s@CA load true" %(model, model, central_residue))
rc("scale 1.02 100")
rc("wait %s" % (frame_count))
rc("movie stop")
rc("movie status")
rc("movie encode %s.mp4 wait true" % (out_file))
rc("wait")
rc("scale 0.7")

##############################################
# set position and initial scale for movie R1 fixed
rc("roll")
rc("scale 0.2")
rc("turn x 90")
# Record movie 7
i = 6
out_file = "%s/%s_%s_%s_%s_%s_viewR1_%i" % (out_folder, rID, qID, rep, attr, metric, i)
rc("movie record")
rc("coordset #%s 1, holdSteady #%s:%s@CA load true" %(model, model, central_residue))
rc("scale 1.02 100")
rc("wait %s" % (frame_count))
rc("movie stop")
rc("movie status")
rc("movie encode %s.mp4 wait true" % (out_file))
rc("wait")
rc("scale 0.7")

##############################################
# set position and initial scale for movie R2 fixed
rc("scale 0.2")
rc("turn x 90")
# Record movie 8
i = 7
out_file = "%s/%s_%s_%s_%s_%s_viewR2_%i" % (out_folder, rID, qID, rep, attr, metric, i)
rc("movie record")
rc("coordset #%s 1, holdSteady #%s:%s@CA load true" %(model, model, central_residue))
rc("scale 1.02 100")
rc("wait %s" % (frame_count))
rc("movie stop")
rc("movie status")
rc("movie encode %s.mp4 wait true" % (out_file))
rc("wait")
rc("scale 0.7")

###############################################
rc("stop")

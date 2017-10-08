# input file names in command line, can be hard coded or kept the same

import sys
import os

import Tkinter

from Tkinter import *

import gi
gi.require_version('Gst', '1.0')
from gi.repository import Gst, GObject

# Needed for set_window_handle():
gi.require_version('GstVideo', '1.0')
from gi.repository import GstVideo

def set_frame_handle(bus, message, frame_id):
    if not message.get_structure() is None:
        if message.get_structure().get_name() == 'prepare-window-handle':
            display_frame = message.src
            display_frame.set_property('force-aspect-ratio', True)
            display_frame.set_window_handle(frame_id)

# Only argument number checked, not validity.
number_of_file_names_given = len(sys.argv) - 1
if number_of_file_names_given < 1:
    print('Give at least one video file name.')
    sys.exit()
file_names = list()
for index in range(number_of_file_names_given):
    file_names.append(sys.argv[index + 1])

# Frame number and size
NUMBER_OF_FRAMES = len(file_names)
relative_height = 2 / float(NUMBER_OF_FRAMES)

# Initialize window
window = Tkinter.Tk()
window.title('DROIDS Output Movie')
window.geometry('450x500')

Gst.init(None)
GObject.threads_init()

# Play callback
def play():
    
    for number in range(NUMBER_OF_FRAMES):
    
        # To display in two columns
        display_frame = Tkinter.Frame(window, bg='')
        if number < NUMBER_OF_FRAMES / 2:
            relative_y = number * relative_height
            relative_x = 0.0
        elif number >= NUMBER_OF_FRAMES / 2:
            relative_y = (number - (NUMBER_OF_FRAMES / 2)) * relative_height
            relative_x = 0.5
            
        # Format each frame
        display_frame.place(relx = relative_x, rely = relative_y,
            anchor = Tkinter.NW, relwidth = 0.5, relheight = relative_height)
        frame_id = display_frame.winfo_id()
        
        
        w = Tkinter.Label(display_frame, text = file_names[number], bg = "Black", fg = "White")
        w.place(relx = 0, rely = 0)
    
        player = Gst.ElementFactory.make('playbin', None)
        fullname = os.path.abspath(file_names[number % len(file_names)])
        player.set_property('uri', 'file://%s' % fullname)
        player.set_state(Gst.State.PLAYING)

        bus = player.get_bus()
        bus.enable_sync_message_emission()
        bus.connect('sync-message::element', set_frame_handle, frame_id)
        
        
    # Add start/stop button
    B = Tkinter.Button(window, text = "Restart Video", command = play, bg = "LightGray")
    B.pack()
    B.place(anchor = SE, relx = 1, rely = 1)
    
      
play()
window.mainloop()

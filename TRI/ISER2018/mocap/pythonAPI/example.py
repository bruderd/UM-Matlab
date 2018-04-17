#!/usr/bin/python

#
# Copyright (c) PhaseSpace, Inc 2017
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# PHASESPACE, INC BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

import owl
import sys

# takes the server IP as a command line input, formatted as a string
SERVER = sys.argv[1]

# instantiate context
o = owl.Context()
# connect to server with timeout of 10000000 microseconds
o.open(SERVER, "timeout=10000000")
# initialize session
o.initialize("streaming=1")

# main loop
evt = None
while evt or (o.isOpen() and o.property("initialized")):

    # poll for events with a timeout (microseconds)
    evt = o.nextEvent(1000000)
    # nothing received, keep waiting
    if not evt: continue
    else: print evt

    # process event
    if evt.type_id == owl.Type.FRAME:
        # print markers
        if "markers" in evt:
            for m in evt.markers: print m
        # print rigids
        if "rigids" in evt:
            for r in evt.rigids: print r
    elif evt.type_id == owl.Type.ERROR:
        # handle errors
        print evt.name, evt.data
        if event.name == "fatal":
            break
    elif evt.name == "done":
        # done event is sent when master connection stops session
        print "done"
        break
# end main loop

# end session
o.done()
# close socket
o.close()

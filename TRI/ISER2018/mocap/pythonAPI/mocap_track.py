#!/usr/bin/env python

from owl import Context, Type
import rospy
from geometry_msgs.msg import Pose


def main():
    rospy.init_node('mocap_system')

    pub = rospy.Publisher('mocap', Pose, queue_size=10)

    address = "192.168.1.11"
    owl = Context()
    owl.open(address)
    owl.initialize()
    owl.frequency(50)

    tracker_id = 0
    owl.createTracker(tracker_id, "rigid", "myrigid")

    owl.assignMarker(tracker_id, 88, "88", "pos=0.0,0.0,215.9")
    owl.assignMarker(tracker_id, 93, "93", "pos=0.0,0.0,-215.9")
    owl.assignMarker(tracker_id, 90, "90", "pos=158.8,0.0,0")
    owl.assignMarker(tracker_id, 95, "95", "pos=-158.8,0.0,0")

    owl.assignMarker(tracker_id, 89, "89", "pos=158.8,0.0,215.9")
    owl.assignMarker(tracker_id, 92, "92", "pos=-158.8,0.0,-215.9")
    owl.assignMarker(tracker_id, 94, "94", "pos=-158.8,0.0,215.9")
    owl.assignMarker(tracker_id, 91, "91", "pos=158.8,0.0,-215.9")

    owl.streaming(1)

    while owl.isOpen() and owl.property("initialized"):
        event = owl.nextEvent(1000)
        if not event:
            continue
        if event.type_id == Type.ERROR:
            pass
            # print event.name, ": ", event.str
        elif event.type_id == Type.FRAME:
            if "rigids" in event:
                # print "rigids=", len(event.rigids))
                for r in event.rigids:
                    if r.cond > 0:
                        print(r.pose)
                        ret = Pose()
                        ret.position.x = r.pose[0]
                        ret.position.y = r.pose[1]
                        ret.position.z = r.pose[2]
                        ret.orientation.w = r.pose[3]
                        ret.orientation.x = r.pose[4]
                        ret.orientation.y = r.pose[5]
                        ret.orientation.z = r.pose[6]
                        # now = rospy.get_rostime()
                        # ret.header(now)
                        pub.publish(ret)

    owl.done()
    owl.close()


if __name__ == "__main__":
    main()

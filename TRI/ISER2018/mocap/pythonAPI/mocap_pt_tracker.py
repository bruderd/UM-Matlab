#!/usr/bin/env python

from owl import Context, Type
import rospy
from geometry_msgs.msg import PointStamped
from geometry_msgs.msg import Point32
import std_msgs.msg


def main():
    rospy.init_node('mocap_system')

    pub = rospy.Publisher('mocap',PointStamped,queue_size=10)

    # connect to NETGEAR48 and it'l be this IP
    address = "192.168.1.11"
    owl = Context()
    owl.open(address)
    owl.initialize()

    # how often the mocap system will send a "frame" that contains the body location
    # you can play with this, saturates the router around 1 kHz
    owl.frequency(50)

    # once you hit this point, the program will connect with the receiver and the lights should turn on
    owl.streaming(1)

    while owl.isOpen() and owl.property("initialized"):
        event = owl.nextEvent(1000)
        if not event:
            continue
        if event.type_id == Type.ERROR:
            pass
            # print event.name, ": ", event.str
        elif event.type_id == Type.FRAME:
            if "markers" in event:
                for r in event.markers:
                    if r.cond > 0:
			pt = PointStamped()

			header = std_msgs.msg.Header()
			header.stamp = rospy.Time.now()
			header.frame_id = str(r.id)

			pt.header = header

			pt.point.x = r.x
			pt.point.y = r.y
			pt.point.z = r.z

                        pub.publish(pt)

    owl.done()
    owl.close()


if __name__ == "__main__":
    main()

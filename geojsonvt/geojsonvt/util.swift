import Foundation

struct Util {

    static var activities = Dictionary<String, NSDate>()

    static func time(activity: String) {

        Util.activities[activity] = NSDate()
    }

    static func timeEnd(activity: String) {

        NSLog("\(activity): %fms", (NSDate().timeIntervalSince1970 - Util.activities[activity]!.timeIntervalSince1970) * 1000)
    }

}

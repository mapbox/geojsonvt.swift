import UIKit

class ViewController: UIViewController {

    var vt: GeoJSONVT!
    var imageView: UIImageView!
    var z: Int = 0
    var x: Int = 0
    var y: Int = 0

    override func viewDidLoad() {
        super.viewDidLoad()

        let json = NSString(contentsOfFile:
            NSBundle.mainBundle().pathForResource("threestates", ofType: "geojson")!,
            encoding: NSUTF8StringEncoding, error: nil)
        NSLog("loaded up feature JSON of %i bytes", json!.length)

        dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), { [unowned self] in
            self.vt = GeoJSONVT(data: json! as String, debug: true)
            dispatch_async(dispatch_get_main_queue(), { [unowned self] in
                self.drawTile()
            })
        });

        let size = self.view.bounds.size.width

        self.imageView = UIImageView(frame: CGRect(x: 0, y: (self.view.bounds.size.height - size) / 2, width: size, height: size))
        self.imageView.userInteractionEnabled = true
        self.imageView.addGestureRecognizer(UITapGestureRecognizer(target: self, action: "singleTap:"))
        self.imageView.addGestureRecognizer({
            let gesture = UITapGestureRecognizer(target: self, action: "twoFingerTap:")
            gesture.numberOfTouchesRequired = 2
            return gesture
            }())
        self.view.addSubview(self.imageView)
    }

    func drawTile() {
        let size = self.view.bounds.size.width

        UIGraphicsBeginImageContext(CGSize(width: size, height: size))
        let c = UIGraphicsGetCurrentContext()

        CGContextSetFillColorWithColor(c, UIColor.whiteColor().CGColor)
        CGContextFillRect(c, CGRect(x: 0, y: 0, width: size, height: size))

        CGContextSetStrokeColorWithColor(c, UIColor.redColor().CGColor)
        CGContextSetFillColorWithColor(c, UIColor.redColor().colorWithAlphaComponent(0.05).CGColor)

        let tile = self.vt.getTile(self.z, x: self.x, y: self.y)

        if (tile != nil) {
            let extent: Double = 4096
            for feature in tile!.features {
                for geometry in feature.geometry {
                    if (feature.type == .Point) {
                        let radius: CGFloat = 1
                        let point = geometry as! TilePoint
                        let x = CGFloat((Double(point.x) / extent) * Double(size))
                        let y = CGFloat((Double(point.y) / extent) * Double(size))
                        let dot = CGRect(x: (x - radius), y: (y - radius), width: (radius * 2), height: (radius * 2))
                        CGContextAddEllipseInRect(c, dot);
                    } else {
                        var pointCount = 0
                        let ring = geometry as! TileRing
                        for point in ring.points {
                            let x = CGFloat((Double(point.x) / extent) * Double(size))
                            let y = CGFloat((Double(point.y) / extent) * Double(size))
                            if (pointCount == 0) {
                                CGContextMoveToPoint(c, x, y)
                            } else {
                                CGContextAddLineToPoint(c, x, y)
                            }
                            pointCount++
                        }
                    }
                }
                if (feature.type == .Polygon) {
                    let p = CGContextCopyPath(c)
                    CGContextEOFillPath(c)
                    CGContextAddPath(c, p)
                }
                CGContextStrokePath(c)
            }

            CGContextSetStrokeColorWithColor(c, UIColor.greenColor().CGColor)
            CGContextSetLineWidth(c, 1)
            CGContextStrokeRect(c, CGRect(x: 0, y: 0, width: size, height: size))
            CGContextMoveToPoint(c, size / 2, 0)
            CGContextAddLineToPoint(c, size / 2, size)
            CGContextMoveToPoint(c, 0, size / 2)
            CGContextAddLineToPoint(c, size, size / 2)
            CGContextStrokePath(c)

            self.imageView.image = UIGraphicsGetImageFromCurrentImageContext()
            UIGraphicsEndImageContext()
        } else {
            self.zoomOut()
        }
    }

    func zoomOut() {
        self.z--
        self.x = x / 2
        self.y = y / 2

        if (z < 0) {
            z = 0
        }

        if (x < 0) {
            x = 0
        }

        if (y < 0) {
            y = 0
        }
    }

    func singleTap(gesture: UITapGestureRecognizer) {
        let left = (gesture.locationInView(gesture.view).x / self.view.bounds.size.width < 0.5)
        let top  = (gesture.locationInView(gesture.view).y / self.view.bounds.size.width < 0.5)

        self.z++
        self.x *= 2
        self.y *= 2
        if (!left) {
            self.x++
        }
        if (!top) {
            self.y++
        }

        self.drawTile()
    }

    func twoFingerTap(gesture: UITapGestureRecognizer) {
        self.zoomOut()
        self.drawTile()
    }

}

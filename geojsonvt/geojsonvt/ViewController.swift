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

        self.vt = GeoJSONVT(data: json!, debug: true)

        let size = self.view.bounds.size.width

        self.imageView = UIImageView(frame: CGRect(x: 0, y: (self.view.bounds.size.height - size) / 2, width: size, height: size))
        self.imageView.userInteractionEnabled = true
        self.imageView.addGestureRecognizer(UITapGestureRecognizer(target: self, action: "tap:"))
        self.view.addSubview(self.imageView)

        self.drawTile()
    }

    func drawTile() {
        let size = self.view.bounds.size.width

        UIGraphicsBeginImageContext(CGSize(width: size, height: size))
        let c = UIGraphicsGetCurrentContext()

        CGContextSetFillColorWithColor(c, UIColor.whiteColor().CGColor)
        CGContextFillRect(c, CGRect(x: 0, y: 0, width: size, height: size))

        let tile = self.vt.getTile(self.z, x: self.x, y: self.y)

        if (tile != nil) {
            let extent: Double = 4096
            for feature in tile!.features {
                for geometry in feature.geometry {
                    var pointCount = 0
                    let ring = geometry as TileRing
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
                    CGContextStrokePath(c)
                }
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
            self.z--
            self.x = x / 2
            self.y = y / 2
        }
    }

    func tap(gesture: UITapGestureRecognizer) {
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

}

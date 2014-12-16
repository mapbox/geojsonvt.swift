import Foundation

let json = "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[-119.53125,44.02442151965934],[-118.564453125,39.639537564366684],[-115.13671875,38.685509760012],[-113.466796875,42.74701217318067],[-115.6640625,33.94335994657882]]}},{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[-110.21484375,42.61779143282346],[-108.6328125,33.50475906922606],[-97.3828125,33.7243396617476],[-99.84374999999999,44.02442151965934],[-109.072265625,43.51668853502906]]}},{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[-94.833984375,43.96119063892024],[-94.482421875,34.161818161230386],[-88.154296875,34.74161249883172]]}},{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[-85.517578125,43.89789239125797],[-84.19921875,33.94335994657882],[-74.35546875,35.17380831799959],[-76.640625,44.77793589631623],[-84.72656249999999,44.08758502824516]]}}]}"

let extent = 4096
let padding = 8 / 512
let minPx = Int(round(-Double(padding) * Double(extent)))
let maxPx = Int(round(Double(1 + padding) * Double(extent)))

class GeoJSONVT {

    var baseZoom: Int
    var maxZoom: Int
    var maxPoints: Int
    var tolerance: Int
    var debug: Bool

    var tiles: [Int: Tile]

    var stats: [Int]

    var total = 0

    init(data: String, baseZoom: Int = 14, maxZoom: Int = 4, maxPoints: Int = 100, tolerance: Int = 3, debug: Bool = false) {
        self.baseZoom = baseZoom
        self.maxZoom = maxZoom
        self.maxPoints = maxPoints
        self.tolerance = tolerance
        self.debug = debug

        if (self.debug) {
            // time 'preprocess data'
        }

        let z2 = 1 << baseZoom

        let deserializedData = NSJSONSerialization.JSONObjectWithData(data.dataUsingEncoding(NSUTF8StringEncoding,
            allowLossyConversion: false)!, options: nil, error: nil) as Dictionary<String, AnyObject>

        let features = Convert.convert(data: deserializedData, tolerance: self.tolerance / (z2 * extent))

        if (self.debug) {
            // time end 'preprocess data'
            // time 'generate tiles up to z' + maxZoom
        }

        self.splitTile(features: features, z: 0, x: 0, y: 0)

        if (self.debug) {
            // log "features: self.tiles[0].numFeatures, points: self.tiles[0].numPoints
            // time end 'generate tiles up to z' + maxZoom
            // log "tiles generated: self. total, self.stats
        }
    }

    func splitTile(#features: [Feature], z: Int, x: Int, y: Int, cz: Int = 0, cx: Int = 0, cy: Int = 0) {

        var stack: [Any] = [features, z, x, y]

        while (stack.count > 0) {
            let features = stack.removeAtIndex(0) as [Feature]
            let z = stack.removeAtIndex(0) as Int
            let x = stack.removeAtIndex(0) as Int
            let y = stack.removeAtIndex(0) as Int

            let z2 = 1 << z
            let id = toID(z: z, x: x, y: y)
            var tile: Tile
            let tileTolerance = (z == self.baseZoom ? 0 : self.tolerance / (z2 * extent))

            if (self.tiles.indexForKey(id) != nil) {
                tile = self.tiles[id]!
            } else {
                if (self.debug) {
                    // time 'creation'
                }

                tile = Tile.createTile(features: features, z2: z2, tx: x, ty: y, tolerance: tileTolerance,
                    extent: extent, noSimplify: (z == self.baseZoom))

                self.tiles[id] = tile

                if (self.debug) {
                    // log "tile z" + z, x, y + " (features: tile.numFeatures, points: tile.numPoints, simplified: tile.numSimplified)
                    // time end 'creation'

                    if (self.stats.count - 1 >= z) {
                        self.stats[z] += 1
                    } else {
                        self.stats[z] = 1
                    }
                    self.total++
                }
            }

            if (cz == 0 && (z == self.maxZoom || tile.numPoints <= self.maxPoints ||
                self.isClippedSquare(features: tile.features)) || z == self.baseZoom || z == cz) {
                tile.source = features
                continue
            }

            if (cz > 0) {
                tile.source = features
            } else {
                tile.source = []
            }

            if (self.debug) {
                // time 'clipping'
            }

            let k1 = 0.5 * Double(padding)
            let k2 = 0.5 - k1
            let k3 = 0.5 + k1
            let k4 = 1 + k1

            var tl: [Feature]
            var bl: [Feature]
            var tr: [Feature]
            var br: [Feature]
            var left: [Feature]
            var right: [Feature]
            var m: Int
            var goLeft: Bool
            var goTop: Bool

            if (cz > 0) {
                m = 1 << (cz - z)
                goLeft = Double(cx) / Double(m) - Double(x) < 0.5
                goTop = Double(cy) / Double(m) - Double(y) < 0.5
            }

            if (cz == 0 || goLeft) {
                left = Clip.clip(features: features, scale: z2, k1: x - Int(k1), k2: x + Int(k3), axis: 0, intersect: self.intersectX)
            } else if (cz == 0 || !goLeft) {
                right = Clip.clip(features: features, scale: z2, k1: x + Int(k2), k2: x + Int(k4), axis: 0, intersect: self.intersectX)
            }

            if (left.count > 0) {
                if (cz == 0 || goTop) {
                    tl = Clip.clip(features: left, scale: z2, k1: y - Int(k1), k2: y + Int(k3), axis: 1, intersect: self.intersectY)
                } else if (cz == 0 || !goTop) {
                    bl = Clip.clip(features: left, scale: z2, k1: y + Int(k2), k2: y + Int(k4), axis: 1, intersect: self.intersectY)
                }
            }

            if (right.count > 0) {
                if (cz == 0 || goTop) {
                    tr = Clip.clip(features: right, scale: z2, k1: y - Int(k1), k2: y + Int(k3), axis: 1, intersect: self.intersectY)
                } else if (cz == 0 || !goTop) {
                    br = Clip.clip(features: right, scale: z2, k1: y + Int(k2), k2: y + Int(k4), axis: 1, intersect: self.intersectY)
                }
            }

            if (self.debug) {
                // time end 'clipping'
            }

            if (tl.count > 0) {
                stack.append(tl)
                stack.append(z + 1)
                stack.append(x * 2)
                stack.append(y * 2)
            }

            if (bl.count > 0) {
                stack.append(bl)
                stack.append(z + 1)
                stack.append(x * 2)
                stack.append(y * 2 + 1)
            }

            if (tr.count > 0) {
                stack.append(tr)
                stack.append(z + 1)
                stack.append(x * 2 + 1)
                stack.append(y * 2)
            }

            if (br.count > 0) {
                stack.append(br)
                stack.append(z + 1)
                stack.append(x * 2 + 1)
                stack.append(y * 2 + 1)
            }
        }
    }

    func getTile(z: Int, x: Int, y: Int) -> Tile {

        let id = self.toID(z: z, x: x, y: y)
        if (self.tiles.indexForKey(id) != nil) {
            return self.tiles[id]!
        }

        if (self.debug) {
            // log 'drilling down to z%d-%d-%d', z, x, y
        }

        var z0 = z
        var x0 = x
        var y0 = y
        var parent: Tile? = nil

        while (parent == nil && z0 > 0) {
            z0--
            x0 = Int(floor(Double(x0) / 2))
            y0 = Int(floor(Double(y0) / 2))
            let checkID = self.toID(z: z0, x: x0, y: y0)
            parent = (self.tiles.indexForKey(checkID) != nil ? self.tiles[checkID]! : nil)
        }

        if (self.debug) {
            // log ''found parent tile z%d-%d-%d', z0, x0, y0
        }

        if (parent!.source.count > 0) {
            if (self.isClippedSquare(features: parent!.features)) {
                return parent!
            }

            if (self.debug) {
                // time 'drilling down'
            }

            self.splitTile(features: parent!.source, z: z0, x: x0, y: y0, cz: z, cx: x, cy: y)

            if (self.debug) {
                // time end 'drilling down
            }

            return self.tiles[id]!
        }
    }

    func isClippedSquare(#features: [Feature]) -> Bool {

        if (features.count != 1) {
            return false
        }

        let feature = features.first!

        if (feature.type != .PolygonType || feature.geometry.count > 1) {
            return false
        }

        for i in 0...(feature.geometry.first! as Line).points.count {
            let p = (feature.geometry.first! as Line).points[i] as Point
            if ((p.x != minPx && p.x != maxPx) ||
                (p.y != minPx && p.y != maxPx)) {
                    return false
            }
        }

        return true
    }

    func toID(#z: Int, x: Int, y: Int) -> Int {
        return (((1 << z) * y + x) * 32) + z;
    }

    func intersectX(a: Point, b: Point, x: Int) -> [Int] {
        let r1 = x
        let r2: Int = {
            let e1 = x - a.x
            let e2 = b.y - a.y
            let e3 = b.x - a.x
            let e4 = a.y
            return e1 * e2 / e3 + e4
            }()
        let r3 = 1

        return [r1, r2, r3];
    }

    func intersectY(a: Point, b: Point, y: Int) -> [Int] {
        let r1: Int = {
            let e1 = y - a.y
            let e2 = b.x - a.x
            let e3 = b.y - a.y
            let e4 = a.x
            return e1 * e2 / e3 + e4
            }()
        let r2 = y
        let r3 = 1

        return [r1, r2, r3]
    }

    func extend(var dest: [Point], src: [Point]) -> [Point] {
        for i in 0...src.count {
            dest[i] = src[i]
        }
        return dest
    }

}

class Geometry {

}

class Point: Geometry {

    var x: Int = -1
    var y: Int = -1

    func isValid() -> Bool {
        return x >= 0 && y >= 0
    }

}

class Line: Geometry {

    var points = [Point]()
    var dist = 0
    var area = 0

}

enum FeatureType: Int {
    case PointType = 1
    case LineStringType
    case PolygonType
}

struct Feature {

    var geometry: [Geometry]
    var type: FeatureType
    var tags: [Int]
    var minPoint: Point
    var maxPoint: Point

}

class Tile {

    var features: [Feature]
    var numPoints: Int
    var numSimplified: Int
    var numFeatures: Int
    var source: [Feature]

    init() {

    }

    class func createTile(#features: [Feature], z2: Int, tx: Int, ty: Int, tolerance: Int, extent: Int, noSimplify: Bool) -> Tile {

        var tile = Tile()

        for i in 0...features.count {
            tile.numFeatures++
            Tile.addFeature(tile: tile, feature: features[i], z2: z2, tx: tx, ty: ty, tolerance: tolerance,
                extent: extent, noSimplify: noSimplify)
        }

        return tile
    }

    private class func addFeature(#tile: Tile, feature: Feature, z2: Int, tx: Int, ty: Int, tolerance: Int, extent: Int, noSimplify: Bool) {

        let geom = feature.geometry
        let type = feature.type
        var transformed: [Geometry]
        let sqTolerance = tolerance * tolerance
        var ring: Line

        if (type == .PointType) {
            for i in 0...geom.count {
                transformed.append(Tile.transformPoint(p: (geom[i] as Point), z2: z2, tx: tx, ty: ty, extent: extent))
                tile.numPoints++
                tile.numSimplified++
            }
        } else {
            for i in 0...geom.count {
                ring = geom[i] as Line

                if (!noSimplify && ((type == .LineStringType && ring.dist < tolerance) ||
                    (type == .PolygonType && ring.area < sqTolerance))) {
                        tile.numPoints += ring.points.count
                        continue
                }

                var transformedRing: Line

                for j in 0...ring.points.count {
                    let p = ring.points[j]
                    if (noSimplify || p[2] > sqTolerance) {
                        transformedRing.append(Tile.transformPoint(p: p, z2: z2, tx: tx, ty: ty, extent: extent))
                        tile.numSimplified++
                    }
                    tile.numPoints++
                }

                transformed.append(transformedRing)
            }
        }

        if (transformed.count > 0) {
            tile.features.append(Feature(geometry: transformed, type: type, tags: feature.tags, min: Point(), max: Point()))
        }
    }

    private class func transformPoint(#p: Point, z2: Int, tx: Int, ty: Int, extent: Int) -> Point {

        let x = round(Double(extent * (p.x * z2 - tx)))
        let y = round(Double(extent * (p.y * z2 - ty)))


        var t = Point()
        t.x = Int(x)
        t.y = Int(y)

        return t
    }

}

class Clip {

    class func clip(#features: [Feature], scale: Int, var k1: Int, var k2: Int, axis: Int, intersect: (Point, Point, Int) -> [Int]) -> [Feature] {

        k1 /= scale
        k2 /= scale

        var clipped = [Feature]()

        for i in 0...features.count {

            let feature = features[i]
            let geometry = feature.geometry
            let type = feature.type
            var min: Int
            var max: Int

            if (feature.min.isValid()) {
                min = (axis == 0 ? feature.min.x : feature.min.y)
                max = (axis == 0 ? feature.max.x : feature.max.y)

                if (min >= k1 && max <= k2) {
                    clipped.append(feature)
                    continue
                } else if (min > k2 || max < k1) {
                    continue
                }

                var slices = (type == FeatureType.PointType ?
                    Clip.clipPoints(geometry: geometry, k1: k1, k2: k2, axis: axis) :
                    Clip.clipGeometry(geometry: geometry, k1: k1, k2: k2, axis: axis, intersect: intersect, closed: (type == FeatureType.PolygonType)))

                if (slices.count > 0) {
                    clipped.append(Feature(geometry: slices, type: type, tags: features[i].tags, min: Point(), max: Point()))
                }
            }

            return (clipped.count > 0 ? clipped : [])
        }
    }

    class func clipPoints(geometry: [Point], k1: Double, k2: Double, axis: Int) -> [Point] {

        var slice = [Point]()

        for i in 0...geometry.count {
            let a = geometry[i]
            let ak = (axis == 0 ? a.x : a.y)

            if (ak >= Int(k1) && ak <= Int(k2)) {
                slice.append(a)
            }
        }

        return slice
    }

    class func clipGeometry(geometry: [Line], k1: Double, k2: Double, axis: Int, intersect: ([Int], [Int], Int -> [Int]), closed: Bool) -> [Slice] {

        var slices = [Slice]()

        for i in 0...geometry.count {

            var ak = 0
            var bk = 0
            var b: Point
            let points = geometry[i]
            let area = points.area
            let dist = points.dist
            let len = points.points.count
            var a: Point

            var slice = Slice()

            for j in 0...(len - 1) {
                a = (b.isValid() ? b : points.points[j])
                b = points.points[j + 1]
                ak = (bk > 0 ? bk : (axis == 0 ? a.x : a.y))
                bk = (axis == 0 ? b.x : b.y)

                if (ak < Int(k1)) {
                    if (bk > Int(k2)) {
                        slice.append(intersect(a, b, Int(k1)), intersect(a, b, Int(k2)))
                        if (!closed) {
                            slice = newSlice(slices: slices, slice: slice, area: area, dist: dist)
                        }
                    } else if (bk <= Int(k2)) {
                        slice.append(intersect(a, b, k2))
                    }
                } else {
                    slice.append(a)

                    if (bk < Int(k1)) {
                        slice.append(intersect(a, b, k1))
                        if (!closed) {
                            slice = newSlice(slices: slices, slice: slice, area: area, dist: dist)
                        }
                    } else if (bk > Int(k2)) {
                        slice.append(intersect(a, b, k2))
                        if (!closed) {
                            slice = newSlice(slices: slices, slice: slice, area: area, dist: dist)
                        }
                    }
                }
            }

            a = points.points[len - 1]
            ak = (axis == 0 ? a.x : a.y)

            let sliceLen = slice.count

            if (closed && slice[0] != slice[sliceLen - 1]) {
                slice.append(slice[0])
            }

            newSlice(slices: slices, slice: slice, area: area, dist: dist)

            return slices
        }
    }

    class func newSlice(var slices: [Slice], var slice: Slice, area: Int, dist: Int) -> Point {

        if (slice.points.count > 0) {
            slice.area = area
            slice.dist = dist
            slices.append(slice)
        }

        return Point()
    }
}


struct Slice {

    var points: [Point]
    var area: Int
    var dist: Int

}



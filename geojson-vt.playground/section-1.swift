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

    var stats: [AnyObject]

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

        let deserializedData = NSJSONSerialization.JSONObjectWithData(data.dataUsingEncoding(NSUTF8StringEncoding, allowLossyConversion: false)!,
            options: nil, error: nil) as NSDictionary

        let dictionary = [String: AnyObject]()

        for key: String in deserializedData.allKeys {
            dictionary[key] = deserializedData.objectForKey(key)
        }

        let features = Convert.convert(data: dictionary, tolerance: self.tolerance / (z2 * extent))

        if (self.debug) {
            // time end 'preprocess data'
            // time 'generate tiles up to z' + maxZoom
        }

        self.splitTile(features, 0, 0, 0)

        if (self.debug) {
            // log "features: self.tiles[0].numFeatures, points: self.tiles[0].numPoints
            // time end 'generate tiles up to z' + maxZoom
            // log "tiles generated: self. total, self.stats
        }
    }

    func splitTile(features: [Feature], z: Int, x: Int, y: Int, cz: Int = 0, cx: Int = 0, cy: Int = 0) {

        var stack: [AnyObject] = [features, z, x, y]

        while (stack.count > 0) {
            let features: Feature = stack.removeAtIndex(0)
            let z: Int = stack.removeAtIndex(0)
            let x: Int = stack.removeAtIndex(0)
            let y: Int = stack.removeAtIndex(0)

            let z2 = 1 << z
            let id = toID(z: z, x: x, y: y)
            var tile: Tile
            let tileTolerance = (z == self.options[OptionType.baseZoom] ? 0 : self.tolerance / (z2 * extent))

            if (self.tiles.indexForKey(id) != nil) {
                tile = self.tiles[id]!
            } else {
                if (self.debug) {
                    // time 'creation'
                }

                tile = Tile.createTile(features, z2: z2, tx: x, ty: y, tolerance: tileTolerance, extent: extent, noSimplify: (z == self.baseZoom))
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
                    isClippedSquare(tile.features)) || z == self.baseZoom || z == cz) {
                tile.source = features
                continue
            }

            if (cz > 0) {
                tile.source = features
            } else {
                tile.source = nil
            }

            if (self.debug) {
                // time 'clipping'
            }

            let k1 = 0.5 * Double(padding)
            let k2 = 0.5 - k1
            let k3 = 0.5 + k1
            let k4 = 1 + k1

            var tl: Double
            var bl: Double
            var tr: Double
            var br: Double
            var left: Double
            var right: Double
            var m: Int
            var goLeft: Bool
            var goTop: Bool

            if (cz > 0) {
                m = 1 << (cz - z)
                goLeft = Double(cx) / Double(m) - Double(x) < 0.5
                goTop = Double(cy) / Double(m) - Double(y) < 0.5
            }

            tl = bl = tr = br = left = right = -1

            if (cz == 0 || goLeft) {
                left = Clip.clip(features, z2, x - k1, x + k3, 0, intersectX)
            } else if (cz == 0 || !goLeft) {
                right = Clip.clip(features, z2, x + k2, x + k4, 0, intersectX)
            }

            if (left >= 0) {
                if (cz == 0 || goTop) {
                    tl = Clip.clip(left, z2, y - k1, y + k3, 1, intersectY)
                } else if (cz == 0 || !goTop) {
                    bl = Clip.clip(left, z2, y + k2, y + k4, 1, intersectY)
                }
            }

            if (right >= 0) {
                if (cz == 0 || goTop) {
                    tr = clip(right, z2, y - k1, y + k3, 1, intersectY)
                } else if (cz == 0 || !goTop) {
                    br = clip(right, z2, y + k2, y + k4, 1, intersectY)
                }
            }

            if (self.debug()) {
                // time end 'clipping'
            }

            if (tl >= 0) {
                stack.append(tl)
                stack.append(z + 1)
                stack.append(x * 2)
                stack.append(y * 2)
            }

            if (bl >= 0) {
                stack.append(bl)
                stack.append(z + 1)
                stack.append(x * 2)
                stack.append(y * 2 + 1)
            }

            if (tr >= 0) {
                stack.append(tr)
                stack.append(z + 1)
                stack.append(x * 2 + 1)
                stack.append(y * 2)
            }

            if (br >= 0) {
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
            x0 = floor(Double(x0) / 2)
            y0 = floor(Double(y0) / 2)
            let checkID = self.toID(z: z0, x: x0, y: y0)
            parent = (self.tiles.indexForKey(checkID) ? self.tiles[checkID]! : nil)
        }

        if (self.debug) {
            // log ''found parent tile z%d-%d-%d', z0, x0, y0
        }

        if (parent!.source != nil) {
            if (self.isClippedSquare(features: parent!.features)) {
                return parent
            }

            if (self.debug) {
                // time 'drilling down'
            }

            self.splitTile(parent!.source, z: z0, x: x0, y: y0, cz: z, cx: x, cy: y)

            if (self.debug) {
                // time end 'drilling down
            }

            return self.tiles[id]
        }
    }

    func isClippedSquare(features: [Feature]) -> Bool {

        if (features.count != 1) {
            return false
        }

        let feature = features.first!

        if (feature.type != .Polygon || feature.geometry.count > 1) {
            return false
        }

        for i in 0...feature.geometry.first!.count {
            let p = feature.geometry.first![i]
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

    func intersectX(a: [Int], b: [Int], x: Int) -> [Int] {
        return [x, (x - a[0]) * (b[1] - a[1]) / (b[0] - a[0]) + a[1], 1];
    }

    func intersectY(a: [Int], b: [Int], x: Int) -> [Int] {
        return [(y - a[1]) * (b[0] - a[0]) / (b[1] - a[1]) + a[0], y, 1];
    }

    func extend(var dest: [Point], src: [Point]) -> [Point] {
        for i in 0...src.count {
            dest[i] = src[i]
        }
        return dest
    }

}

struct Point {

    var x: Int
    var y: Int

}

struct Feature {

    enum Type: Int {
        case Point = 1
        case LineString
        case Polygon
    }

    var geometry: [AnyObject]
    var type: Type
    var tags: [Int]
    var min: Point
    var max: Point

}

class Tile {

    var features: [AnyObject]
    var numPoints: Int
    var numSimplified: Int
    var numFeatures: Int
    var source: AnyObject?

    init() {

    }

    class func createTile(features: [AnyObject], z2: Int, tx: Int, ty: Int, tolerance: Int, extent: Int, noSimplify: Bool) -> Tile {

        var tile = Tile()

        for i in 0...features.count {
            tile.numFeatures++
            addFeature(tile, features[i], z2, tx, ty, tolerance, extent, noSimplify)
        }

        return tile
    }

    private class func addFeature(tile: Tile, feature: Feature, z2: Int, tx: Int, ty: Int, tolerance: Int, extent: Int, noSimplify: Bool) {

        let geom = feature.geometry
        let type = feature.type
        var transformed: [Point]
        let sqTolerance = tolerance * tolerance
        var ring: [Point]

        if (type == .Point) {
            for i in 0...geom.count {
                transformed.append(Tile.transformPoint(geom[i], z2: z2, tx: tx, ty: ty, extent: extent))
                tile.numPoints++
                tile.numSimplified++
            }
        } else {
            for i in 0...geom.count {
                ring = geom[i]

                if (!noSimplify && ((type == .LineString && ring.dist < tolerance) ||
                    (type == .Polygon && ring.area < sqTolerance))) {
                        tile.numPoints += ring.count
                        continue
                }

                var transformedRing: [Point]

                for j in 0...ring.count {
                    let p = ring[j]
                    if (noSimplify || p[2] > sqTolerance) {
                        transformedRing.append(Tile.transformPoint(p, z2: z2, tx: tx, ty: ty, extent: extent))
                        tile.numSimplified++
                    }
                    tile.numPoints++
                }

                transformed.append(transformedRing)
            }
        }

        if (transformed.count > 0) {
            tile.features.append(Feature(geometry: transformed, type: type, tags: feature.tags))
        }
    }

    private class func transformPoint(p: Point, z2: Int, tx: Int, ty: Int, extent: Int) -> Point {

        let x = round(Double(extent * (p.x * z2 - tx)))
        let y = round(Double(extent * (p.y * z2 - ty)))

        return Point(x: Int(x), y: Int(y))
    }

}

class Clip {

    class func clip(features: [Feature], scale: Int, var k1: Int, var k2: Int, axis: Int, intersect: [Int]) -> [Feature] {

        k1 /= scale
        k2 /= scale

        var clipped = [Feature]()

        for i in 0...features.count {

            let feature = features[i]
            let geometry = feature.geometry
            let type = feature.type
            var min: Int
            var max: Int

            if (feature.min) {
                min = (axis == 0 ? feature.min.x : feature.min.y)
                max = (axis == 0 ? feature.max.x : feature.max.y)

                if (min >= k1 && max <= k2) {
                    clipped.append(feature)
                    continue
                } else if (min > k2 || max < k1) {
                    continue
                }

                var slices = (type == .Point ?
                    clipPoints(geometry, k1: k1, k2: k2, axis: axis) :
                    clipGeometry(geometry, k1: k1, k2: k2, axis: axis, intersect: intersect, closed: (type == .Polygon)))

                if (slices.count > 0) {
                    clipped.append(Feature(geometry: slices, type: type, tags: features.[i].tags, min: nil, max: nil))
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

            if (ak >= k1 && ak <= k2) {
                slice.append(a)
            }
        }

        return slice
    }

    class func clipGeometry(geometry: [Point], k1: Double, k2: Double, axis: Int, intersect: [Int], closed: Bool) -> [Slice] {

        var slices = [Slice]()

        for i in 0...geometry.count {

            var ak = 0
            var bk = 0
            var b: Point
            let points = geometry[i]
            let area = points.area
            let dist = points.dist
            let len = points.count
            var a: Point

            var slice = [Point]

            for j in 0...(len - 1) {
                a = (b ? b : points[j])
                b = points[j + 1]
                ak = (bk ? bk : (axis == 0 ? a.x : a.y))
                bk = (axis == 0 ? b.x b.y)

                if (ak < k1) {
                    if (bk > k2) {
                        slice.append(intersect(a, b, k1), intersect(a, b, k2))
                        if (!closed) {
                            slice = newSlice(slices: slices, slice: slice, area: area, dist: dist)
                        }
                    } else if (bk <= k2) {
                        slice.append(intersect(a, b, k2))
                    }
                } else {
                    slice.append(a)

                    if (bk < k1) {
                        slice.append(intersect(a, b, k1))
                        if (!closed) {
                            slice = newSlice(slices: slices, slice: slice, area: area, dist: dist)
                        }
                    } else if (bk > k2) {
                        slice.append(intersect(a, b, k2))
                        if (!closed) {
                            slice = newSlice(slices: slices, slice: slice, area: area, dist: dist)
                        }
                    }
                }
            }

            a = points[len - 1]
            ak = (axis == 0 ? a.x : a.y)

            let sliceLen = slice.count

            if (closed && slice[0] != slice[sliceLen - 1]) {
                slice.append(slice[0])
            }

            newSlice(slices: slices, slice: slice, area: area, dist: dist)

            return slices
        }
    }

    class func newSlice(slices: [Slice], slice: Slice, area: Int, dist: Int) -> [Int] {

        if (slice.points.count > 0) {
            slice.area = area
            slice.dist = dist
            slices.append(slice)
        }

        return [Int]()
    }
}

class Simplify {

    class func simplify(points: [Point], tolerance: Int) {

        let sqTolerance = tolerance * tolerance
        let len = points.count
        var first = 0
        var last = len - 1
        var stack = [Int]()
        var maxSqDist: Int
        var sqDist: Int
        var index: Int

        points[first][2] = 1
        points[last][2] = 1

        while (last > 0) {

            maxSqDist = 0

            for i in (first + 1)...last {
                sqDist = getSqSegDist(points[i], a: points[first], b: points[last])

                if (sqDist > maxSqDist) {
                    index = i
                    maxSqDist = sqDist
                }
            }

            if (maxSqDist > sqTolerance) {
                points[index][2] = maxSqDist
                stack.append(first)
                stack.append(index)
                stack.append(index)
                stack.append(last)
            }

            last = stack.removeLast()
            first = stack.removeLast()
        }
    }

    class func getSqSegDist(p: Point, a: Point, b: Point) -> Int {

        var x = a.x
        var y = a.y
        var dx = b.x - a.x
        var dy = b.y - a.y

        if (dx != 0 || dy != 0) {

            let t = ((p.x - a.x) * dx + (p.y - a.y) * dy) / (dx * dx + dy * dy)

            if (t > 1) {
                x = b.x
                y = b.y
            } else if (t > 0) {
                x += dx * t
                y += dy * t
            }
        }

        dx = p.x - x
        dy = p.y - y

        return dx * dx + dy * dy
    }
}

class Slice {

    var points: [Point]
    var area: Int
    var dist: Int

}

class LonLat {

    var lon: Double
    var lat: Double

}

class Convert {

    class func convert(data: [String: AnyObject], tolerance: Int) -> [Feature] {

        var features = [Feature]()

        if ((data["type"] as String) == "FeatureCollection") {
            for i in 0...(data["features"] as [Feature]).count {
                convertFeature(features: features, data: data["features"][i], tolerance: tolerance)
            }
        } else if (data["type"] == "Feature") {
            convertFeature(features: features, data: data, tolerance: tolerance)
        } else {
            convertFeature(features: features, [geometry: data], tolerance)
        }

        return features
    }

    class func convertFeature(features: [Feature], feature: Feature, tolerance: Int) {

        let geom = feature.geometry
        let type = geom.type
        let coords = geom.coordinates
        let tags = feature.properties
        var rings = [Feature]()

        if (type == .Point) {
            features.append(create(tags: tags, type: .Point, geometry: [projectPoint(coords)])
        } else if (type == .MultiPoint) {
            features.append(create(tags: tags, type: .Point, project(lonlats: coords)))
        } else if (type == .LineString) {
            features.append(create(tags: tags, type: .LineString, [project(lonlats: coords, tolerance: tolerance)]))
        } else if (type == .MultiLineString || type == .Polygon) {
            for i in 0...coords.count {
                rings.append(project(lonlats: coords[i], tolerance: tolerance)
            }
            features.append(create(tags: tags, (type == .Polygon ? .Polygon : .LineString), geometry: rings))
        } else if (type == .MultiPolygon) {
            for i in 0...coords.count {
                for j in 0...coords[i].count {
                    rings.append(project(lonlats: coords[i][j], tolerance))
                }
            }
            features.append(create(tags: tags, type: .Polygon, features: rings))
        } else if (type == .GeometryCollection) {
            for i in 0...geom.geometries.count {
                convertFeature(features: features, feature: Feature(geometry: geom.geometries[i], properties: tags) tolerance: tolerance)
            }
        } else {
            prinln("Unsupported GeoJSON type: " + geom.type)
        }
    }

    class func create(tags: [Int], type: Feature.Type, geometry: [Feature]) -> Feature {

        var feature = Feature(geometry: geometry, type: type, tags: tags, min: Point(x: 1, y: 1), max: Point(x: 0, y: 0))
        calcBBox(feature)

        return feature
    }

    class func project(lonlats: [LonLat], tolerance: Int = 0) {

        var projected = [Point]()
        for i in 0...lonlats.count {
            projected.append(projectPoint(lonlats[i]))
        }
        if (tolerance > 0) {
            Simplify.simplify(projected, tolerance: tolerance)
            calcSize(projected)
        }

        return projected
    }

    class func projectPoint(p: Point) -> [Int] {

        let sin = sin(Double(p.y) * M_PI / 180)
        let x = Double(p.x) / 360 + 0.5
        let y = 0.5 - 0.25 * log((1 + sin) / (1 - sin)) / M_PI

        return [x, y, 0]
    }

    class func calcSize(points: [Point]) {

        var area = 0
        var dist = 0
        var a: Point
        var b: Point

        for i in 0...points.count {
            a = (b ? b : points[i])
            b = point[i + 1]

            area += a.x * b.y - b.x * a.y

            dist += abs(b.x - a.x) + abs(b.y - a.y)
        }

        points.area = abs(area / 2)
        points.dist = dist
    }

    class func calcBBox(var feature: Feature) -> Feature {

        let geometry = feature.geometry
        let min = feature.min
        let max = feature.max

        if (feature.type == .Point) {
            calcRingBBox(min, max, geometry)
        } else {
            for i in 0...geometry.count {
                calcRingBBox(min, max, geometry[i])
            }
        }

        return feature
    }

    class func calcRingBBox(var min: Point, var max: Point, points: [Point]) {

        for i in 0...points.count {
            let p = points[i]
            min.x = min(p.x, min.x)
            max.x = max(p.x, max.x)
            min.y = min(p.y, min.y)
            max.y = max(p.y, max.y)
        }
    }

}

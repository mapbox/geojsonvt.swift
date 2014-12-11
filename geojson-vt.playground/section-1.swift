import Foundation

let json = "{\"type\":\"FeatureCollection\",\"features\":[{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[-119.53125,44.02442151965934],[-118.564453125,39.639537564366684],[-115.13671875,38.685509760012],[-113.466796875,42.74701217318067],[-115.6640625,33.94335994657882]]}},{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[-110.21484375,42.61779143282346],[-108.6328125,33.50475906922606],[-97.3828125,33.7243396617476],[-99.84374999999999,44.02442151965934],[-109.072265625,43.51668853502906]]}},{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[-94.833984375,43.96119063892024],[-94.482421875,34.161818161230386],[-88.154296875,34.74161249883172]]}},{\"type\":\"Feature\",\"properties\":{},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[[-85.517578125,43.89789239125797],[-84.19921875,33.94335994657882],[-74.35546875,35.17380831799959],[-76.640625,44.77793589631623],[-84.72656249999999,44.08758502824516]]}}]}"

let extent = 4096
let padding = 8 / 512
let minPx = Int(round(-Double(padding) * Double(extent)))
let maxPx = Int(round(Double(1 + padding) * Double(extent)))

class GeoJSONVT {

    enum OptionType {
        case maxZoom
        case baseZoom
        case maxPoints
        case tolerance
        case debug
    }

    let options = [
        OptionType.maxZoom: 4,
        OptionType.baseZoom: 14,
        OptionType.maxPoints: 100,
        OptionType.tolerance: 3,
        OptionType.debug: 0
    ]

    var tiles: [Int: Tile]

    var stats: [AnyObject]

    var total = 0

    func debug() -> Bool {
        return (self.options[OptionType.debug] > 0)
    }

    init(data: String, options: (name: OptionType, value: Int)...) {
        for option in options {
            self.options[option.name] = option.value
        }

        if (self.debug()) {
            // time 'preprocess data'
        }

        let z2 = 1 << self.options[OptionType.baseZoom]!

        let features = Convert(data, self.options[OptionType.tolerance] / (z2 * extent))

        if (self.debug()) {
            // time end 'preprocess data'
            // time 'generate tiles up to z' + maxZoom
        }

        splitTile(features, 0, 0, 0)

        if (self.debug()) {
            // log "features: self.tiles[0].numFeatures, points: self.tiles[0].numPoints
            // time end 'generate tiles up to z' + maxZoom
            // log "tiles generated: self. total, self.stats
        }
    }

    func splitTile(features: [Feature], z: Int, x: Int, y: Int, cz: Int, cx: Int, cy: Int) {

        var stack: [AnyObject] = [features, z, x, y]

        while (stack.count > 0) {
            let features: Feature = stack.removeAtIndex(0)
            let z: Int = stack.removeAtIndex(0)
            let x: Int = stack.removeAtIndex(0)
            let y: Int = stack.removeAtIndex(0)

            let z2 = 1 << z
            let id = toID(z, x, y)
            var tile: Tile
            let tileTolerance = (z == self.options[OptionType.baseZoom] ? 0 : self.options[OptionType.tolerance] / (z2 * extent))

            if (self.tiles.indexForKey(id) != nil) {
                tile = self.tiles[id]!
            } else {
                if (self.debug()) {
                    // time 'creation'
                }

                tile = Tile.createTile(features, z2: z2, tx: x, ty: y, tolerance: tileTolerance, extent: extent, noSimplify: (z == self.options[OptionType.baseZoom]))
                self.tiles[id] = tile

                if (self.debug()) {
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

            if (!cz && (z == self.options[OptionType.maxZoom] || tile.numPoints <= self.options[OptionType.maxPoints] ||
                    isClippedSquare(tile.features)) || z == self.options[OptionType.baseZoom] || z == cz) {
                tile.source = features
                continue
            }

            if (cz) {
                tile.source = features
            } else {
                tile.source = nil
            }

            if (self.debug()) {
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
                left = clip(features, z2, x - k1, x + k3, 0, intersectX)
            } else if (cz == 0 || !goLeft) {
                right = clip(features, z2, x + k2, x + k4, 0, intersectX)
            }

            if (left >= 0) {
                if (cz == 0 || goTop) {
                    tl = clip(left, z2, y - k1, y + k3, 1, intersectY)
                } else if (cz == 0 || !goTop) {
                    bl = clip(left, z2, y + k2, y + k4, 1, intersectY)
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

//            var k1 = 0.5 * padding,
//            k2 = 0.5 - k1,
//            k3 = 0.5 + k1,
//            k4 = 1 + k1,
//
//            tl, bl, tr, br, left, right,
//            m, goLeft, goTop;
//
//            if (cz) { // if we have a specific tile to drill down to, calculate where to go
//                m = 1 << (cz - z);
//                goLeft = cx / m - x < 0.5;
//                goTop = cy / m - y < 0.5;
//            }
//
//            tl = bl = tr = br = left = right = null;
//
//            if (!cz ||  goLeft) left  = clip(features, z2, x - k1, x + k3, 0, intersectX);
//            if (!cz || !goLeft) right = clip(features, z2, x + k2, x + k4, 0, intersectX);
//
//            if (left) {
//                if (!cz ||  goTop) tl = clip(left, z2, y - k1, y + k3, 1, intersectY);
//                if (!cz || !goTop) bl = clip(left, z2, y + k2, y + k4, 1, intersectY);
//            }
//
//            if (right) {
//                if (!cz ||  goTop) tr = clip(right, z2, y - k1, y + k3, 1, intersectY);
//                if (!cz || !goTop) br = clip(right, z2, y + k2, y + k4, 1, intersectY);
//            }
//
//            if (debug > 1) console.timeEnd('clipping');
//
//            if (tl) stack.push(tl, z + 1, x * 2,     y * 2);
//            if (bl) stack.push(bl, z + 1, x * 2,     y * 2 + 1);
//            if (tr) stack.push(tr, z + 1, x * 2 + 1, y * 2);
//            if (br) stack.push(br, z + 1, x * 2 + 1, y * 2 + 1);

        }
    }

    func getTile(z: Int, x: Int, y: Int) -> Tile {

//        var id = toID(z, x, y);
//        if (this.tiles[id]) return this.tiles[id];
//
//        var debug = this.options.debug;
//
//        if (debug > 1) console.log('drilling down to z%d-%d-%d', z, x, y);
//
//        var z0 = z,
//        x0 = x,
//        y0 = y,
//        parent;
//
//        while (!parent && z0 > 0) {
//            z0--;
//            x0 = Math.floor(x0 / 2);
//            y0 = Math.floor(y0 / 2);
//            parent = this.tiles[toID(z0, x0, y0)];
//        }
//
//        if (debug > 1) console.log('found parent tile z%d-%d-%d', z0, x0, y0);
//
//        if (parent.source) {
//            if (isClippedSquare(parent.features)) return parent;
//
//            if (debug) console.time('drilling down');
//            this.splitTile(parent.source, z0, x0, y0, z, x, y);
//            if (debug) console.timeEnd('drilling down');
//        }
//        
//        return this.tiles[id];

    }

    func isClippedSquare(features: [Feature]) -> Bool {

//        if (features.length !== 1) return false;
//
//        var feature = features[0];
//        if (feature.type !== 3 || feature.geometry.length > 1) return false;
//
//        for (var i = 0; i < feature.geometry[0].length; i++) {
//            var p = feature.geometry[0][i];
//            if ((p[0] !== minPx && p[0] !== maxPx) ||
//            (p[1] !== minPx && p[1] !== maxPx)) return false;
//        }
//        return true;

    }

    func toID(z: Int, x: Int, y: Int) -> Int {
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
    var source: AnyObject

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

//        var geom = feature.geometry,
//        type = feature.type,
//        transformed = [],
//        sqTolerance = tolerance * tolerance,
//        i, j, ring, p;
//
//        if (type === 1) {
//            for (i = 0; i < geom.length; i++) {
//                transformed.push(transformPoint(geom[i], z2, tx, ty, extent));
//                tile.numPoints++;
//                tile.numSimplified++;
//            }
//
//        } else {
//
//            // simplify and transform projected coordinates for tile geometry
//            for (i = 0; i < geom.length; i++) {
//                ring = geom[i];
//
//                // filter out tiny polylines & polygons
//                if (!noSimplify && ((type === 2 && ring.dist < tolerance) ||
//                    (type === 3 && ring.area < sqTolerance))) {
//                        tile.numPoints += ring.length;
//                        continue;
//                }
//
//                var transformedRing = [];
//
//                for (j = 0; j < ring.length; j++) {
//                    p = ring[j];
//                    // keep points with significance > tolerance and points introduced by clipping
//                    if (noSimplify || p[2] > sqTolerance) {
//                        transformedRing.push(transformPoint(p, z2, tx, ty, extent));
//                        tile.numSimplified++;
//                    }
//                    tile.numPoints++;
//                }
//
//                transformed.push(transformedRing);
//            }
//        }
//        
//        if (transformed.length) {
//            tile.features.push({
//                geometry: transformed,
//                type: type,
//                tags: feature.tags || null
//            });
//        }

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

//        k1 /= scale;
//        k2 /= scale;
//
//        var clipped = [];
//
//        for (var i = 0; i < features.length; i++) {
//
//            var feature = features[i],
//            geometry = feature.geometry,
//            type = feature.type,
//            min, max;
//
//            if (feature.min) {
//                min = feature.min[axis];
//                max = feature.max[axis];
//
//                if (min >= k1 && max <= k2) { // trivial accept
//                    clipped.push(feature);
//                    continue;
//                } else if (min > k2 || max < k1) continue; // trivial reject
//            }
//
//            var slices = type === 1 ?
//                clipPoints(geometry, k1, k2, axis) :
//                clipGeometry(geometry, k1, k2, axis, intersect, type === 3);
//
//            if (slices.length) {
//                clipped.push({
//                    geometry: slices,
//                    type: type,
//                    tags: features[i].tags || null
//                });
//            }
//        }
//        
//        return clipped.length ? clipped : null;

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

//        var slices = [];
//
//        for (var i = 0; i < geometry.length; i++) {
//
//            var ak = 0,
//            bk = 0,
//            b = null,
//            points = geometry[i],
//            area = points.area,
//            dist = points.dist,
//            len = points.length,
//            a, j;
//
//            var slice = [];
//
//            for (j = 0; j < len - 1; j++) {
//                a = b || points[j];
//                b = points[j + 1];
//                ak = bk || a[axis];
//                bk = b[axis];
//
//                if (ak < k1) {
//
//                    if ((bk > k2)) { // ---|-----|-->
//                        slice.push(intersect(a, b, k1), intersect(a, b, k2));
//                        if (!closed) slice = newSlice(slices, slice, area, dist);
//
//                    } else if (bk >= k1) slice.push(intersect(a, b, k1)); // ---|-->  |
//
//                } else if (ak > k2) {
//
//                    if ((bk < k1)) { // <--|-----|---
//                        slice.push(intersect(a, b, k2), intersect(a, b, k1));
//                        if (!closed) slice = newSlice(slices, slice, area, dist);
//
//                    } else if (bk <= k2) slice.push(intersect(a, b, k2)); // |  <--|---
//
//                } else {
//
//                    slice.push(a);
//
//                    if (bk < k1) { // <--|---  |
//                        slice.push(intersect(a, b, k1));
//                        if (!closed) slice = newSlice(slices, slice, area, dist);
//
//                    } else if (bk > k2) { // |  ---|-->
//                        slice.push(intersect(a, b, k2));
//                        if (!closed) slice = newSlice(slices, slice, area, dist);
//                    }
//                    // | --> |
//                }
//            }
//
//            a = points[len - 1];
//            ak = a[axis];
//            if (ak >= k1 && ak <= k2) slice.push(a);
//            
//            var sliceLen = slice.length;
//            
//            // close the polygon if its endpoints are not the same after clipping
//            if (closed && slice[0] !== slice[sliceLen - 1]) slice.push(slice[0]);
//            
//            // add the final slice
//            newSlice(slices, slice, area, dist);
//        }
//        
//        return slices;

    }

    class func newSlice(slices: [Slice], slice: Slice, area: Int, dist: Int) {

        if (slice.points.count > 0) {
            slice.area = area
            slice.dist = dist
            slices.append(slice)
        }

    }

//        if (slice.length) {
//                slice.area = area;
//                slice.dist = dist;
//                slices.push(slice);
//            }
//            return [];

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

//        var sqTolerance = tolerance * tolerance,
//        len = points.length,
//        first = 0,
//        last = len - 1,
//        stack = [],
//        i, maxSqDist, sqDist, index;
//
//        // always retain the endpoints (1 is the max value)
//        points[first][2] = 1;
//        points[last][2] = 1;
//
//        // avoid recursion by using a stack
//        while (last) {
//
//            maxSqDist = 0;
//
//            for (i = first + 1; i < last; i++) {
//                sqDist = getSqSegDist(points[i], points[first], points[last]);
//
//                if (sqDist > maxSqDist) {
//                    index = i;
//                    maxSqDist = sqDist;
//                }
//            }
//
//            if (maxSqDist > sqTolerance) {
//                points[index][2] = maxSqDist; // save the point importance in squared pixels as a z coordinate
//                stack.push(first, index, index, last);
//            }
//            
//            last = stack.pop();
//            first = stack.pop();
//        }

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

//        var x = a[0], y = a[1],
//        bx = b[0], by = b[1],
//        px = p[0], py = p[1],
//        dx = bx - x,
//        dy = by - y;
//
//        if (dx !== 0 || dy !== 0) {
//
//            var t = ((px - x) * dx + (py - y) * dy) / (dx * dx + dy * dy);
//
//            if (t > 1) {
//                x = bx;
//                y = by;
//
//            } else if (t > 0) {
//                x += dx * t;
//                y += dy * t;
//            }
//        }
//        
//        dx = px - x;
//        dy = py - y;
//        
//        return dx * dx + dy * dy;

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

//        var features = [];
//
//        if (data.type === 'FeatureCollection') {
//            for (var i = 0; i < data.features.length; i++) {
//                convertFeature(features, data.features[i], tolerance);
//            }
//        } else if (data.type === 'Feature') {
//            convertFeature(features, data, tolerance);
//
//        } else {
//            // single geometry or a geometry collection
//            convertFeature(features, {geometry: data}, tolerance);
//        }
//        return features;

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

//        var geom = feature.geometry,
//        type = geom.type,
//        coords = geom.coordinates,
//        tags = feature.properties,
//        i, j, rings;
//
//        if (type === 'Point') {
//            features.push(create(tags, 1, [projectPoint(coords)]));
//
//        } else if (type === 'MultiPoint') {
//            features.push(create(tags, 1, project(coords)));
//
//        } else if (type === 'LineString') {
//            features.push(create(tags, 2, [project(coords, tolerance)]));
//
//        } else if (type === 'MultiLineString' || type === 'Polygon') {
//            rings = [];
//            for (i = 0; i < coords.length; i++) {
//                rings.push(project(coords[i], tolerance));
//            }
//            features.push(create(tags, type === 'Polygon' ? 3 : 2, rings));
//
//        } else if (type === 'MultiPolygon') {
//            rings = [];
//            for (i = 0; i < coords.length; i++) {
//                for (j = 0; j < coords[i].length; j++) {
//                    rings.push(project(coords[i][j], tolerance));
//                }
//            }
//            features.push(create(tags, 3, rings));
//
//        } else if (type === 'GeometryCollection') {
//            for (i = 0; i < geom.geometries.length; i++) {
//                convertFeature(features, {
//                    geometry: geom.geometries[i],
//                    properties: tags
//                    }, tolerance);
//            }
//            
//        } else {
//            console.warn('Unsupported GeoJSON type: ' + geom.type);
//        }

    }

    class func create(tags: [Int], type: Feature.Type, geometry: [Feature]) -> Feature {

        var feature = Feature(geometry: geometry, type: type, tags: tags, min: Point(x: 1, y: 1), max: Point(x: 0, y: 0))
        calcBBox(feature)
        return feature

//    function create(tags, type, geometry) {
//        var feature = {
//            geometry: geometry,
//            type: type,
//            tags: tags || null,
//            min: [1, 1], // initial bbox values;
//            max: [0, 0]  // note that all coords are in [0..1] range
//        };
//        calcBBox(feature);
//        return feature;
//    }

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

//    function project(lonlats, tolerance) {
//        var projected = [];
//            for (var i = 0; i < lonlats.length; i++) {
//                projected.push(projectPoint(lonlats[i]));
//            }
//            if (tolerance) {
//                simplify(projected, tolerance);
//                calcSize(projected);
//            }
//            return projected;
//    }

    }

    class func projectPoint(p: Point) -> [Int] {

        let sin = sin(Double(p.y) * M_PI / 180)
        let x = Double(p.x) / 360 + 0.5
        let y = 0.5 - 0.25 * log((1 + sin) / (1 - sin)) / M_PI

        return [x, y, 0]

//    function projectPoint(p) {
//        var sin = Math.sin(p[1] * Math.PI / 180),
//        x = (p[0] / 360 + 0.5),
//        y = (0.5 - 0.25 * Math.log((1 + sin) / (1 - sin)) / Math.PI);
//    return [x, y, 0];
//    }

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

//    // calculate area and length of the poly
//    function calcSize(points) {
//        var area = 0,
//        dist = 0;
//
//        for (var i = 0, a, b; i < points.length - 1; i++) {
//            a = b || points[i];
//            b = points[i + 1];
//
//            area += a[0] * b[1] - b[0] * a[1];
//
//            // use Manhattan distance instead of Euclidian one to avoid expensive square root computation
//            dist += Math.abs(b[0] - a[0]) + Math.abs(b[1] - a[1]);
//        }
//        points.area = Math.abs(area / 2);
//        points.dist = dist;
//    }

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

//    // calculate the feature bounding box for faster clipping later
//    function calcBBox(feature) {
//        var geometry = feature.geometry,
//        min = feature.min,
//        max = feature.max;
//
//        if (feature.type === 1) {
//            calcRingBBox(min, max, geometry);
//        } else {
//            for (var i = 0; i < geometry.length; i++) {
//                calcRingBBox(min, max, geometry[i]);
//            }
//        }
//        return feature;
//    }

    }

    class func calcRingBBox(var min: Point, var max: Point, points: [Point]) {

        for i in 0...points.count {
            let p = points[i]
            min.x = min(p.x, min.x)
            max.x = max(p.x, max.x)
            min.y = min(p.y, min.y)
            max.y = max(p.y, max.y)
        }

//    function calcRingBBox(min, max, points) {
//        for (var i = 0, p; i < points.length; i++) {
//            p = points[i];
//            min[0] = Math.min(p[0], min[0]);
//            max[0] = Math.max(p[0], max[0]);
//            min[1] = Math.min(p[1], min[1]);
//            max[1] = Math.max(p[1], max[1]);
//        }
//    }

    }

}






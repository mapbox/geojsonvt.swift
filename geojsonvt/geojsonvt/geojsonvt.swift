import Foundation

class GeoJSONVT {

    var baseZoom: Int
    var maxZoom: Int
    var maxPoints: Int
    var tolerance: Double
    var extent: Double = 4096
    var buffer: Double = 64
    var debug: Bool = false

    var tiles = [Int: Tile]()

    var stats = [String: Int]()

    var total = 0

    init(data: String, baseZoom: Int = 14, maxZoom: Int = 4, maxPoints: Int = 100, tolerance: Double = 3, debug: Bool = false) {

        self.baseZoom = baseZoom
        self.maxZoom = maxZoom
        self.maxPoints = maxPoints
        self.tolerance = tolerance
        self.debug = debug

        if (self.debug) {
            Util.time("preprocess data")
        }

        let z2 = 1 << baseZoom

        let deserializedData = NSJSONSerialization.JSONObjectWithData(data.dataUsingEncoding(NSUTF8StringEncoding,
            allowLossyConversion: false)!, options: nil, error: nil) as JSON

        let features = Convert.convert(data: deserializedData, tolerance: self.tolerance / (Double(z2) * self.extent))

        if (self.debug) {
            Util.timeEnd("preprocess data")
            Util.time("generate tiles up to z\(maxZoom)")
        }

        self.splitTile(features: features, z: 0, x: 0, y: 0)

        if (self.debug) {
            NSLog("features: %i, points: %i", self.tiles[0]!.numFeatures, self.tiles[0]!.numPoints)
            Util.timeEnd("generate tiles up to z\(maxZoom)")
            NSLog("tiles generated: %i %@", self.total, self.stats)
        }
    }

    func splitTile(#features: [ProjectedFeature], z: Int, x: Int, y: Int, cz: Int = -1, cx: Int = -1, cy: Int = -1) {

        var stack: [Any] = [features, z, x, y]

        while (stack.count > 0) {
            let features = stack.removeAtIndex(0) as [ProjectedFeature]
            let z = stack.removeAtIndex(0) as Int
            let x = stack.removeAtIndex(0) as Int
            let y = stack.removeAtIndex(0) as Int

            let z2 = 1 << z
            let id = toID(z: z, x: x, y: y)
            var tile: Tile
            let tileTolerance = (z == self.baseZoom ? 0 : self.tolerance / (Double(z2) * self.extent))

            if (self.tiles.indexForKey(id) != nil) {
                tile = self.tiles[id]!
            } else {
                if (self.debug) {
                    Util.time("creation")
                }

                tile = Tile.createTile(features: features, z2: z2, tx: x, ty: y, tolerance: tileTolerance,
                    extent: extent, noSimplify: (z == self.baseZoom))

                self.tiles[id] = tile

                if (self.debug) {
                    NSLog("tile z%i-%i-%i (features: %i, points: %i, simplified: %i", z, x, y,
                        tile.numFeatures, tile.numPoints, tile.numSimplified)
                    Util.timeEnd("creation")

                    let key = "z\(z):"
                    if (self.stats.count - 1 >= z) {
                        var value = self.stats[key]! + 1
                        self.stats[key] = value
                    } else {
                        self.stats[key] = 1
                    }
                    self.total++
                }
            }

            if (cz < 0 && (z == self.maxZoom || tile.numPoints <= self.maxPoints ||
                self.isClippedSquare(features: tile.features, extent: self.extent, buffer: self.buffer)) || z == self.baseZoom || z == cz) {
                    tile.source = features
                    continue
            }

            if (cz >= 0) {
                tile.source = features
            } else {
                tile.source = []
            }

            if (self.debug) {
                Util.time("clipping")
            }

            let k1 = 0.5 * self.buffer / self.extent
            let k2 = 0.5 - k1
            let k3 = 0.5 + k1
            let k4 = 1 + k1

            var tl = [ProjectedFeature]()
            var bl = [ProjectedFeature]()
            var tr = [ProjectedFeature]()
            var br = [ProjectedFeature]()
            var left = [ProjectedFeature]()
            var right = [ProjectedFeature]()
            var m: Int = 0
            var goLeft = false
            var goTop = false

            if (cz >= 0) {
                m = 1 << (cz - z)
                goLeft = Double(cx) / Double(m) - Double(x) < 0.5
                goTop  = Double(cy) / Double(m) - Double(y) < 0.5
            }

            if (cz < 0 || goLeft) {
                left = Clip.clip(features: features, scale: z2, k1: Double(x) - k1, k2: Double(x) + k3, axis: 0, intersect: intersectX)
            }

            if (cz < 0 || !goLeft) {
                right = Clip.clip(features: features, scale: z2, k1: Double(x) + k2, k2: Double(x) + k4, axis: 0, intersect: intersectX)
            }

            if (left.count > 0) {
                if (cz < 0 || goTop) {
                    tl = Clip.clip(features: left, scale: z2, k1: Double(y) - k1, k2: Double(y) + k3, axis: 1, intersect: intersectY)
                }

                if (cz < 0 || !goTop) {
                    bl = Clip.clip(features: left, scale: z2, k1: Double(y) + k2, k2: Double(y) + k4, axis: 1, intersect: intersectY)
                }
            }

            if (right.count > 0) {
                if (cz < 0 || goTop) {
                    tr = Clip.clip(features: right, scale: z2, k1: Double(y) - k1, k2: Double(y) + k3, axis: 1, intersect: intersectY)
                }

                if (cz < 0 || !goTop) {
                    br = Clip.clip(features: right, scale: z2, k1: Double(y) + k2, k2: Double(y) + k4, axis: 1, intersect: intersectY)
                }
            }

            if (self.debug) {
                Util.timeEnd("clipping")
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
            NSLog("drilling down to z%i-%i-%i", z, x, y)
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
            NSLog("found parent tile z%i-%i-%i", z0, x0, y0)
        }

        if (parent!.source.count > 0) {
            if (self.isClippedSquare(features: parent!.features, extent: self.extent, buffer: self.buffer)) {
                return parent!
            }

            if (self.debug) {
                Util.time("drilling down")
            }

            self.splitTile(features: parent!.source, z: z0, x: x0, y: y0, cz: z, cx: x, cy: y)

            if (self.debug) {
                Util.timeEnd("drilling down")
            }
        }

        return self.tiles[id]!
    }

    func isClippedSquare(#features: [TileFeature], extent: Double, buffer: Double) -> Bool {

        if (features.count != 1) {
            return false
        }

        let feature = features.first!

        if (feature.type != .Polygon || feature.geometry.count > 1) {
            return false
        }

        for i in 0...(feature.geometry.first! as TileRing).points.count {
            let p = (feature.geometry.first! as TileRing).points[i] as TilePoint
            if ((Double(p.x) != -self.buffer && Double(p.x) != self.extent + self.buffer) ||
                (Double(p.y) != -self.buffer && Double(p.y) != self.extent + self.buffer)) {
                    return false
            }
        }

        return true
    }

    func toID(#z: Int, x: Int, y: Int) -> Int {

        return (((1 << z) * y + x) * 32) + z;
    }

    func intersectX(a: ProjectedPoint, b: ProjectedPoint, x: Double) -> ProjectedPoint {

        let r1 = x
        let r2 = (x - a.x) * (b.y - a.y) / (b.x - a.x) + a.y
        let r3: Double = 1

        return ProjectedPoint(x: r1, y: r2, z: r3);
    }
    
    func intersectY(a: ProjectedPoint, b: ProjectedPoint, y: Double) -> ProjectedPoint {

        let r1 = (y - a.y) * (b.x - a.x) / (b.y - a.y) + a.x
        let r2 = y
        let r3: Double = 1
        
        return ProjectedPoint(x: r1, y: r2, z: r3);
    }
    
    func extend(var #dest: [ProjectedPoint], src: [ProjectedPoint]) -> [ProjectedPoint] {

        for i in 0..<src.count {
            dest[i] = src[i]
        }
        return dest
    }
    
}

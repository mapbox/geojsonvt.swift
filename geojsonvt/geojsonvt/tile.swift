import Foundation

class Tile {

    var features = [TileFeature]()
    var numPoints: Int = 0
    var numSimplified: Int = 0
    var numFeatures: Int = 0
    var source = [ProjectedFeature]()

    init() {

    }

    class func createTile(#features: [ProjectedFeature], z2: Int, tx: Int, ty: Int,
        tolerance: Double, extent: Double, noSimplify: Bool) -> Tile {

        var tile = Tile()

        for i in 0..<features.count {
            tile.numFeatures++
            Tile.addFeature(tile: &tile, feature: features[i], z2: z2, tx: tx, ty: ty, tolerance: tolerance,
                extent: extent, noSimplify: noSimplify)
        }

        return tile
    }

    class func addFeature(inout #tile: Tile, feature: ProjectedFeature, z2: Int, tx: Int, ty: Int,
        tolerance: Double, extent: Double, noSimplify: Bool) {

        let geom = feature.geometry as! ProjectedGeometryContainer
        let type = feature.type
        var transformed = [TileGeometry]()
        let sqTolerance = tolerance * tolerance
        var ring: ProjectedGeometryContainer

        if (type == .Point) {
            for i in 0..<geom.members.count {
                let p = geom.members[i] as! ProjectedPoint
                transformed.append(Tile.transformPoint(p, z2: z2, tx: tx, ty: ty, extent: extent))
                tile.numPoints++
                tile.numSimplified++
            }
        } else {
            for i in 0..<geom.members.count {
                ring = geom.members[i] as! ProjectedGeometryContainer

                if (!noSimplify && ((type == .LineString && ring.dist < tolerance) ||
                    (type == .Polygon && ring.area < sqTolerance))) {
                    tile.numPoints += ring.members.count
                    continue
                }

                var transformedRing = TileRing()

                for j in 0..<ring.members.count {
                    let p = ring.members[j] as! ProjectedPoint
                    if (noSimplify || p.z > Double(sqTolerance)) {
                        let transformedPoint = Tile.transformPoint(p, z2: z2, tx: tx, ty: ty, extent: extent)
                        transformedRing.points.append(transformedPoint)
                        tile.numSimplified++
                    }
                    tile.numPoints++
                }

                transformed.append(transformedRing)
            }
        }

        if (transformed.count > 0) {
            tile.features.append(TileFeature(geometry: transformed, type: type, tags: feature.tags))
        }
    }

    class func transformPoint(p: ProjectedPoint, z2: Int, tx: Int, ty: Int, extent: Double) -> TilePoint {

        let x = Int(extent * (p.x * Double(z2) - Double(tx)))
        let y = Int(extent * (p.y * Double(z2) - Double(ty)))

        return TilePoint(x: x, y: y)
    }
    
}

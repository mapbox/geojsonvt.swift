class LonLat {

    let lon: Double
    let lat: Double

    init(coordinates: [Double]) {
        self.lon = coordinates[0]
        self.lat = coordinates[1]
    }

}

class ProjectedGeometry {
    // abstract
}

class ProjectedPoint: ProjectedGeometry {

    var x: Double
    var y: Double
    var z: Double

    init(x: Double = -1, y: Double = -1, z: Double = -1) {
        self.x = x
        self.y = y
        self.z = z
    }

    func isValid() -> Bool {
        return (x >= 0 && y >= 0 && z >= 0)
    }

    func isEqualToPoint(p2: ProjectedPoint) -> Bool {
        return (self.x == p2.x && self.y == p2.y && self.z == p2.z)
    }

}

class ProjectedGeometryContainer: ProjectedGeometry {

    var members: [ProjectedGeometry]
    var isPolyline = false
    var area: Double = 0
    var dist: Double = 0

    init(members: [ProjectedGeometry] = []) {
        self.members = members
    }

}

enum ProjectedFeatureType: Int {
    case Point = 1
    case LineString
    case Polygon
}

class ProjectedFeature {

    let geometry: ProjectedGeometry
    let type: ProjectedFeatureType
    let tags: [String]
    let minPoint: ProjectedPoint
    let maxPoint: ProjectedPoint

    init(geometry: ProjectedGeometry, type: ProjectedFeatureType, tags: [String]) {
        self.geometry = geometry
        self.type = type
        self.tags = tags
        self.minPoint = ProjectedPoint(x: 1, y: 1, z: 0)
        self.maxPoint = ProjectedPoint(x: 0, y: 0, z: 0)
        Convert.calcBBox(feature: self)
    }
}

class TileGeometry {
    // abstract
}

class TilePoint: TileGeometry {

    var x: Int
    var y: Int

    init(x: Int, y: Int) {
        self.x = x
        self.y = y
    }

}

class TileRing: TileGeometry {

    var points: [TilePoint]

    init(points: [TilePoint] = []) {
        self.points = points
    }

}

typealias TileFeatureType = ProjectedFeatureType

class TileFeature {

    let geometry: [TileGeometry]
    let type: TileFeatureType
    let tags: [String]

    init(geometry: [TileGeometry], type: TileFeatureType, tags: [String]) {
        self.geometry = geometry
        self.type = type
        self.tags = tags
    }

}

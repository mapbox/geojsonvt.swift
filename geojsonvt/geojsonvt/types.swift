import Foundation

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

    var x: Double = -1
    var y: Double = -1
    var z: Double = -1

    init(x: Double = -1, y: Double = -1, z: Double = -1) {
        self.x = x
        self.y = y
        self.z = z
    }

    func isValid() -> Bool {
        return (x >= 0 && y >= 0 && z >= 0)
    }

}

func != (left: ProjectedPoint, right: ProjectedPoint) -> Bool {
    return (left.x != right.x || left.y != right.y || left.z != right.z)
}

class ProjectedGeometryContainer: ProjectedGeometry {

    var members = [ProjectedGeometry]()
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
    let tags: Tags
    let minPoint: ProjectedPoint
    let maxPoint: ProjectedPoint

    init(geometry: ProjectedGeometry, type: ProjectedFeatureType, tags: Tags) {
        self.geometry = geometry
        self.type = type
        self.tags = tags
        self.minPoint = ProjectedPoint(x: 1, y: 1, z: 0)
        self.maxPoint = ProjectedPoint(x: 0, y: 0, z: 0)
    }
}

class TileGeometry {
    // abstract
}

class TilePoint: TileGeometry {

    let x: Int
    let y: Int

    init(x: Int, y: Int) {
        self.x = x
        self.y = y
    }

}

class TileRing: TileGeometry {

    var points = [TilePoint]()

    init(points: [TilePoint] = []) {
        self.points = points
    }

}

typealias TileFeatureType = ProjectedFeatureType

typealias Tags = Dictionary<String, String>

class TileFeature {

    let geometry: [TileGeometry]
    let type: TileFeatureType
    let tags: Tags

    init(geometry: [TileGeometry], type: TileFeatureType, tags: Tags) {
        self.geometry = geometry
        self.type = type
        self.tags = tags
    }

}

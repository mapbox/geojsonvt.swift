import Foundation

typealias JSON = Dictionary<String, AnyObject>

class Convert {

    class func convert(#data: JSON, tolerance: Double) -> [ProjectedFeature] {

        var features = [ProjectedFeature]()

        if ((data["type"] as String) == "FeatureCollection") {

            let count = (data["features"] as [JSON]).count
            NSLog("there are %i total features to convert", count)

            for i in 0..<(data["features"] as [JSON]).count {
                Convert.convertFeature(features: &features, feature: (data["features"] as [JSON])[i], tolerance: tolerance)
            }
        } else if (data["type"] as String == "Feature") {
            Convert.convertFeature(features: &features, feature: data, tolerance: tolerance)
        } else {
            Convert.convertFeature(features: &features, feature: ["geometry": data], tolerance: tolerance)
        }

        return features
    }

    class func convertFeature(inout #features: [ProjectedFeature], feature: JSON, tolerance: Double) {

        let geom = feature["geometry"] as JSON
        let type = geom["type"] as String
        let tags = feature["properties"] as Tags

        if (type == "Point") {

            let coordinates = geom["coordinates"] as [Double]
            let point = Convert.projectPoint(LonLat(coordinates: coordinates))

            let geometry = ProjectedGeometryContainer(members: [point])

            features.append(Convert.create(tags: tags, type: ProjectedFeatureType.Point, geometry: geometry))

        } else if (type == "MultiPoint") {

            let coordinatePairs = geom["coordinates"] as [[Double]]
            var points = [LonLat]()
            for coordinatePair in coordinatePairs {
                points.append(LonLat(coordinates: coordinatePair))
            }

            let geometry = Convert.project(lonlats: points)

            features.append(Convert.create(tags: tags, type: ProjectedFeatureType.Point, geometry: geometry))

        } else if (type == "LineString") {

            let coordinatePairs = geom["coordinates"] as [[Double]]
            var points = [LonLat]()
            for coordinatePair in coordinatePairs {
                points.append(LonLat(coordinates: coordinatePair))
            }

            let geometry = ProjectedGeometryContainer(members: [Convert.project(lonlats: points, tolerance: tolerance)])

            features.append(Convert.create(tags: tags, type: ProjectedFeatureType.LineString, geometry: geometry))

        } else if (type == "MultiLineString" || type == "Polygon") {

            var rings = ProjectedGeometryContainer()
            let lines = geom["coordinates"] as [[[Double]]]
            for line in lines {
                var points = [LonLat]()
                for coordinatePair in line {
                    points.append(LonLat(coordinates: coordinatePair))
                }
                let ring = Convert.project(lonlats: points, tolerance: tolerance)
                rings.members.append(ring)
            }

            let projectedType = (type == "Polygon" ?
                ProjectedFeatureType.Polygon :
                ProjectedFeatureType.LineString)

            let geometry = rings

            features.append(Convert.create(tags: tags, type: projectedType, geometry: geometry))

            NSLog("features now has %i items", features.count)

        } else if (type == "MultiPolygon") {

            let rings = ProjectedGeometryContainer()
            let polygons = geom["coordinates"] as [[[[Double]]]]
            for polygon in polygons {
                for line in polygon {
                    var points = [LonLat]()
                    for coordinatePair in line {
                        points.append(LonLat(coordinates: coordinatePair))
                    }
                    let ring = Convert.project(lonlats: points, tolerance: tolerance)
                    rings.members.append(ring)
                }
            }

            let geometry = rings

            features.append(Convert.create(tags: tags, type: ProjectedFeatureType.Polygon, geometry: geometry))

        } else if (type == "GeometryCollection") {

            let geometries = geom["geometries"] as [JSON]
            for geometry in geometries {
                Convert.convertFeature(features: &features, feature: geometry, tolerance: tolerance)
            }

        } else {

            let geometryType = geom["type"] as String
            println("Unsupported GeoJSON type: \(geometryType)")

        }
    }

    class func create(#tags: Tags, type: ProjectedFeatureType, geometry: ProjectedGeometry) -> ProjectedFeature {

        var feature = ProjectedFeature(geometry: geometry, type: type, tags: tags)
        Convert.calcBBox(feature: &feature)

        return feature
    }

    class func project(#lonlats: [LonLat], tolerance: Double = 0) -> ProjectedGeometryContainer {

        var projected = ProjectedGeometryContainer()
        for i in 0..<lonlats.count {
            projected.members.append(Convert.projectPoint(lonlats[i]))
        }
        if (tolerance > 0) {
            Simplify.simplify(points: projected, tolerance: tolerance)
            Convert.calcSize(geometryContainer: &projected)
        }

        return projected
    }

    class func projectPoint(p: LonLat) -> ProjectedPoint {

        let sine = sin(p.lat * M_PI / 180)
        let x = p.lon / 360 + 0.5
        let y = 0.5 - 0.25 * log((1 + sine) / (1 - sine)) / M_PI

        return ProjectedPoint(x: x, y: y, z: 0)
    }

    class func calcSize(inout #geometryContainer: ProjectedGeometryContainer) {

        var area: Double = 0
        var dist: Double = 0
        var a = ProjectedPoint()
        var b = ProjectedPoint()

        for i in 0..<(geometryContainer.members.count - 1) {
            a = (b.isValid() ? b : geometryContainer.members[i] as ProjectedPoint)
            b = geometryContainer.members[i + 1] as ProjectedPoint

            area += a.x * b.y - b.x * a.y
            dist += abs(b.x - a.x) + abs(b.y - a.y)
        }

        geometryContainer.area = abs(area / 2)
        geometryContainer.dist = dist
    }

    class func calcBBox(inout #feature: ProjectedFeature) {

        let geometry = feature.geometry
        var minPoint = feature.minPoint
        var maxPoint = feature.maxPoint

        if (feature.type == ProjectedFeatureType.Point) {
            Convert.calcRingBBox(minPoint: &minPoint, maxPoint: &maxPoint, geometry: geometry as ProjectedGeometryContainer)
        } else {
            for i in 0..<(geometry as ProjectedGeometryContainer).members.count {
                let featureGeometry = (geometry as ProjectedGeometryContainer).members[i] as ProjectedGeometryContainer
                Convert.calcRingBBox(minPoint: &minPoint, maxPoint: &maxPoint, geometry: featureGeometry)
            }
        }
    }

    class func calcRingBBox(inout #minPoint: ProjectedPoint, inout maxPoint: ProjectedPoint, geometry: ProjectedGeometryContainer) {

        for i in 0..<geometry.members.count {
            let p = geometry.members[i] as ProjectedPoint
            minPoint.x = min(p.x, minPoint.x)
            maxPoint.x = max(p.x, maxPoint.x)
            minPoint.y = min(p.y, minPoint.y)
            maxPoint.y = max(p.y, maxPoint.y)
        }
    }
    
}

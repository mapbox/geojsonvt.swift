import Foundation

class Clip {

    class func clip(#features: [ProjectedFeature], scale: Int, var k1: Double, var k2: Double, axis: Int,
        intersect: (ProjectedPoint, ProjectedPoint, Double) -> ProjectedPoint) -> [ProjectedFeature] {

        k1 /= Double(scale)
        k2 /= Double(scale)

        var clipped = [ProjectedFeature]()

        for i in 0..<features.count {

            let feature = features[i]
            let geometry = feature.geometry
            let type = feature.type
            var min: Double = 0
            var max: Double = 0

            if (feature.minPoint.isValid()) {
                min = (axis == 0 ? feature.minPoint.x : feature.minPoint.y)
                max = (axis == 0 ? feature.maxPoint.x : feature.maxPoint.y)

                if (min >= Double(k1) && max <= Double(k2)) {
                    clipped.append(feature)
                    continue
                } else if (min > Double(k2) || max < Double(k1)) {
                    continue
                }
            }

            var slices = ProjectedGeometryContainer()

            if (type == .Point) {
                slices = Clip.clipPoints(geometry: geometry as ProjectedGeometryContainer, k1: k1, k2: k2, axis: axis)
            } else {
                slices = Clip.clipGeometry(geometry: geometry as ProjectedGeometryContainer, k1: k1, k2: k2, axis: axis,
                    intersect: intersect, closed: (type == ProjectedFeatureType.Polygon))
            }

            if (slices.members.count > 0) {
                clipped.append(ProjectedFeature(geometry: slices, type: type, tags: features[i].tags))
            }
        }

        return clipped
    }

    class func clipPoints(#geometry: ProjectedGeometryContainer, k1: Double, k2: Double, axis: Int) -> ProjectedGeometryContainer {

        var slice = ProjectedGeometryContainer()

        for i in 0..<geometry.members.count {
            let a = geometry.members[i] as ProjectedPoint
            let ak = (axis == 0 ? a.x : a.y)

            if (ak >= k1 && ak <= k2) {
                slice.members.append(a)
            }
        }

        return slice
    }

    class func clipGeometry(#geometry: ProjectedGeometryContainer, k1: Double, k2: Double, axis: Int,
        intersect: (ProjectedPoint, ProjectedPoint, Double) -> ProjectedPoint, closed: Bool) -> ProjectedGeometryContainer {

        var slices = ProjectedGeometryContainer()

        for i in 0..<geometry.members.count {

            var ak: Double = 0
            var bk: Double = 0
            var b = ProjectedPoint()
            let points = geometry.members[i] as ProjectedGeometryContainer
            let area = points.area
            let dist = points.dist
            let len = points.members.count
            var a = ProjectedPoint()

            var slice = ProjectedGeometryContainer()

            for j in 0..<(len - 1) {
                a = (b.isValid() ? b : points.members[j] as ProjectedPoint)
                b = points.members[j + 1] as ProjectedPoint
                ak = (bk > 0 ? bk : (axis == 0 ? a.x : a.y))
                bk = (axis == 0 ? b.x : b.y)

                if (ak < k1) {
                    if (bk > k2) {
                        slice.members.append(intersect(a, b, k1))
                        slice.members.append(intersect(a, b, k2))
                        if (!closed) {
                            slice = newSlice(slices: &slices, slice: slice, area: area, dist: dist)
                        }
                    } else if (bk >= k1) {
                        slice.members.append(intersect(a, b, k1))
                    }
                } else if (ak > k2) {
                    if (bk < k1) {
                        slice.members.append(intersect(a, b, k2))
                        slice.members.append(intersect(a, b, k1))
                        if (!closed) {
                            slice = newSlice(slices: &slices, slice: slice, area: area, dist: dist)
                        }
                    } else if (bk <= k2) {
                        slice.members.append(intersect(a, b, k2))
                    }
                } else {
                    slice.members.append(a)

                    if (bk < k1) {
                        slice.members.append(intersect(a, b, k1))
                        if (!closed) {
                            slice = newSlice(slices: &slices, slice: slice, area: area, dist: dist)
                        }
                    } else if (bk > k2) {
                        slice.members.append(intersect(a, b, k2))
                        if (!closed) {
                            slice = newSlice(slices: &slices, slice: slice, area: area, dist: dist)
                        }
                    }
                }
            }

            a = points.members[len - 1] as ProjectedPoint
            ak = (axis == 0 ? a.x : a.y)

            if (ak >= k1 && ak <= k2) {
                slice.members.append(a)
            }

            if (closed && slice.members.count > 0) {
                let first = slice.members[0] as ProjectedPoint
                let last  = slice.members[slice.members.count - 1] as ProjectedPoint
                if (!first.isEqualToPoint(last)) {
                    slice.members.append(ProjectedPoint(x: first.x, y: first.y, z: first.z))
                }
            }

            newSlice(slices: &slices, slice: slice, area: area, dist: dist)
        }

        return slices
    }
    
    class func newSlice(inout #slices: ProjectedGeometryContainer, slice: ProjectedGeometryContainer,
        area: Double, dist: Double) -> ProjectedGeometryContainer {
        
        if (slice.members.count > 0) {
            slice.area = area
            slice.dist = dist
            slices.members.append(slice)
        }
        
        return ProjectedGeometryContainer()
    }
}

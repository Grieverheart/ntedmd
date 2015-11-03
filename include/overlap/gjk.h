#ifndef EDMD_OVERLAP_GJK_H
#define EDMD_OVERLAP_GJK_H

#include "clam.h"

struct Particle;
namespace shape{
    class Convex;
}

namespace overlap{
    clam::Vec3d gjk_distance(const Particle&, const shape::Convex&, const Particle&, const shape::Convex&, double feather = 0);
    bool gjk_boolean(const Particle&, const shape::Convex&, const Particle&, const shape::Convex&, double feather = 0);
    bool gjk_raycast(const Particle&, const shape::Convex&, const Particle&, const shape::Convex&, const clam::Vec3d& ray_dir, double& t, clam::Vec3d& normal);
    //TODO: Think about returning closest point on A and distance vector.
    clam::Vec3d gjk_closest_points(const Particle&, const shape::Convex&, const Particle&, const shape::Convex&, clam::Vec3d& pa, clam::Vec3d& pb);

}

#endif

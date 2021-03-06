#include "overlap/overlap.h"
#include "particle.h"
#include "shape/variant.h"
#include "overlap/gjk.h"
#include "overlap/ray_casting.h"

namespace overlap{

    namespace{

        class ShapeDistanceVisitor: public boost::static_visitor<clam::Vec3d> {
        public:
            ShapeDistanceVisitor(const Transform& pa, const Transform& pb):
                pa_(pa), pb_(pb)
            {}

            template<typename T, typename U>
            clam::Vec3d operator()(const T& a, const U& b)const{
                return gjk_distance(pa_, a, pb_, b);
            }

        private:
            const Transform& pa_;
            const Transform& pb_;
        };

        template<>
        inline clam::Vec3d ShapeDistanceVisitor::operator()(const shape::Sphere& a, const shape::Sphere& b)const{
            return pb_.pos_ - pa_.pos_ - pa_.size_ * a.radius() - pb_.size_ * b.radius();
        }

        //TODO: Move to separate header file
        template<typename T>
        static inline T sqr(T val){
            return val * val;
        }

        class ShapeOverlapVisitor: public boost::static_visitor<bool> {
        public:
            ShapeOverlapVisitor(const Transform& pa, const Transform& pb, double feather):
                pa_(pa), pb_(pb), feather_(feather)
            {}

            template<typename T, typename U>
            bool operator()(const T& a, const U& b)const{
                return gjk_boolean(pa_, a, pb_, b, feather_);
            }

        private:
            const Transform& pa_;
            const Transform& pb_;
            double feather_;
        };

        template<>
        inline bool ShapeOverlapVisitor::operator()(const shape::Sphere& a, const shape::Sphere& b)const{
            return (pb_.pos_ - pa_.pos_).length2() < sqr(pa_.size_ * a.radius() + pb_.size_ * b.radius() + feather_);
        }

        class ShapeRaycastVisitor: public boost::static_visitor<bool> {
        public:
            ShapeRaycastVisitor(const Transform& pa, const Transform& pb, const clam::Vec3d& ray_dir, double& distance, clam::Vec3d& normal):
                pa_(pa), pb_(pb), ray_dir_(ray_dir), dist_(distance), normal_(normal)
            {}

            template<typename T, typename U>
            bool operator()(const T& a, const U& b)const{
                return gjk_raycast(pa_, a, pb_, b, ray_dir_, dist_, normal_);
            }

        private:
            const Transform& pa_;
            const Transform& pb_;
            const clam::Vec3d& ray_dir_;
            double& dist_;
            clam::Vec3d& normal_;
        };

        template<>
        inline bool ShapeRaycastVisitor::operator()(const shape::Sphere& a, const shape::Sphere& b)const{
            return sphere_raycast(pa_.size_ * a.radius() + pb_.size_ * b.radius(), pa_.pos_ - pb_.pos_, ray_dir_, dist_, &normal_);
        }
    }//namespace detail

    clam::Vec3d shape_distance(const Transform& pa, const shape::Variant& a, const Transform& pb, const shape::Variant& b){
        return boost::apply_visitor(ShapeDistanceVisitor(pa, pb), a, b);
    }

    bool shape_overlap(const Transform& pa, const shape::Variant& a, const Transform& pb, const shape::Variant& b, double feather){
        return boost::apply_visitor(ShapeOverlapVisitor(pa, pb, feather), a, b);
    }

    bool shape_raycast(const Transform& pa, const shape::Variant& a, const Transform& pb, const shape::Variant& b, const clam::Vec3d& ray_dir, double& distance, clam::Vec3d& normal){
        return boost::apply_visitor(ShapeRaycastVisitor(pa, pb, ray_dir, distance, normal), a, b);
    }

}

#include "configuration.h"
#include "shape/variant.h"
#include "serialization/common.h"
#include "serialization/vector.h"

class SerializationVisitor: public boost::static_visitor<>{
public:
    SerializationVisitor(Archive& ar):
        ar_(ar)
    {}

    void operator()(const shape::Polyhedron& poly)const{
        int shape_type = shape::POLYHEDRON;
        ar_.write(&shape_type, sizeof(shape_type));
        serialize(ar_, poly);
    }

    void operator()(const shape::Sphere& sph)const{
        int shape_type = shape::SPHERE;
        ar_.write(&shape_type, sizeof(shape_type));
        serialize(ar_, sph);
    }

private:
    Archive& ar_;
};

Configuration::~Configuration(void){
    for(auto shape: shapes_) delete shape;
}

Configuration::Configuration(const RectangularPBC& pbc, const std::vector<Particle>&  particles, const std::vector<shape::Variant*>& shapes):
    particles_(particles),
    pbc_(pbc)
{
    shapes_.resize(shapes.size());
    for(size_t i = 0; i < shapes.size(); ++i) shapes_[i] = new shape::Variant(*shapes[i]);
}

//TODO: Unsafe when config is not initialized
Configuration::Configuration(const Configuration& other):
    particles_(other.particles_),
    pbc_(other.pbc_)
{
    shapes_.resize(other.shapes_.size());
    for(size_t i = 0; i < other.shapes_.size(); ++i) shapes_[i] = new shape::Variant(*other.shapes_[i]);
}

void serialize(Archive& ar, const Configuration& config){
    serialize(ar, config.particles_);
    serialize(ar, config.shapes_.size());
    for(auto shape: config.shapes_) boost::apply_visitor(SerializationVisitor(ar), *shape);
    serialize(ar, config.pbc_);
}

void deserialize(Archive& ar, Configuration* config){
    using size_type = decltype(config->shapes_)::size_type;

    deserialize(ar, &config->particles_);
    size_type shapes_size = 0;
    deserialize(ar, &shapes_size);
    config->shapes_.resize(shapes_size);
    for(auto& shape_ptr: config->shapes_){
        int shape_type;
        deserialize(ar, &shape_type);
        switch(shape_type){
            case shape::POLYHEDRON:{
                shape::Polyhedron poly;
                deserialize(ar, &poly);
                shape_ptr = new shape::Variant(poly);
            } break;
            case shape::SPHERE:{
                shape::Sphere sph;
                deserialize(ar, &sph);
                shape_ptr = new shape::Variant(sph);
            } break;
            //TODO: Implement panic handling.
            default: break;
        }
    }
    deserialize(ar, &config->pbc_);
}

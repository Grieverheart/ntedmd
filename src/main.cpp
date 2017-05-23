#include <cstdio>
#include <cstring>
#include <unistd.h>
#include "shape/variant.h"
#include "Simulation.h"
#include "obj_loader.h"
#include "io/config_xml.h"
#include <boost/variant/polymorphic_get.hpp>

//TODO: move these to separate file.
inline void stream_position(Particle& particle, double time){
    particle.xform.pos_ += particle.vel * (time - particle.time);
}

inline void stream_rotation(Particle& particle, double time){
    clam::Vec3d ha = ((time - particle.time) * 0.5) * particle.ang_vel;
    double l = ha.length(); // magnitude
    if(l > 0.0){
        double sl = sin(l);
        double cl = cos(l);
        particle.xform.rot_ = clam::Quatd(ha * (sl / l), cl) * particle.xform.rot_;
    }
}

inline void update_particle(Particle& particle, double time, const RectangularPBC& pbc){
    if(particle.time < time){
        stream_position(particle, time);
        stream_rotation(particle, time);
        particle.xform.pos_ = pbc.apply(particle.xform.pos_);
        particle.time = time;
    }
}

double packing_fraction(const Configuration& config){
    auto box_size = config.pbc_.getSize();
    double volume = box_size[0] * box_size[1] * box_size[2];

    double particle_volume = 0.0;
    for(auto p: config.particles_){
        particle_volume += pow(p.xform.size_, 3.0) * boost::polymorphic_get<shape::Convex>(config.shapes_[p.shape_id])->volume();
    }

    //NOTE: Assume particle unit volume!
    return particle_volume / volume;
}

int main(int argc, char *argv[]){

    Simulation* sim;
    const double check_delta = 0.1;
    const double output_delta = 10.0;
    double start_time = 0.0;

    const auto& directory = argv[3];

    if(strcmp(argv[1], "-r") == 0){
        sim = new Simulation();
        sim->load(argv[2]);

        start_time = sim->time();
    }
    else{
        Configuration config;
        xml_load_config(argv[1], config);

        //Scale box and positions for reaching target_pf
        {
            double target_pf = atof(argv[2]);
            double pf = packing_fraction(config);
            double scale = pow(pf / target_pf, 1.0 / 3.0);
            config.pbc_.setSize(scale * config.pbc_.getSize());
            for(auto& particle: config.particles_) particle.xform.pos_ = scale * particle.xform.pos_;
        }

        sim = new Simulation(config);
    }

    double pf = packing_fraction(sim->configuration());

    PeriodicCallback output(start_time + check_delta);
    output.setNextFunction([check_delta](double time){
        return time + check_delta;
    });

    FILE* pressure_fp;
    {
        char fp_buff[64];
        sprintf(fp_buff, "%s/pressure.pf%.3f.pid%u.dat", directory, pf, getpid());
        pressure_fp = fopen(fp_buff, "w");
    }

    output.setCallback([sim, pf, pressure_fp, directory, output_delta, start_time](double time) -> bool {
        static int nFiles = 0;
        static double next_output = start_time;

        printf("%f\n", time);

        //Check for overlaps
        if(sim->check_overlaps()){
            printf("Overlaps found!\n");

            //Kind of a hack >_<
            sim->~Simulation();
            new (sim) Simulation();

            char buff[128];
            sprintf(buff, "%s/archive.pf%.3f.pid%u.bin", directory, pf, getpid());
            printf("Resetting...\n");

            sim->load(buff);
            sim->restart();

            printf("Starting at: %f\n", sim->time());

            return false;
        }

        Configuration config = sim->configuration();

        for(auto& particle: config.particles_){
            update_particle(particle, time, config.pbc_);
        }

        if(time >= next_output){
            char buff[128];
            sprintf(buff, "%s/pid%u.pf%.3f.step%06u.xml", directory, getpid(), pf, nFiles);
            xml_save_config(buff, config);

            fprintf(pressure_fp, "%e\t%e\n", sim->time(), sim->average_pressure());
            fflush(pressure_fp);

            sim->reset_statistics();

            ++nFiles;
            next_output = time + output_delta;
        }

        char buff[128];
        sprintf(buff, "%s/archive.pf%.3f.pid%u.bin", directory, pf, getpid());
        sim->save(buff);

        return true;
    });

    sim->run(100000.0, output);

    fclose(pressure_fp);

    return 0;
}

/*
 * Modified version of code at:
 * https://gist.github.com/Mathiasb17/5fb84d6be9745270f66a4dec6be9b9da
 */

#include <boost/filesystem.hpp>
#include <cstdio>
#include <fstream>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <string>
#include <vector>
#include <openvdb/openvdb.h>
#include <openvdb/tools/ParticlesToLevelSet.h>
#include <openvdb/Exceptions.h>
#include <openvdb/Types.h>
#include <openvdb/tree/LeafNode.h>
#include <openvdb/tools/Filter.h>
#include <openvdb/tools/LevelSetFilter.h>
#include <openvdb/tools/VolumeToMesh.h>

DEFINE_string(input, "", "Input file containing particle data");
DEFINE_string(output, "", "Output file prefix to save signed distance to");

namespace {
  typedef openvdb::Real Real;
  typedef openvdb::Vec3R  Vec3R;
  typedef openvdb::Vec3d  Vec3d;
  typedef openvdb::Coord Coord;
  typedef openvdb::CoordBBox CoordBBox;
}  // namespace anonymous

class ParticlesList {
  public:
    typedef Vec3R PosType;

  protected:
    // Each particle has a position, velocity and radius.
    struct Particle
    {
      Vec3R p, v;
      Real r;
    };

    // Parameters to control output.
    Real radius_scale_  = 1.0;
    Real velocity_scale_ = 1.0;

    // Particle list stored in a vector.
    std::vector<Particle> particle_list_;

  public:
    ParticlesList(Real r_scale=1, Real v_scale=1)
      : radius_scale_(r_scale), velocity_scale_(v_scale) {}

    // Add particle with given position and radius.
    void AddParticle(const Vec3R &p, const Vec3R &v, const Real r = 3.0) {
      Particle part;
      part.p = p;
      part.v = v;
      part.r = r;
      particle_list_.push_back(part);
    }

    // Methods required by OpenVDB
    size_t size() const { return particle_list_.size(); }
    void getPos(size_t n,  Vec3R& pos) const { pos = particle_list_[n].p; }
    void getPosRad(size_t n,  Vec3R& pos, Real& rad) const {
      pos = particle_list_[n].p;
      rad = radius_scale_*particle_list_[n].r;
    }
    void getPosRadVel(size_t n,  Vec3R& pos, Real& rad, Vec3R& vel) const {
      pos = particle_list_[n].p;
      rad = radius_scale_*particle_list_[n].r;
      vel = velocity_scale_*particle_list_[n].v;
    }
    void getAtt(size_t n, openvdb::Index32& att) const { att = openvdb::Index32(n); }
};

int main(int argc, char *argv[]) {
  google::InitGoogleLogging(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  CHECK_NE(FLAGS_input, "") << "Must supply an input file name";
  CHECK_NE(FLAGS_output, "") << "Must supply an output file name";

  openvdb::initialize();

  const std::string input = FLAGS_input;
  const std::string output = FLAGS_output;

  /* Point cloud loading */
  CHECK(boost::filesystem::exists(input)) << "File " << input << " not found";
  FILE *fp = fopen(input.c_str(), "r");
  CHECK(fp != NULL) << "Error reading file " << input;

  // First read the header.
  Coord start, end, dims;
  int num_particles;
  CHECK(fscanf(fp, "%d %d %d\n", &start[0], &start[1], &start[2]) > 0);
  CHECK(fscanf(fp, "%d %d %d\n", &end[0], &end[1], &end[2]) > 0);
  dims = end - start + Coord(1);
  CHECK(fscanf(fp, "%d\n", &num_particles) > 0);

  // Now read particle data and add particles.
  ParticlesList pl(1, 1);
  for (int i = 0; i < num_particles; ++i) {
    Vec3R pos, vel;
    int p;
    CHECK(fscanf(fp, "%d : %lf %lf %lf, %lf %lf %lf\n", &p,
                 &pos[0], &pos[1], &pos[2], &vel[0], &vel[1], &vel[2]) > 0);
    pl.AddParticle(pos, vel);
  }
  fclose(fp);

  LOG(INFO) << "Loaded " << num_particles << " particles";
  CHECK_EQ(int(pl.size()), num_particles);

  /* OpenVDB processing */

  const float voxelSize = 1.f;
  openvdb::FloatGrid::Ptr ls = openvdb::createLevelSet<openvdb::FloatGrid>(voxelSize);

  openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid> raster(*ls);
  //raster.rasterizeTrails(pl, 0.75);
  raster.rasterizeSpheres(pl);

  openvdb::io::File unfiltered_output(output+"_unfiltered.vdb");
  openvdb::GridPtrVec unfiltered;
  unfiltered.push_back(ls);
  unfiltered_output.write(unfiltered);
  unfiltered_output.close();

  openvdb::tools::LevelSetFilter<openvdb::FloatGrid> filterer(*ls);
  filterer.mean();

  openvdb::io::File filtered_output(output+"_filtered.vdb");
  openvdb::GridPtrVec filtered;
  filtered.push_back(ls);
  filtered_output.write(unfiltered);
  filtered_output.close();

  /*
  std::vector<openvdb::Vec3s> points;
  std::vector<openvdb::Vec4I> quads;
  std::vector<openvdb::Vec3I> triangles;

  openvdb::tools::volumeToMesh(*ls, points, triangles, quads, 0.005, 0);

  std::cout << "points " << points.size() << std::endl;
  std::cout << "triangles " << triangles.size() << std::endl;
  std::cout << "quads " << quads.size() << std::endl;
  */

  return 0;
}

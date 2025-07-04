#pragma once

#include <Kin/kin.h>
#include <KOMO/objective.h>
#include <Optim/MathematicalProgram.h>
#include "fcl/broadphase/broadphase_dynamic_AABB_tree.h"
#include "fcl/broadphase/broadphase_dynamic_AABB_tree_array.h"
#include "fcl/broadphase/default_broadphase_callbacks.h"

#include <unordered_map>

#include "Animation.h"
typedef fcl::CollisionObject<float> CollObject;
typedef fcl::CollisionGeometry<float> CollGeom;

typedef fcl::Vector3<float> Vec3f;
typedef fcl::Quaternionf Quaternionf;
typedef fcl::BroadPhaseCollisionManager<float> BroadPhaseCollisionManager;
typedef fcl::DynamicAABBTreeCollisionManager<float> DynamicAABBTreeCollisionManager;
typedef fcl::CollisionRequest<float> CollisionRequest;
typedef fcl::CollisionResult<float> CollisionResult;
typedef fcl::DistanceRequest<float> DistanceRequest;
typedef fcl::DistanceResult<float> DistanceResult;
typedef fcl::Box<float> Box;
typedef fcl::Sphere<float> Sphere;
typedef fcl::Capsule<float> Capsule;
typedef fcl::Cylinder<float> Cylinder;


struct ConfigurationProblem;

struct QueryResult{
  // TODO: remove this again
  rai::Array<bool> computed;

  //goal features
  arr goal_y, goal_J;

  //collision features
  uintA collisions;
  arr coll_y, coll_J;
  arr normal_y, normal_J;
  arr side_J;

  bool isGoal=true;
  bool isFeasible=true;

  //optional a 3D coordinate for display
  arr disp3d;

  void getViolatedContacts(arr& y, arr& J, double margin=0.);
  arr getSideStep();
  arr getForwardStep();
  arr getBackwardStep(double relativeStepLength=1.1, double margin=0., const arr& nullRef=NoArr); //1.1 means overstepping the constraint IK a bit

  void write(ostream& os) const;
  void writeDetails(ostream& os, const ConfigurationProblem& P, double margin=0.) const;
};
stdOutPipe(QueryResult)

struct ConfigurationProblem {
  double PENETRATION_TOLERANCE = 0.;
  // const double PENETRATION_TOLERANCE = 1e-3;
  
  rai::Configuration C;
  arr q0, limits, max_step;
  rai::Array<shared_ptr<Objective>> objectives;

  int display=0;
  uint evals=0;
  bool computeCollisions;

  // do we want anything besides a boolean result?
  bool computeCollisionOnly = true;

  // ignore collision-pairs
  bool activeOnly = true;
  bool isActuated(const rai::Frame *f);
  std::unordered_map<uint, bool> actuated;

  ConfigurationProblem(const rai::Configuration& _C, bool _computeCollisions=true);

  shared_ptr<Objective> addObjective(const FeatureSymbol& feat, const StringA& frames, ObjectiveType type, const arr& scale=NoArr, const arr& target=NoArr);

  shared_ptr<QueryResult> queryUsingSplitFcl(const arr& x, const bool robot, const bool obs, const bool env);

  shared_ptr<QueryResult> query(const arr& x, const bool setJoints=true, const bool compProxies=true);
  shared_ptr<QueryResult> queryUsingFramestate(const arr &X);

  bool checkConnection(const arr &start, const arr &end, const uint disc=3, const bool binary=false);
  arr sample(const arr &start={}, const arr &goal={}, const double c_max=0, const double c_min=0);
};

struct TimedConfigurationProblem : ConfigurationProblem{
  rai::Animation A;
  fcl::DynamicAABBTreeCollisionManager<float> safe_interval_collision_manager;
  fcl::DynamicAABBTreeCollisionManager<float> static_collision_manager;
  std::vector<fcl::CollisionObject<float>*> collision_objects;
  std::vector<fcl::CollisionObject<float>*> static_collision_objects;
  std::vector<std::pair<rai::Frame*,CollObject*>> robot_link_objects;
  std::vector<double*> time_marks;
  double min_time;
  double max_time;
  bool collision_manager_was_initialised = false;
  void init_safe_interval_collisison_check(const arr &start, const double& t_start, const double& t_max);
  std::pair<rai::String,CollObject*> create_obstacle_fcl(const rai::Frame* frame);
  std::vector<std::pair<rai::String,CollObject*>> fcl_obstacles(const rai::Configuration& C, const bool only_static);
  std::vector<std::pair<rai::Frame*,CollObject*>> fcl_robots(const rai::Configuration& C);
  static bool BroadphaseCallback(CollObject* o1, CollObject* o2, void* cdata_);
  static bool staticCallback(CollObject* o1, CollObject* o2, void* cdata_);
  std::vector<double> collision_moments;
  TimedConfigurationProblem(const rai::Configuration& _C, const rai::Animation &_A);
  ConfigurationProblem getConfigurationProblemAtTime(const double t);

  ptr<QueryResult> query(const arr& x, const double t, const double tMin=0);
  ptr<QueryResult> query(const arr& x, const std::vector<double> times, const double tMin=0);
  std::vector<std::pair<double,double>> get_safe_intervals(const arr& x);
  bool checkEdge(const arr& x0, const double t0, const arr &x1, const double t1, const uint discretization=3);
  bool checkEdgeStatic(const arr& x0, const arr &x1);
  arr steer_to_obstacle(const arr& x0, const arr &x1);
  arr sample(const arr &start={}, const arr &goal={}, const double c_max=0, const double c_min=0);

  void get_animated_ids();
  std::set<int> animated_ids;
  std::set<uint> actuated_ids;
  bool static_collided = false;
  int robot_base_links_offset = 0;
};

struct GoalStateProblem : MathematicalProgram {
  ConfigurationProblem &P;
  double scaleCollisions = 1e1;
  int nEq=-1, nIneq=-1;

  GoalStateProblem(ConfigurationProblem& _P) : P(_P) {}

  virtual uint getDimension(){ return P.q0.N; }
  virtual void getFeatureTypes(ObjectiveTypeA& tt);
  virtual void evaluate(arr& phi, arr& J, const arr& x);
};

bool makePoseFeasible(arr& x, ConfigurationProblem& P, double IKstepSize=1.1, double maxQStepSize=.1, uint trials=3);



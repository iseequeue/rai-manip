#include <Kin/proxy.h>
#include <Optim/constrained.h>
#include <fcl/fcl.h>
#include <fcl/config.h>
// #include <rapidjson/document.h>
// #include <rapidjson/writer.h>
// #include <rapidjson/stringbuffer.h>
// #include <rapidjson/prettywriter.h>
// #include <fstream>
// #include <filesystem>

#include "ConfigurationProblem.h"

#include <Geo/fclInterface.h>
typedef fcl::CollisionObject<float> CollObject;
typedef fcl::CollisionGeometry<float> CollGeom;

typedef fcl::Vector3<float> Vec3f;
typedef fcl::Quaternionf Quaternionf;
typedef fcl::BroadPhaseCollisionManager<float> BroadPhaseCollisionManager;
typedef fcl::DynamicAABBTreeCollisionManager<float> DynamicAABBTreeCollisionManager;
typedef fcl::NaiveCollisionManager<float> NaiveCollisionManager;
typedef fcl::CollisionRequest<float> CollisionRequest;
typedef fcl::CollisionResult<float> CollisionResult;
typedef fcl::DistanceRequest<float> DistanceRequest;
typedef fcl::DistanceResult<float> DistanceResult;
typedef fcl::Box<float> Box;
typedef fcl::Sphere<float> Sphere;
typedef fcl::Capsule<float> Capsule;
typedef fcl::Cylinder<float> Cylinder;

ConfigurationProblem::ConfigurationProblem(const rai::Configuration& _C, bool _computeCollisions)
  : C(_C), computeCollisions(_computeCollisions) {

  q0 = C.getJointState();
  limits = C.getLimits();
  max_step = zeros(limits.d0);

  for(rai::Joint *j: C.activeJoints) {
    uint i=j->qIndex;
    uint d=j->qDim();
    if(d){
      switch(j->type) {
        case rai::JT_transXY:
        case rai::JT_transXYPhi:
        case rai::JT_free:
        case rai::JT_transX:
        case rai::JT_transZ:
        case rai::JT_trans3:
        case rai::JT_quatBall: 
        case rai::JT_hingeX: 
        case rai::JT_hingeY: 
        case rai::JT_hingeZ: 
        case rai::JT_rigid: {
          for(uint k=0; k<d; k++) max_step(i+k) = 1.;
        } break;

        default: NIY
      };
    }
  }
}

shared_ptr<Objective> ConfigurationProblem::addObjective(const FeatureSymbol& feat, const StringA& frames, ObjectiveType type, const arr& scale, const arr& target){
  shared_ptr<Feature> f = symbols2feature(feat, frames, C, scale, target, 0);

  shared_ptr<Objective> ob = make_shared<Objective>(f, type, f->shortTag(C));

  objectives.append(ob);
  return ob;
}

bool ConfigurationProblem::isActuated(const rai::Frame *f){
  if(actuated.count(f->ID) > 0){
    return actuated[f->ID];
  }

  if (f->joint && f->joint->active) {
    actuated[f->ID] = true;
    return true;
  }
  while(f->parent){
    if (f->joint && f->joint->active) {
      actuated[f->ID] = true;
      return true;
    }
    f = f->parent;
  }

  return false;
}

shared_ptr<QueryResult> ConfigurationProblem::queryUsingSplitFcl(const arr& x, const bool robot, const bool obs, const bool env){
  bool limitsRespected = true;
  constexpr double tol = 0.2;
  if(limits.N){
    for(uint i=0;i<x.N;i++){
      if(limits(i,1)>limits(i,0) && 
          (x.elem(i)<limits(i,0)*(1+tol) || x.elem(i)>limits(i,1)*(1+tol))){
        //LOG(-1) <<"QUERY OUT OF LIMIT: joint " <<i <<": " <<x.elem(i) <<' ' <<limits[i];
        limitsRespected = false;
        break;
      }
    }
  }

  shared_ptr<QueryResult> qr = make_shared<QueryResult>();

  if (!limitsRespected){
    qr->isFeasible = false;
    return qr;
  }

  C.setJointState(x);
  C.collideSplitFcl(robot, obs, env);
  evals++;

  if (C.proxies.N > 0){
    qr->isFeasible = false; return qr;
  }

  //C.watch();


  //collision features
  uint N = C.proxies.N;
  qr->collisions.resize(N, 2);
  qr->coll_y.resize(N, 1);
  qr->coll_J.resize(N, 1, x.N);
  qr->normal_y.resize(N, 3);
  qr->normal_J.resize(N, 3, x.N);
  qr->side_J.resize(N, 3, x.N);

  qr->computed.resize(N);
  qr->computed = false;

  // to make sure that we do not trigger something on accident
  qr->coll_y = 1;

  uint i=0;
  for(rai::Proxy& p:C.proxies){
    bool hasActive = true;
    if (activeOnly){
      hasActive = isActuated(p.a) | isActuated(p.b);

      if (!hasActive) {++i; continue;}
    }
    p.ensure_coll();

    qr->computed(i) = true;
    qr->collisions[i] =  TUP(p.a->ID, p.b->ID);
    arr Jp1, Jp2, Jx1, Jx2;
    if (!computeCollisionOnly){
      C.jacobian_pos(Jp1, C(p.a->ID), p.collision->p1);
      C.jacobian_pos(Jp2, C(p.b->ID), p.collision->p2);
      C.jacobian_angular(Jx1, C(p.a->ID));
      C.jacobian_angular(Jx2, C(p.b->ID));
    }
    if (computeCollisionOnly){ p.collision->kinDistance(qr->coll_y[i](), NoArr, NoArr, NoArr);}
    else{p.collision->kinDistance(qr->coll_y[i](), qr->coll_J[i](), Jp1, Jp2);}
    if (!computeCollisionOnly) {p.collision->kinNormal(qr->normal_y[i](), qr->normal_J[i](), Jp1, Jp2, Jx1, Jx2);}

    //std::cout << qr->coll_y[i](0) << std::endl;
    //if (computeCollisionOnly && qr->coll_y[i](0)<=-PENETRATION_TOLERANCE) {qr->isFeasible = false; return qr;}

    if (!computeCollisionOnly){
      arr a, b, Ja, Jb;
      C.kinematicsPos(a, Ja, C(p.a->ID));
      C.kinematicsPos(b, Jb, C(p.b->ID));
#if 0
      arr z = a-b;
      z /= length(z);
#else
      arr z = qr->normal_y[i];
#endif
      qr->side_J[i] = (eye(3) - (z^z))* (Ja - Jb);
    }

    i++;
  }
  CHECK_EQ(i, N, "");
  qr->coll_J.reshape(qr->coll_y.N, x.N);

  //is feasible?
  qr->isFeasible = (!qr->coll_y.N || min(qr->coll_y)>=-PENETRATION_TOLERANCE);
  if (!limitsRespected){
    qr->isFeasible = false;
  }

  if (computeCollisionOnly){
    return qr;
  }

  //goal features
  N=0;
  for(shared_ptr<Objective>& ob : objectives) N += ob->feat->__dim_phi(C);
  qr->goal_y.resize(N);
  qr->goal_J.resize(N, x.N);

  i=0;
  arr z, Jz;
  for(shared_ptr<Objective>& ob : objectives){
    ob->feat->__phi(z, Jz, C);
    for(uint j=0;j<z.N;j++){
      qr->goal_y(i+j) = z(j);
      qr->goal_J[i+j] = Jz[j];
    }
    i += z.N;
  }
  CHECK_EQ(i, N, "");

  //is goal?
  qr->isGoal= (absMax(qr->goal_y)<1e-2);

  //display (link of last joint)
  qr->disp3d = C.activeJoints.last()->frame->getPosition();

  if(display) C.watch(display>1, STRING("ConfigurationProblem query:\n" <<*qr));

  return qr;

}

shared_ptr<QueryResult> ConfigurationProblem::queryUsingFramestate(const arr &X){
  shared_ptr<QueryResult> qr = make_shared<QueryResult>();

  if(computeCollisions){
    C.stepFclUsingFramestate(X);
    for(rai::Proxy& p:C.proxies) {
      p.ensure_coll();
    }
  }
  evals++;

  //collision features
  const uint N = C.proxies.N;
  const uint xn = C.getJointState().N;

  qr->collisions.resize(N, 2);
  qr->coll_y.resize(N, 1);
  qr->coll_J.resize(N, 1, xn);
  qr->normal_y.resize(N, 3);
  qr->normal_J.resize(N, 3, xn);
  qr->side_J.resize(N, 3, xn);

  qr->computed.resize(N);
  qr->computed = false;

  // to make sure that we do not trigger something on accident
  qr->coll_y = 1;

  uint i=0;
  for(const rai::Proxy& p:C.proxies){
    bool hasActive = true;
    if (activeOnly){
      hasActive = isActuated(p.a) | isActuated(p.b);
      // std::cout<<p.a->name<< " "<<p.b->name<<std::endl;

      if (!hasActive) {++i; continue;}
    }

    qr->computed(i) = true;
    qr->collisions[i] =  TUP(p.a->ID, p.b->ID);
    arr Jp1, Jp2, Jx1, Jx2;
    if (!computeCollisionOnly){
      C.jacobian_pos(Jp1, C(p.a->ID), p.collision->p1);
      C.jacobian_pos(Jp2, C(p.b->ID), p.collision->p2);
      C.jacobian_angular(Jx1, C(p.a->ID));
      C.jacobian_angular(Jx2, C(p.b->ID));
    }
    if (computeCollisionOnly){ p.collision->kinDistance(qr->coll_y[i](), NoArr, NoArr, NoArr);}
    else{p.collision->kinDistance(qr->coll_y[i](), qr->coll_J[i](), Jp1, Jp2);}
    if (!computeCollisionOnly) {p.collision->kinNormal(qr->normal_y[i](), qr->normal_J[i](), Jp1, Jp2, Jx1, Jx2);}

    //std::cout << qr->coll_y[i](0) << std::endl;
    //if (computeCollisionOnly && qr->coll_y[i](0)<=-PENETRATION_TOLERANCE) {qr->isFeasible = false; return qr;}

    if (!computeCollisionOnly){
      arr a, b, Ja, Jb;
      C.kinematicsPos(a, Ja, C(p.a->ID));
      C.kinematicsPos(b, Jb, C(p.b->ID));
#if 0
      arr z = a-b;
      z /= length(z);
#else
      arr z = qr->normal_y[i];
#endif
      qr->side_J[i] = (eye(3) - (z^z))* (Ja - Jb);
    }

    i++;
  }
  CHECK_EQ(i, N, "");
  qr->coll_J.reshape(qr->coll_y.N, xn);

  //is feasible?
  qr->isFeasible = (!qr->coll_y.N || min(qr->coll_y)>=-PENETRATION_TOLERANCE);

  if (computeCollisionOnly){
    return qr;
  }

  //goal features
  uint Ngoal = 0;
  for(shared_ptr<Objective>& ob : objectives) Ngoal += ob->feat->__dim_phi(C);
  qr->goal_y.resize(Ngoal);
  qr->goal_J.resize(Ngoal, xn);

  i=0;
  arr z, Jz;
  for(shared_ptr<Objective>& ob : objectives){
    ob->feat->__phi(z, Jz, C);
    for(uint j=0;j<z.N;j++){
      qr->goal_y(i+j) = z(j);
      qr->goal_J[i+j] = Jz[j];
    }
    i += z.N;
  }
  CHECK_EQ(i, Ngoal, "");

  //is goal?
  qr->isGoal= (absMax(qr->goal_y)<1e-2);

  //display (link of last joint)
  qr->disp3d = C.activeJoints.last()->frame->getPosition();

  if(display) C.watch(display>1, STRING("ConfigurationProblem query:\n" <<*qr));

  return qr;
}

shared_ptr<QueryResult> ConfigurationProblem::query(const arr& x, const bool setJoints, const bool compProxies){
  bool limitsRespected = true;
  constexpr double tol = 0.2;
  limits = C.getLimits();
  if(limits.N){
    for(uint i=0;i<x.N;i++){
      if(limits(i,1)>limits(i,0) && 
          (x.elem(i)<limits(i,0)*(1+tol) || x.elem(i)>limits(i,1)*(1+tol))){
        LOG(-1) <<"QUERY OUT OF LIMIT: joint " <<i <<": " <<x.elem(i) <<' ' <<limits[i];
        std::cout<<"QUERY OUT OF LIMIT"<<std::endl;
        limitsRespected = false;
        break;
      }
    }
  }

  // C.fcl()->stopEarly = true;

  shared_ptr<QueryResult> qr = make_shared<QueryResult>();
  qr->disp3d = C.activeJoints.last()->frame->getPosition();

  if (!limitsRespected){
    qr->isFeasible = false;
    return qr;
  }

  if (setJoints){
    C.setJointState(x);
    qr->disp3d = C.activeJoints.last()->frame->getPosition();
  }
  if(computeCollisions && compProxies){
    //C.stepSwift();
    C.stepFcl();
    for(rai::Proxy& p:C.proxies) {
      p.ensure_coll();
    }
  }
  evals++;

  // C.watch(true);

  //collision features
  uint N = C.proxies.N;
  qr->collisions.resize(N, 2);
  qr->coll_y.resize(N, 1);
  qr->coll_J.resize(N, 1, x.N);
  qr->normal_y.resize(N, 3);
  qr->normal_J.resize(N, 3, x.N);
  qr->side_J.resize(N, 3, x.N);

  qr->computed.resize(N);
  qr->computed = false;

  // to make sure that we do not trigger something on accident
  qr->coll_y = 1;

  uint i=0;
  for(const rai::Proxy& p:C.proxies){
    bool hasActive = true;
    if (activeOnly){
      hasActive = isActuated(p.a) | isActuated(p.b);
      // std::cout<<p.a->name<< " "<<p.b->name<<std::endl;

      if (!hasActive) {++i; continue;}
    }

    qr->computed(i) = true;
    qr->collisions[i] =  TUP(p.a->ID, p.b->ID);
    arr Jp1, Jp2, Jx1, Jx2;
    if (!computeCollisionOnly){
      C.jacobian_pos(Jp1, C(p.a->ID), p.collision->p1);
      C.jacobian_pos(Jp2, C(p.b->ID), p.collision->p2);
      C.jacobian_angular(Jx1, C(p.a->ID));
      C.jacobian_angular(Jx2, C(p.b->ID));
    }
    if (computeCollisionOnly){ p.collision->kinDistance(qr->coll_y[i](), NoArr, NoArr, NoArr);}
    else{p.collision->kinDistance(qr->coll_y[i](), qr->coll_J[i](), Jp1, Jp2);}
    if (!computeCollisionOnly) {p.collision->kinNormal(qr->normal_y[i](), qr->normal_J[i](), Jp1, Jp2, Jx1, Jx2);}

    //std::cout << qr->coll_y[i](0) << std::endl;
    //if (computeCollisionOnly && qr->coll_y[i](0)<=-PENETRATION_TOLERANCE) {qr->isFeasible = false; return qr;}

    if (!computeCollisionOnly){
      arr a, b, Ja, Jb;
      C.kinematicsPos(a, Ja, C(p.a->ID));
      C.kinematicsPos(b, Jb, C(p.b->ID));
#if 0
      arr z = a-b;
      z /= length(z);
#else
      arr z = qr->normal_y[i];
#endif
      qr->side_J[i] = (eye(3) - (z^z))* (Ja - Jb);
    }

    i++;
  }
  CHECK_EQ(i, N, "");
  qr->coll_J.reshape(qr->coll_y.N, x.N);

  //is feasible?
  qr->isFeasible = (!qr->coll_y.N || min(qr->coll_y)>=-PENETRATION_TOLERANCE);
  // if (qr->coll_y.N > 0) std::cout << min(qr->coll_y) << std::endl;
  if (!limitsRespected){
    qr->isFeasible = false;
  }

  if (computeCollisionOnly){
    return qr;
  }

  //goal features
  N=0;
  for(shared_ptr<Objective>& ob : objectives) N += ob->feat->__dim_phi(C);
  qr->goal_y.resize(N);
  qr->goal_J.resize(N, x.N);

  i=0;
  arr z, Jz;
  for(shared_ptr<Objective>& ob : objectives){
    ob->feat->__phi(z, Jz, C);
    for(uint j=0;j<z.N;j++){
      qr->goal_y(i+j) = z(j);
      qr->goal_J[i+j] = Jz[j];
    }
    i += z.N;
  }
  CHECK_EQ(i, N, "");

  //is goal?
  qr->isGoal= (absMax(qr->goal_y)<1e-2);

  //display (link of last joint)
  qr->disp3d = C.activeJoints.last()->frame->getPosition();

  if(display) C.watch(display>1, STRING("ConfigurationProblem query:\n" <<*qr));

  return qr;
}

TimedConfigurationProblem::TimedConfigurationProblem(const rai::Configuration &_C, const rai::Animation &_A)
  : ConfigurationProblem(_C), A(_A){
}

ConfigurationProblem TimedConfigurationProblem::getConfigurationProblemAtTime(const double t){
  A.setToTime(C, t);
  return ConfigurationProblem(C);
}

ptr<QueryResult> TimedConfigurationProblem::query(const arr& x, const std::vector<double> times, const double tMin){
  ptr<QueryResult> res;
  for (const auto t: times){
    res = this->query(x, t, tMin);

    if (!res->isFeasible){
      return res;
    }
  }

  return res;
}


ptr<QueryResult> TimedConfigurationProblem::query(const arr& x, const double t, const double tMin){
  // this seems more efficient than getConfigurationProblemAtTime(t).query(x);

  A.setToTime(C, t, tMin);

  //for (auto f: A.prePlannedFrames) std::cout << f->name << std::endl;

  //std::cout << t << std::endl;
  //std::cout << "setting to time" << std::endl;
  //if (A.prePlannedFrames.N != 0 && t < 95) {C.watch(true);}

  // take the frame-state of the ones that are fixed
  arr poses(0, 7);
  if (A.prePlannedFrames.N != 0 && t < A.tPrePlanned){
    //std::cout << C["b0"]->getPose() << std::endl;

    for (auto f: A.prePlannedFrames) poses.append(C[f->name]->getPose());
    //poses = C.getFrameState(A.prePlannedFrames) * 1.;
    //std::cout << C["b0"]->getPose() << std::endl;
  }

    //std::cout << C["b0"]->ID << std::endl;

  //std::cout << poses << std::endl;

  //std::cout << "setting to x" << std::endl;
  C.setJointState(x);
    //std::cout << C["b0"]->getPose() << std::endl;
  //if (A.prePlannedFrames.N != 0 && t<95) {C.watch(true);}
  //std::cout << poses << std::endl;

  // set frame state of the ones that are fixed again
  if (A.prePlannedFrames.N != 0 && t < A.tPrePlanned){
    // C.watch(true);

    //C.setFrameState(poses, A.prePlannedFrames);
    uint cnt=0;
    for (const auto &f: A.prePlannedFrames) C[f->name]->setPose(poses[cnt++]);
    //std::cout << "setting back to time" << std::endl;
    //if (A.prePlannedFrames.N != 0 && t<95) {C.watch(true);}
    // std::cout << "corrected" << std::endl;
    // C.watch(true);
  }

  const auto res = ConfigurationProblem::query(x, false);

  //std::cout << res->isFeasible << std::endl;

  // if (!res->isFeasible) C.watch(false);
  //std::cout << "final check" << std::endl;
  // if (A.prePlannedFrames.N != 0) {std::cout << t << std::endl; C.watch(true);}

  return res;
}

void GoalStateProblem::getFeatureTypes(ObjectiveTypeA& tt){
  if(nEq==-1){
    NIY;
  }

  tt.resize(nEq+nIneq);
  tt({0,nEq-1}) = OT_eq;
  tt({nEq,-1}) = OT_ineq;

}

void GoalStateProblem::evaluate(arr& phi, arr& J, const arr& x){
  auto qr = P.query(x);
  phi = qr->goal_y;
  if(!!J) J = qr->goal_J;
  nEq = phi.N;

  if(scaleCollisions>0.){
    arr _y, _J;
    accumulateInequalities(_y, _J, -qr->coll_y, -qr->coll_J);
    phi.append( scaleCollisions*_y );
    if(!!J) J.append( scaleCollisions*_J );
    nIneq = _y.N;
  }else{
    nIneq = 0;
  }

  CHECK_EQ(phi.N, J.d0, "");
}

arr goalOptim(ConfigurationProblem& P){
  GoalStateProblem G(P);

  G.scaleCollisions=0.;
  arr x=P.q0;
  {
    arr dual;
    OptConstrained opt(x, dual, G);
    opt.run();
  }

  G.scaleCollisions=1e0;
  {
    arr dual;
    OptConstrained opt(x, dual, G);
    opt.run();
  }

  return x;
}

// TODO: This can probably be made more efficient for our use-case
// Produce the whole sequence in one go
double corput(uint n, uint base){
  double q=0;
  double bk=(double)1/base;

  while (n > 0) {
    q += (n % base)*bk;
    n /= base;
    bk /= base;
  }
  return q;
}

bool ConfigurationProblem::checkConnection(
    const arr &start,
    const arr &end,
    const uint disc,
    const bool binary){
  if (binary){
    for (uint i=1; i<disc; ++i){
      double ind = corput(i, 2);
      arr p = start + ind * (end-start);

      // TODO: change to check feasibility properly (with path constraints)
      if(!query(p)->isFeasible){
        return false;
      }

      //C.watch(false);
    }
  }
  else{
    for (uint i=1; i<disc-1; ++i){
      arr p = start + 1.0 * i / (disc-1) * (end-start);

      // TODO: change to check feasibility properly (with path constraints)
      if(!query(p)->isFeasible){
        return false;
      }
    }
  }
  return true;
}

void cp_localRndUniform(rai::Array<double>& a, double low=0., double high=1., bool add=false){
  if(!add) for(uint i=0; i<a.N; i++) a.p[i] =(double)rnd2_.uni(low, high);
  else     for(uint i=0; i<a.N; i++) a.p[i]+=(double)rnd2_.uni(low, high);
}

arr ConfigurationProblem::sample(const arr &start, const arr &goal, const double c_max, const double c_min){
  const arr limits = C.getLimits();
  const uint dim = limits.dim(0);
  arr sample(dim);

  if (start.d0 == 0){
    // sample uniformly between 0,1
    cp_localRndUniform(sample,0,1,false);

    // scale sample
    for (uint i=0; i<sample.d0; ++i){
      if(limits(i,1) > limits(i,0)){
        sample(i) = limits(i, 0) + sample(i) * (limits(i, 1) - limits(i, 0));
      }
      else {
        // default: [-5, 5]
        sample(i) = sample(i) * 10 - 5;
      }
    }
  }
  else{
    // generate ball: dropped coordinates method
    arr tmp(dim+2);
    rndGauss(tmp, 1, false);
    const double norm = length(tmp);
    for(uint i=0; i<sample.d0; ++i){ sample(i) = tmp(i) / norm;}
    // std::cout << sample << std::endl;

    // const double c_min = length(start - goal);
    const double val = sqrt(c_max*c_max - c_min*c_min)/2;

    // std::cout << c_min << std::endl;
    // std::cout << c_max << std::endl;

    arr r(dim);
    r(0) = c_max / 2;
    for (uint i=1;i<dim; ++i){ r(i) = val; }

    arr u = zeros(dim);
    u(0) = 1;
    arr C;
    rotationFromAtoB(C, (goal-start) / c_min, u);
    transpose(C);

    // std::cout << C << std::endl;

    arr L;
    L.setDiag(r);

    // std::cout << L << std::endl;

    const arr center = (start + goal) / 2.;
    sample = C * L * sample + center;

    // std::cout << sample << std::endl;

    // make sure that the sample is in the joint-limits
    // TODO: the version below is simply truncating, leading to non-uniform sampling
    // possible solutions: rejection sampling? (probably inefficient in lots of cases)
    for (uint i=0; i<sample.d0; ++i){
      if (abs(limits(i, 0) - limits(i, 1)) < 1e-6) {continue;} // no limit in this dim

      if(sample(i) > limits(i,1)){sample(i) = limits(i, 1);}
      if(sample(i) < limits(i,0)){sample(i) = limits(i, 0);}

      // the manual limits
      if(sample(i) < -5){sample(i) = -5;}
      if(sample(i) > 5){sample(i) = 5;}
    }
  }
  return sample;
}

std::vector<std::pair<rai::String,CollObject*>> TimedConfigurationProblem::fcl_obstacles(const rai::Configuration& C) {
  const arr X = C.getFrameState();
  
  std::vector<std::pair<rai::String,CollObject*>> result;
    std::vector<rai::Shape*> geometries;
    std::vector<uint> geometries_id;
    std::vector<rai::Frame*> frs;


    for(rai::Frame* f:C.frames) {
      if(isActuated(f)){
        // std::cout<<"Skipped actuated frame"<<std::endl;
        continue;
      }
      bool cant_collide = false;
      for (std::pair<rai::Frame*,CollObject*> &b:this->robot_link_objects){
        // std::cout<<"robot joint name in can't  collide"<<b.first->name<<std::endl;
        if(f->name.contains("table")){
            continue;
          }
        rai::Frame* c = C.getFrames(TUP(b.first->ID))(0);
        if (!(f->getShape().canCollideWith(c))) {
          // std::cout<<"robot joint name in can't collide"<<b.first->name<<" and "<<f->name<<std::endl;
          cant_collide = true;
          break;
        }

        // Skip if one of the objects is an object we move.
        // Otherwise some collision with the table will be ignored.
        // if (f->name.contains("obj") || c->name.contains("obj")){
        //   cant_collide = true;
        //   break;
        // }

        if (c == f->parent || f == c->parent) {
          if (f->name.contains("table") && c->name.contains("obj")){
            break;
          }
          cant_collide = true;
          break;
        }
        if (c == f->getUpwardLink() || f == c->getDownwardLink()) {
          if (f->name.contains("table") && c->name.contains("obj")){
            break;
          }
          cant_collide = true;
          break;
        }
    }
    if (cant_collide){
      // std::cout<<"name, that cant collide: "<<f->name<<std::endl;
      continue;
    }
      if(f->shape && f->shape->cont) {
        if(!f->shape->mesh().V.N) f->shape->createMeshes();

        frs.push_back(f);
        geometries.push_back(frs.back()->shape);
        geometries_id.push_back(frs.back()->ID);
      }
    }

    for(long int i=0; i<geometries.size(); i++) {
      rai::Shape* shape = geometries[i];
  
      if(shape) {
        std::shared_ptr<CollGeom> geom;
        if(shape->type()==rai::ST_capsule) {
          geom = make_shared<Capsule>(shape->size(-1), shape->size(-2));
        } else if(shape->type()==rai::ST_cylinder) {
          geom = make_shared<Cylinder>(shape->size(-1), shape->size(-2));
        } else if(shape->type()==rai::ST_sphere) {
          geom = make_shared<Sphere>(shape->size(-1));
        }else if(shape->type()==rai::ST_box){
          geom = make_shared<Box>(shape->size(0), shape->size(1), shape->size(2));
        } else {
          std::cout << shape->type() << std::endl;
          std::cout << shape->frame.name << std::endl;
          rai::Mesh& mesh = shape->mesh();

          mesh.computeNormals();
  
          auto verts = make_shared<std::vector<fcl::Vector3<float>>>(mesh.V.d0);
          auto faces = make_shared<std::vector<int>>(mesh.T.N);
          for(uint i=0; i<verts->size(); i++){(*verts)[i] = {(float)mesh.V(i, 0), (float)mesh.V(i, 1), (float)mesh.V(i, 2)};
          std::cout << mesh.V(i, 0) << " " << mesh.V(i, 1) << " " << mesh.V(i, 2) << std::endl;}
          for(uint i=0; i<faces->size(); i++){(*faces)[i] = mesh.T.elem(i) * 1.;
          std::cout << mesh.T.elem(i) << std::endl;}
          geom = make_shared<fcl::Convex<float>>(verts, mesh.T.d0, faces, true);

          const auto model = make_shared<fcl::Sphere<float>>(mesh.getRadius());
        }
        CollObject* obj = new CollObject(geom);
        obj->setTranslation(Vec3f(X(geometries_id[i], 0), X(geometries_id[i], 1), X(geometries_id[i], 2)));
        obj->setQuatRotation(Quaternionf(X(geometries_id[i], 3), X(geometries_id[i], 4), X(geometries_id[i], 5), X(geometries_id[i], 6)));
        obj->computeAABB();
        
        result.push_back(std::pair<rai::String,CollObject*>(frs[i]->name,obj));
      }
    }
  return result;
}

std::vector<std::pair<rai::Frame*,CollObject*>> TimedConfigurationProblem::fcl_robots(const rai::Configuration& C) {
  const arr X = C.getFrameState();
  
  std::vector<std::pair<rai::Frame*,CollObject*>> result;
    std::vector<rai::Shape*> geometries;
    std::vector<uint> geometries_id;
    std::vector<rai::Frame*> frs;


    for(rai::Frame* f:C.frames) {
      if(!isActuated(f)){
        // std::cout<<"Skipped unactuated frame"<<std::endl;

        continue;
      }
      // if (f->name.contains("obj")) {
      //   continue;
      // }

      if(f->shape && f->shape->cont) {
        // std::cout<<"robot joint name "<<f->name<<std::endl;
        
        if(!f->shape->mesh().V.N) f->shape->createMeshes();
        geometries.push_back(f->shape);
        geometries_id.push_back(f->ID);
        frs.push_back(f);
      }
    }

    for(long int i=0; i<geometries.size(); i++) {
      rai::Shape* shape = geometries[i];
  
      if(shape) {
        std::shared_ptr<CollGeom> geom;
        if(shape->type()==rai::ST_capsule) {
          geom = make_shared<Capsule>(shape->size(-1), shape->size(-2));
        } else if(shape->type()==rai::ST_cylinder) {
          geom = make_shared<Cylinder>(shape->size(-1), shape->size(-2));
        } else if(shape->type()==rai::ST_sphere) {
          geom = make_shared<Sphere>(shape->size(-1));
        }else if(shape->type()==rai::ST_box){
          geom = make_shared<Box>(shape->size(0), shape->size(1), shape->size(2));
        } else {
          std::cout << shape->type() << std::endl;
          std::cout << shape->frame.name << std::endl;
          rai::Mesh& mesh = shape->mesh();

          mesh.computeNormals();
  
          auto verts = make_shared<std::vector<fcl::Vector3<float>>>(mesh.V.d0);
          auto faces = make_shared<std::vector<int>>(mesh.T.N);
          for(uint i=0; i<verts->size(); i++){(*verts)[i] = {(float)mesh.V(i, 0), (float)mesh.V(i, 1), (float)mesh.V(i, 2)};
          std::cout << mesh.V(i, 0) << " " << mesh.V(i, 1) << " " << mesh.V(i, 2) << std::endl;}
          for(uint i=0; i<faces->size(); i++){(*faces)[i] = mesh.T.elem(i) * 1.;
          std::cout << mesh.T.elem(i) << std::endl;}
          geom = make_shared<fcl::Convex<float>>(verts, mesh.T.d0, faces, true);

          const auto model = make_shared<fcl::Sphere<float>>(mesh.getRadius());
        }
        CollObject* obj = new CollObject(geom);
        obj->setTranslation(Vec3f(X(geometries_id[i], 0), X(geometries_id[i], 1), X(geometries_id[i], 2)));
        obj->setQuatRotation(Quaternionf(X(geometries_id[i], 3), X(geometries_id[i], 4), X(geometries_id[i], 5), X(geometries_id[i], 6)));
        obj->computeAABB();
        
        result.push_back(std::pair<rai::Frame*,CollObject*>(frs[i],obj));
      }
    }
  return result;
}

std::pair<rai::String,CollObject*> TimedConfigurationProblem::create_obstacle_fcl(const rai::Frame* frame){
  
  rai::Shape* shape = frame->shape;

  assert(shape);

  std::shared_ptr<CollGeom> geom;
  if(shape->type()==rai::ST_capsule) {
    geom = make_shared<Capsule>(shape->size(-1), shape->size(-2));
  } else if(shape->type()==rai::ST_cylinder) {
    geom = make_shared<Cylinder>(shape->size(-1), shape->size(-2));
  } else if(shape->type()==rai::ST_sphere) {
    geom = make_shared<Sphere>(shape->size(-1));
  }else if(shape->type()==rai::ST_box){
    geom = make_shared<Box>(shape->size(0), shape->size(1), shape->size(2));
  } else {
    std::cout << shape->type() << std::endl;
    std::cout << shape->frame.name << std::endl;
    rai::Mesh& mesh = shape->mesh();

    mesh.computeNormals();

    auto verts = make_shared<std::vector<fcl::Vector3<float>>>(mesh.V.d0);
    auto faces = make_shared<std::vector<int>>(mesh.T.N);
    for(uint i=0; i<verts->size(); i++){(*verts)[i] = {(float)mesh.V(i, 0), (float)mesh.V(i, 1), (float)mesh.V(i, 2)};
    std::cout << mesh.V(i, 0) << " " << mesh.V(i, 1) << " " << mesh.V(i, 2) << std::endl;}
    for(uint i=0; i<faces->size(); i++){(*faces)[i] = mesh.T.elem(i) * 1.;
    std::cout << mesh.T.elem(i) << std::endl;}
    geom = make_shared<fcl::Convex<float>>(verts, mesh.T.d0, faces, true);

    const auto model = make_shared<fcl::Sphere<float>>(mesh.getRadius());
  }
  CollObject* obj = new CollObject(geom);
  obj->setTranslation(Vec3f(frame->X.pos.x, frame->X.pos.y, frame->X.pos.z));
  obj->setQuatRotation(Quaternionf(frame->X.rot.w, frame->X.rot.x, frame->X.rot.y, frame->X.rot.z));
  obj->computeAABB();
  
  return std::pair<rai::String,CollObject*>(frame->name,obj);
}

void TimedConfigurationProblem::init_safe_interval_collisison_check(const arr &start, const double& t_start, const double& t_max){
  //Hardcode: 1 fps due to animation realisation.

  // std::cout<<"Frames: "<<std::endl;
  // for (rai::Frame* fr:C.frames){
  //   if(fr->parent){
  // std::cout<<"Frame: "<<fr->name<<" parent: "<<fr->parent->name<< " parentLink: "<<fr->getUpwardLink()->name<< " isActuated: " <<isActuated(fr) <<std::endl;

  //   }
  //   else{
  // std::cout<<"Frame: "<<fr->name<<" parent: None"<<std::endl;
  //   }
  // }


  //destructor
  if(collision_manager_was_initialised){
    // std::cout<<"destructor"<<std::endl;
    for(int col_id =0; col_id<this->collision_objects.size();col_id++){
      delete (std::pair<double*,rai::String>*) (this->collision_objects[col_id]->getUserData());
      delete this->collision_objects[col_id];
    }
    for(int col_id =0; col_id<this->robot_link_objects.size();col_id++){
      delete (std::pair<double*,rai::Frame*>*) (this->robot_link_objects[col_id].second->getUserData());
      delete this->robot_link_objects[col_id].second;
    }
    for(int time_id =0; time_id<this->time_marks.size();time_id++){
      delete this->time_marks[time_id];
    }

    
  }
  this->time_marks.clear();
  this->collision_objects.clear();
  this->robot_link_objects.clear();
  this->safe_interval_collision_manager.clear();

  C.setJointState(start);//set joint pos
   // get robot links collision objects
   this->robot_link_objects = this->fcl_robots(this->C);

  uint number_of_frames = this->max_time-this->min_time+1;

  // std::cout<<"number_of_frames "<<number_of_frames<<std::endl;


  // for(int i=A.d0-1; i>=0; --i){
  //   const auto &a = A(i);
  //   const double &start = a.start;

  //   if (t_max < start || start < 0 || a.X.d0 == 0) {
  //     LOG(-1) <<"Warning: start < 0";
  //     continue;
  //   }
  //   if (t_start > start + a.X.d0){
  //     continue;
  //   }

    
  //   for (int i=0; i< a.X.d0;i++){
  //     t_moment = start + i;
  //     for (int j=0; j< a.frameIDs.d0;j++)
  //     rai::Frame* curr_frame = this->C.getFrames(TUP(a.frameIDs(j)))(0);
  //     curr_frame.X = a.X(i)(j);

  //     this->time_marks.push_back(new double(t_moment));
  //     std::pair<rai::String,CollObject*> result = this->create_obstacle_fcl(curr_frame)
  //     this->collision_objects.push_back(result.second);
  //     this->collision_objects.back()->setUserData((void*)(new std::pair<double*,rai::String>(this->time_marks.back(),result.first)));
  //     this->collision_objects.back()->computeAABB();
  //   }
  // }
  

  
  for (int frame_number=0;frame_number<number_of_frames;frame_number++){
    ConfigurationProblem conf = this->getConfigurationProblemAtTime(this->min_time+frame_number);
    std::vector<std::pair<rai::String,CollObject*>> obstacles =  this->fcl_obstacles(conf.C);
    // for( auto obs:obstacles){
    //     std::cout<<"name, to be added to fcl "<<obs.first<<std::endl;
    //     // std::cout<<obs.second->getTranslation()<<std::endl;
        
    // }
    this->time_marks.push_back(new double(this->min_time+frame_number));
    // std::cout<<"obstacles "<<obstacles.size()<<std::endl;
    
    for (int col_id=0;col_id<obstacles.size();col_id++){  
      // if(!obstacles[col_id].second){
      //   // std::cout<<"obstacles NULLPTR DETECTED! "<<col_id<<std::endl;
      // }
      this->collision_objects.push_back(obstacles[col_id].second);
      this->collision_objects.back()->setUserData((void*)(new std::pair<double*,rai::String>(this->time_marks.back(),obstacles[col_id].first)));
      this->collision_objects.back()->computeAABB();
    }
  }
  for (int frame_number=0;frame_number<this->collision_objects.size();frame_number++){
    // if(!collision_objects[frame_number]){
    //   // std::cout<<"NULLPTR DETECTED! "<<number_of_frames<<std::endl;
    // }
  }
  this->safe_interval_collision_manager.registerObjects(this->collision_objects);
  this->safe_interval_collision_manager.setup();
  this->collision_manager_was_initialised = true;


 

  //to identify link in collision callback:
  this->time_marks.push_back(new double(-1));
  for(int robot_joint_id = 0; robot_joint_id < this->robot_link_objects.size(); robot_joint_id++){
    this->robot_link_objects[robot_joint_id].second->setUserData((void*)(new std::pair<double*,rai::Frame*>(this->time_marks.back(),this->robot_link_objects[robot_joint_id].first)));
  }

  // // Export FCL objects to JSON for Blender visualization using rapidjson
  // rapidjson::Document document;
  // document.SetObject();
  // rapidjson::Document::AllocatorType& allocator = document.GetAllocator();
  
  // // Create obstacle and robot arrays
  // rapidjson::Value obstaclesArray(rapidjson::kArrayType);
  // rapidjson::Value robotsArray(rapidjson::kArrayType);
  
  // // Export obstacles
  // for(size_t i = 0; i < this->collision_objects.size(); i++) {
  //   CollObject* obj = this->collision_objects[i];
  //   std::pair<double*, rai::String>* userData = static_cast<std::pair<double*, rai::String>*>(obj->getUserData());
  //   double time = *(userData->first);
  //   std::string name = userData->second.p;
    
  //   rapidjson::Value obstacle(rapidjson::kObjectType);
    
  //   // Add name
  //   rapidjson::Value nameValue;
  //   nameValue.SetString(name.c_str(), name.length(), allocator);
  //   obstacle.AddMember("name", nameValue, allocator);
    
  //   // Add time
  //   obstacle.AddMember("time", time, allocator);
    
  //   // Export translation
  //   Vec3f translation = obj->getTranslation();
  //   rapidjson::Value position(rapidjson::kArrayType);
  //   position.PushBack(translation[0], allocator);
  //   position.PushBack(translation[1], allocator);
  //   position.PushBack(translation[2], allocator);
  //   obstacle.AddMember("position", position, allocator);
    
  //   // Export rotation
  //   Quaternionf rotation = obj->getQuatRotation();
  //   rapidjson::Value rotationArray(rapidjson::kArrayType);
  //   rotationArray.PushBack(rotation.w(), allocator);
  //   rotationArray.PushBack(rotation.x(), allocator);
  //   rotationArray.PushBack(rotation.y(), allocator);
  //   rotationArray.PushBack(rotation.z(), allocator);
  //   obstacle.AddMember("rotation", rotationArray, allocator);
    
  //   // Export geometry type
  //   auto geom = obj->collisionGeometry();
  //   obstacle.AddMember("geometry_type", geom->getNodeType(), allocator);
    
  //   // Handle different geometry types
  //   rapidjson::Value geometry(rapidjson::kObjectType);
    
  //   if (geom->getNodeType() == fcl::GEOM_BOX) {
  //     const Box* box = static_cast<const Box*>(geom.get());
      
  //     rapidjson::Value typeValue;
  //     typeValue.SetString("box", allocator);
  //     geometry.AddMember("type", typeValue, allocator);
      
  //     rapidjson::Value sizeArray(rapidjson::kArrayType);
  //     sizeArray.PushBack(box->side[0], allocator);
  //     sizeArray.PushBack(box->side[1], allocator);
  //     sizeArray.PushBack(box->side[2], allocator);
  //     geometry.AddMember("size", sizeArray, allocator);
  //   }
  //   else if (geom->getNodeType() == fcl::GEOM_SPHERE) {
  //     const Sphere* sphere = static_cast<const Sphere*>(geom.get());
      
  //     rapidjson::Value typeValue;
  //     typeValue.SetString("sphere", allocator);
  //     geometry.AddMember("type", typeValue, allocator);
      
  //     geometry.AddMember("radius", sphere->radius, allocator);
  //   }
  //   else if (geom->getNodeType() == fcl::GEOM_CAPSULE) {
  //     const Capsule* capsule = static_cast<const Capsule*>(geom.get());
      
  //     rapidjson::Value typeValue;
  //     typeValue.SetString("capsule", allocator);
  //     geometry.AddMember("type", typeValue, allocator);
      
  //     geometry.AddMember("radius", capsule->radius, allocator);
  //     geometry.AddMember("height", capsule->lz, allocator);
  //   }
  //   else if (geom->getNodeType() == fcl::GEOM_CYLINDER) {
  //     const Cylinder* cylinder = static_cast<const Cylinder*>(geom.get());
      
  //     rapidjson::Value typeValue;
  //     typeValue.SetString("cylinder", allocator);
  //     geometry.AddMember("type", typeValue, allocator);
      
  //     geometry.AddMember("radius", cylinder->radius, allocator);
  //     geometry.AddMember("height", cylinder->lz, allocator);
  //   }
  //   else if (geom->getNodeType() == fcl::GEOM_CONVEX) {
  //     rapidjson::Value typeValue;
  //     typeValue.SetString("convex", allocator);
  //     geometry.AddMember("type", typeValue, allocator);
  //     // Would need to extract vertices and faces for complete export
  //   }
    
  //   obstacle.AddMember("geometry", geometry, allocator);
  //   obstaclesArray.PushBack(obstacle, allocator);
  // }
  
  // // Export robot links
  // for(size_t i = 0; i < this->robot_link_objects.size(); i++) {
  //   CollObject* obj = this->robot_link_objects[i].second;
  //   rai::Frame* frame = this->robot_link_objects[i].first;
    
  //   rapidjson::Value robot(rapidjson::kObjectType);
    
  //   // Add name
  //   rapidjson::Value nameValue;
  //   nameValue.SetString(frame->name.p, frame->name.N, allocator);
  //   robot.AddMember("name", nameValue, allocator);
    
  //   // Add time (-1 for robot links as they aren't time-dependent)
  //   robot.AddMember("time", -1, allocator);
    
  //   // Export position
  //   Vec3f translation = obj->getTranslation();
  //   rapidjson::Value position(rapidjson::kArrayType);
  //   position.PushBack(translation[0], allocator);
  //   position.PushBack(translation[1], allocator);
  //   position.PushBack(translation[2], allocator);
  //   robot.AddMember("position", position, allocator);
    
  //   // Export rotation
  //   Quaternionf rotation = obj->getQuatRotation();
  //   rapidjson::Value rotationArray(rapidjson::kArrayType);
  //   rotationArray.PushBack(rotation.w(), allocator);
  //   rotationArray.PushBack(rotation.x(), allocator);
  //   rotationArray.PushBack(rotation.y(), allocator);
  //   rotationArray.PushBack(rotation.z(), allocator);
  //   robot.AddMember("rotation", rotationArray, allocator);
    
  //   // Export geometry type
  //   auto geom = obj->collisionGeometry();
  //   robot.AddMember("geometry_type", geom->getNodeType(), allocator);
    
  //   // Handle different geometry types (same as for obstacles)
  //   rapidjson::Value geometry(rapidjson::kObjectType);
    
  //   if (geom->getNodeType() == fcl::GEOM_BOX) {
  //     const Box* box = static_cast<const Box*>(geom.get());
      
  //     rapidjson::Value typeValue;
  //     typeValue.SetString("box", allocator);
  //     geometry.AddMember("type", typeValue, allocator);
      
  //     rapidjson::Value sizeArray(rapidjson::kArrayType);
  //     sizeArray.PushBack(box->side[0], allocator);
  //     sizeArray.PushBack(box->side[1], allocator);
  //     sizeArray.PushBack(box->side[2], allocator);
  //     geometry.AddMember("size", sizeArray, allocator);
  //   }
  //   else if (geom->getNodeType() == fcl::GEOM_SPHERE) {
  //     const Sphere* sphere = static_cast<const Sphere*>(geom.get());
      
  //     rapidjson::Value typeValue;
  //     typeValue.SetString("sphere", allocator);
  //     geometry.AddMember("type", typeValue, allocator);
      
  //     geometry.AddMember("radius", sphere->radius, allocator);
  //   }
  //   else if (geom->getNodeType() == fcl::GEOM_CAPSULE) {
  //     const Capsule* capsule = static_cast<const Capsule*>(geom.get());
      
  //     rapidjson::Value typeValue;
  //     typeValue.SetString("capsule", allocator);
  //     geometry.AddMember("type", typeValue, allocator);
      
  //     geometry.AddMember("radius", capsule->radius, allocator);
  //     geometry.AddMember("height", capsule->lz, allocator);
  //   }
  //   else if (geom->getNodeType() == fcl::GEOM_CYLINDER) {
  //     const Cylinder* cylinder = static_cast<const Cylinder*>(geom.get());
      
  //     rapidjson::Value typeValue;
  //     typeValue.SetString("cylinder", allocator);
  //     geometry.AddMember("type", typeValue, allocator);
      
  //     geometry.AddMember("radius", cylinder->radius, allocator);
  //     geometry.AddMember("height", cylinder->lz, allocator);
  //   }
  //   else if (geom->getNodeType() == fcl::GEOM_CONVEX) {
  //     rapidjson::Value typeValue;
  //     typeValue.SetString("convex", allocator);
  //     geometry.AddMember("type", typeValue, allocator);
  //     // Would need to extract vertices and faces for complete export
  //   }
    
  //   robot.AddMember("geometry", geometry, allocator);
  //   robotsArray.PushBack(robot, allocator);
  // }
  
  // // Add arrays to document
  // document.AddMember("obstacles", obstaclesArray, allocator);
  // document.AddMember("robots", robotsArray, allocator);
  
  // // Create the output directory if it doesn't exist
  // std::filesystem::create_directories("fcl_export");
  
  // // Write to file with pretty printing
  // rapidjson::StringBuffer buffer;
  // rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
  // document.Accept(writer);
  
  // std::ofstream file("fcl_export/scene.json");
  // file << buffer.GetString();
  // file.close();
  
  // std::cout << "Exported FCL scene to fcl_export/scene.json" << std::endl;
}

std::vector<std::pair<double,double>> TimedConfigurationProblem::get_safe_intervals(const arr& x){
  std::vector<std::pair<double,double>> results;
  
  C.setJointState(x);//set joint pos
  
  // // Export current scene state to JSON - using the same format as init_safe_interval_collisison_check
  // rapidjson::Document document;
  // document.SetObject();
  // rapidjson::Document::AllocatorType& allocator = document.GetAllocator();
  
  // // Create obstacle and robot arrays (same structure as in init_safe_interval_collisison_check)
  // rapidjson::Value obstaclesArray(rapidjson::kArrayType);
  // rapidjson::Value robotsArray(rapidjson::kArrayType);
  
  // // Get the obstacles (non-actuated frames with shapes)
  // std::vector<std::pair<rai::String,CollObject*>> obstacles = this->fcl_obstacles(C);
  
  // // Export obstacles
  // for(auto& obs : obstacles) {
  //   CollObject* obj = obs.second;
  //   std::string name = obs.first.p;
    
  //   rapidjson::Value obstacle(rapidjson::kObjectType);
    
  //   // Add name
  //   rapidjson::Value nameValue;
  //   nameValue.SetString(name.c_str(), name.length(), allocator);
  //   obstacle.AddMember("name", nameValue, allocator);
    
  //   // Add time (current state)
  //   obstacle.AddMember("time", 0, allocator);
    
  //   // Export translation
  //   Vec3f translation = obj->getTranslation();
  //   rapidjson::Value position(rapidjson::kArrayType);
  //   position.PushBack(translation[0], allocator);
  //   position.PushBack(translation[1], allocator);
  //   position.PushBack(translation[2], allocator);
  //   obstacle.AddMember("position", position, allocator);
    
  //   // Export rotation
  //   Quaternionf rotation = obj->getQuatRotation();
  //   rapidjson::Value rotationArray(rapidjson::kArrayType);
  //   rotationArray.PushBack(rotation.w(), allocator);
  //   rotationArray.PushBack(rotation.x(), allocator);
  //   rotationArray.PushBack(rotation.y(), allocator);
  //   rotationArray.PushBack(rotation.z(), allocator);
  //   obstacle.AddMember("rotation", rotationArray, allocator);
    
  //   // Export geometry type
  //   auto geom = obj->collisionGeometry();
  //   obstacle.AddMember("geometry_type", geom->getNodeType(), allocator);
    
  //   // Handle different geometry types
  //   rapidjson::Value geometry(rapidjson::kObjectType);
    
  //   if (geom->getNodeType() == fcl::GEOM_BOX) {
  //     const Box* box = static_cast<const Box*>(geom.get());
      
  //     rapidjson::Value typeValue;
  //     typeValue.SetString("box", allocator);
  //     geometry.AddMember("type", typeValue, allocator);
      
  //     rapidjson::Value sizeArray(rapidjson::kArrayType);
  //     sizeArray.PushBack(box->side[0], allocator);
  //     sizeArray.PushBack(box->side[1], allocator);
  //     sizeArray.PushBack(box->side[2], allocator);
  //     geometry.AddMember("size", sizeArray, allocator);
  //   }
  //   else if (geom->getNodeType() == fcl::GEOM_SPHERE) {
  //     const Sphere* sphere = static_cast<const Sphere*>(geom.get());
      
  //     rapidjson::Value typeValue;
  //     typeValue.SetString("sphere", allocator);
  //     geometry.AddMember("type", typeValue, allocator);
      
  //     geometry.AddMember("radius", sphere->radius, allocator);
  //   }
  //   else if (geom->getNodeType() == fcl::GEOM_CAPSULE) {
  //     const Capsule* capsule = static_cast<const Capsule*>(geom.get());
      
  //     rapidjson::Value typeValue;
  //     typeValue.SetString("capsule", allocator);
  //     geometry.AddMember("type", typeValue, allocator);
      
  //     geometry.AddMember("radius", capsule->radius, allocator);
  //     geometry.AddMember("height", capsule->lz, allocator);
  //   }
  //   else if (geom->getNodeType() == fcl::GEOM_CYLINDER) {
  //     const Cylinder* cylinder = static_cast<const Cylinder*>(geom.get());
      
  //     rapidjson::Value typeValue;
  //     typeValue.SetString("cylinder", allocator);
  //     geometry.AddMember("type", typeValue, allocator);
      
  //     geometry.AddMember("radius", cylinder->radius, allocator);
  //     geometry.AddMember("height", cylinder->lz, allocator);
  //   }
  //   else if (geom->getNodeType() == fcl::GEOM_CONVEX) {
  //     rapidjson::Value typeValue;
  //     typeValue.SetString("convex", allocator);
  //     geometry.AddMember("type", typeValue, allocator);
  //     // Would need to extract vertices and faces for complete export
  //   }
    
  //   obstacle.AddMember("geometry", geometry, allocator);
  //   obstaclesArray.PushBack(obstacle, allocator);
    
  //   // Clean up the temporary obstacle object
  //   delete obj;
  // }
  
  // // Export robot links (actuated frames)
  // std::vector<std::pair<rai::Frame*,CollObject*>> robotLinks = this->fcl_robots(C);
  
  // for(auto& robot_pair : robotLinks) {
  //   rai::Frame* frame = robot_pair.first;
  //   CollObject* obj = robot_pair.second;
    
  //   rapidjson::Value robot(rapidjson::kObjectType);
    
  //   // Add name
  //   rapidjson::Value nameValue;
  //   nameValue.SetString(frame->name.p, frame->name.N, allocator);
  //   robot.AddMember("name", nameValue, allocator);
    
  //   // Add time (-1 for robot links as they aren't time-dependent)
  //   robot.AddMember("time", -1, allocator);
    
  //   // Export position
  //   Vec3f translation = obj->getTranslation();
  //   rapidjson::Value position(rapidjson::kArrayType);
  //   position.PushBack(translation[0], allocator);
  //   position.PushBack(translation[1], allocator);
  //   position.PushBack(translation[2], allocator);
  //   robot.AddMember("position", position, allocator);
    
  //   // Export rotation
  //   Quaternionf rotation = obj->getQuatRotation();
  //   rapidjson::Value rotationArray(rapidjson::kArrayType);
  //   rotationArray.PushBack(rotation.w(), allocator);
  //   rotationArray.PushBack(rotation.x(), allocator);
  //   rotationArray.PushBack(rotation.y(), allocator);
  //   rotationArray.PushBack(rotation.z(), allocator);
  //   robot.AddMember("rotation", rotationArray, allocator);
    
  //   // Export joint data (additional data not in init_safe_interval_collisison_check)
  //   if(frame->joint && frame->joint->active) {
  //     rapidjson::Value jointData(rapidjson::kObjectType);
      
  //     // Joint index
  //     jointData.AddMember("qIndex", frame->joint->qIndex, allocator);
      
  //     // Joint type
  //     const char* jointTypeName = "";
  //     switch(frame->joint->type) {
  //       case rai::JT_hingeX: jointTypeName = "hingeX"; break;
  //       case rai::JT_hingeY: jointTypeName = "hingeY"; break;
  //       case rai::JT_hingeZ: jointTypeName = "hingeZ"; break;
  //       case rai::JT_transX: jointTypeName = "transX"; break;
  //       case rai::JT_transY: jointTypeName = "transY"; break;
  //       case rai::JT_transZ: jointTypeName = "transZ"; break;
  //       case rai::JT_transXY: jointTypeName = "transXY"; break;
  //       case rai::JT_trans3: jointTypeName = "trans3"; break;
  //       case rai::JT_transXYPhi: jointTypeName = "transXYPhi"; break;
  //       case rai::JT_quatBall: jointTypeName = "quatBall"; break;
  //       case rai::JT_free: jointTypeName = "free"; break;
  //       case rai::JT_rigid: jointTypeName = "rigid"; break;
  //       default: jointTypeName = "unknown";
  //     }
      
  //     rapidjson::Value jointTypeValue;
  //     jointTypeValue.SetString(jointTypeName, allocator);
  //     jointData.AddMember("type", jointTypeValue, allocator);
      
  //     // Joint angle if available
  //     if(frame->joint->qDim() > 0 && frame->joint->qIndex >= 0 && frame->joint->qIndex < x.N) {
  //       double angle = x(frame->joint->qIndex);
  //       jointData.AddMember("angle", angle, allocator);
  //     }
      
  //     robot.AddMember("joint", jointData, allocator);
  //   }
    
  //   // Export geometry type
  //   auto geom = obj->collisionGeometry();
  //   robot.AddMember("geometry_type", geom->getNodeType(), allocator);
    
  //   // Handle different geometry types (same as for obstacles)
  //   rapidjson::Value geometry(rapidjson::kObjectType);
    
  //   if (geom->getNodeType() == fcl::GEOM_BOX) {
  //     const Box* box = static_cast<const Box*>(geom.get());
      
  //     rapidjson::Value typeValue;
  //     typeValue.SetString("box", allocator);
  //     geometry.AddMember("type", typeValue, allocator);
      
  //     rapidjson::Value sizeArray(rapidjson::kArrayType);
  //     sizeArray.PushBack(box->side[0], allocator);
  //     sizeArray.PushBack(box->side[1], allocator);
  //     sizeArray.PushBack(box->side[2], allocator);
  //     geometry.AddMember("size", sizeArray, allocator);
  //   }
  //   else if (geom->getNodeType() == fcl::GEOM_SPHERE) {
  //     const Sphere* sphere = static_cast<const Sphere*>(geom.get());
      
  //     rapidjson::Value typeValue;
  //     typeValue.SetString("sphere", allocator);
  //     geometry.AddMember("type", typeValue, allocator);
      
  //     geometry.AddMember("radius", sphere->radius, allocator);
  //   }
  //   else if (geom->getNodeType() == fcl::GEOM_CAPSULE) {
  //     const Capsule* capsule = static_cast<const Capsule*>(geom.get());
      
  //     rapidjson::Value typeValue;
  //     typeValue.SetString("capsule", allocator);
  //     geometry.AddMember("type", typeValue, allocator);
      
  //     geometry.AddMember("radius", capsule->radius, allocator);
  //     geometry.AddMember("height", capsule->lz, allocator);
  //   }
  //   else if (geom->getNodeType() == fcl::GEOM_CYLINDER) {
  //     const Cylinder* cylinder = static_cast<const Cylinder*>(geom.get());
      
  //     rapidjson::Value typeValue;
  //     typeValue.SetString("cylinder", allocator);
  //     geometry.AddMember("type", typeValue, allocator);
      
  //     geometry.AddMember("radius", cylinder->radius, allocator);
  //     geometry.AddMember("height", cylinder->lz, allocator);
  //   }
  //   else if (geom->getNodeType() == fcl::GEOM_CONVEX) {
  //     rapidjson::Value typeValue;
  //     typeValue.SetString("convex", allocator);
  //     geometry.AddMember("type", typeValue, allocator);
  //     // Would need to extract vertices and faces for complete export
  //   }
    
  //   robot.AddMember("geometry", geometry, allocator);
  //   robotsArray.PushBack(robot, allocator);
    
  //   // Clean up the temporary robot object
  //   delete obj;
  // }
  
  // // Add arrays to document
  // document.AddMember("obstacles", obstaclesArray, allocator);
  // document.AddMember("robots", robotsArray, allocator);
  
  // // Create the output directory if it doesn't exist
  // std::filesystem::create_directories("fcl_export");
  
  // // Generate a unique filename with timestamp
  // std::string timestamp = std::to_string(time(nullptr));
  // std::string filename = "fcl_export/safe_intervals_" + timestamp + ".json";
  
  // // Write to file with pretty printing
  // rapidjson::StringBuffer buffer;
  // rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
  // document.Accept(writer);
  
  // std::ofstream file(filename);
  // file << buffer.GetString();
  // file.close();
  
  // std::cout << "Exported current scene state to " << filename << std::endl;
  
  // check collision with static obstacles
  // if(this->collision_with_static()){
  //   return results;
  // }
  //update robot link objects pos

  // assert(C.activeJoints.N == this->robot_link_objects.size());
    // std::cout<<"Robot joint pos"<<std::endl;

  for(int joint_id=0;joint_id<this->robot_link_objects.size();joint_id++){
    // std::cout<<Vec3f(this->robot_link_objects[joint_id].first->X.pos.x, this->robot_link_objects[joint_id].first->X.pos.y,this->robot_link_objects[joint_id].first->X.pos.z)<<std::endl;
    this->robot_link_objects[joint_id].first->ensure_X();
    this->robot_link_objects[joint_id].second->setTranslation(Vec3f(this->robot_link_objects[joint_id].first->X.pos.x, this->robot_link_objects[joint_id].first->X.pos.y,this->robot_link_objects[joint_id].first->X.pos.z));
    this->robot_link_objects[joint_id].second->setQuatRotation(Quaternionf(this->robot_link_objects[joint_id].first->X.rot.w,this->robot_link_objects[joint_id].first->X.rot.x,this->robot_link_objects[joint_id].first->X.rot.y,this->robot_link_objects[joint_id].first->X.rot.z));
    this->robot_link_objects[joint_id].second->computeAABB();
  }

  // check for the selfcollision
  for (int robot_joint_id = 0; robot_joint_id < this->robot_link_objects.size(); robot_joint_id++)
  {
    for (int another_robot_joint_id = robot_joint_id+1; another_robot_joint_id < this->robot_link_objects.size(); another_robot_joint_id++)
    {
      if ((!this->robot_link_objects[another_robot_joint_id].first->shape->canCollideWith(this->robot_link_objects[robot_joint_id].first))){
        continue;
      }
      CollisionRequest request;
      CollisionResult result;
      fcl::collide(this->robot_link_objects[robot_joint_id].second, this->robot_link_objects[another_robot_joint_id].second, request, result);
      if (result.isCollision())
      {
        // std::cout<<"Collision between "<<this->robot_link_objects[robot_joint_id].first->name<<" and "<<this->robot_link_objects[another_robot_joint_id].first->name<<std::endl;
        return results;
      }
    }
  } 


  //collide
  this->collision_moments.clear();
  
  for (int robot_joint_id = 0; robot_joint_id < this->robot_link_objects.size(); robot_joint_id++)
  {
      this->safe_interval_collision_manager.collide(this->robot_link_objects[robot_joint_id].second, this, this->BroadphaseCallback);
  }
    // for (const int collision_frame : collision_moments)
    // {
    //   std::cout<< "Collision_frame "<<collision_frame<<std::endl;
    // }
    // compute safe intervals

    if (collision_moments.size() == 0)
    {
      // std::cout<<"all safe"<<std::endl;
        results.emplace_back(this->min_time, this->max_time);
        return results;
    }
    // sort and delete duplicates
    std::sort(collision_moments.begin(), collision_moments.end());
    collision_moments.erase(std::unique(collision_moments.begin(), collision_moments.end()), collision_moments.end());

    int last_collision_frame = this->min_time-1;
    for (const int collision_frame : collision_moments)
    {
        if ((collision_frame - last_collision_frame) > 1)
        {
            results.emplace_back(last_collision_frame + 1, collision_frame - 1);
        }

        last_collision_frame = collision_frame;
    }
    if (this->max_time - collision_moments.back()  > 1)
    {
        results.emplace_back(collision_moments.back() + 1,  this->max_time);
    }
    return results;
}

bool TimedConfigurationProblem::BroadphaseCallback(CollObject* o1, CollObject* o2, void* cdata_) {
  TimedConfigurationProblem* self = static_cast<TimedConfigurationProblem*>(cdata_);


  CollisionRequest request;
  CollisionResult result;
  fcl::collide(o1, o2, request, result);
  if(result.isCollision()) {

    double d1,d2,d3;
    d1 = *((double*)(((std::pair<double*,rai::Frame*>*)(o1->getUserData()))->first));
    d2 = *((double*)(((std::pair<double*,rai::Frame*>*)(o2->getUserData()))->first));
    d3 = d1;
    if(d1 < 0){
      d3 = d2;
    }

    // std::cout<<"Collision between"<<std::endl;
    // std::cout<<((std::pair<double*,rai::String>*)(o1->getUserData()))->second<<" "<<d1<<std::endl;
    // std::cout<<((std::pair<double*,rai::Frame*>*)(o2->getUserData()))->second->name<<" "<<d2<<std::endl;
    
    
  
    // fcl::DistanceRequest<float> requestd;
    // fcl::DistanceResult<float> resultd;
    // fcl::distance(o1, o2, requestd, resultd);
    // std::cout<<"Distance: "<<resultd.min_distance<<std::endl;
    self->collision_moments.push_back(d3);


    
  }

  return false;
}

// NOTE: this asumes that we are doing inf-norm sampling
arr TimedConfigurationProblem::sample(const arr &start, const arr &goal, const double c_max, const double c_min){
  const arr limits = C.getLimits();
  const uint dim = limits.dim(0);
  arr sample(dim);

  if (start.d0 == 0){
    // sample uniformly between 0,1
    cp_localRndUniform(sample,0,1,false);

    // scale sample
    for (uint i=0; i<sample.d0; ++i){
      if(limits(i,1) > limits(i,0)){
        sample(i) = limits(i, 0) + sample(i) * (limits(i, 1) - limits(i, 0));
      }
      else {
        // default: [-5, 5]
        sample(i) = sample(i) * 10 - 5;
      }
    }
  }
  else{
    // this approach is sampling from an overapproximation of the actual space that can improve the solution
    // actually, we would wan tto sample from a polytope of a very specific shape. This is hard tho.
    
    // this assumes that we are dealing with an infty norm
    const arr midpoint = (start + goal) / 2.;

    for (uint i=0; i<sample.d0; ++i){
      double lo = midpoint(i) - c_max * 0.5;
      double hi = midpoint(i) + c_max * 0.5;
  
      // if the diff between both limits is smaller than 1e-6, we treat this as not having a limit
      if (abs(limits(i, 0) - limits(i, 1)) > 1e-6) {
        lo = std::max(lo, limits(i, 0));
        hi = std::min(hi, limits(i, 1));
      }
    
      sample.p[i] = (double)rnd2_.uni(lo, hi);
    }
  }
  return sample;
}


double corput(int n, const int base=2){
  double q=0, bk=(double)1/base;

  while (n > 0) {
    q += (n % base)*bk;
    n /= base;
    bk /= base;
  }

  return q;
}

bool TimedConfigurationProblem::checkEdge(const arr& x0, const double t0, const arr &x1, const double t1, const uint discretization){
  const double l = length(x0 - x1);
  uint samples = discretization * l;
  if (samples < 5) { samples = 2;}

  samples = std::ceil(std::abs(t0 - t1));
  
  const double tMin = (t0 < t1) ? t0: t1;

  // #pragma omp parallel for
  for(uint i=0; i<=samples; ++i){
    // const double interp = corput(i);
    const double interp = 1. * (i) / (samples);
    const arr x = x0 + interp * (x1-x0);
    double t = t0 + interp * (t1 - t0);

    const auto res = this->query(x, t, tMin);
    if (!res->isFeasible){ return false;}
  }

  return true;
}

void QueryResult::getViolatedContacts(arr& y, arr& J, double margin){
  uintA violated;
  for(uint i=0;i<coll_y.N;i++) if(coll_y.elem(i)<margin) violated.append(i);

  if(!violated.N){
    y.resize(0);
    J.resize(0,coll_J.d1);
  }else{
    y = coll_y.sub(violated);
    J = coll_J.sub(violated);
  }

}

arr QueryResult::getSideStep(){
  arr s = randn(3);
  s /=length(s);

  arr S(side_J.d0, 3);
  for(uint i=0;i<S.d0;i++) S[i] = s;

  arr J = side_J;

  S.reshape(-1);
  J.reshape(S.N, -1);

#if 0
  arr U, sig ,V;
  svd(U, sig, V, J);
  arr d = ~V * sig % V * randn(V.d1); //random step in input space of J!
#else
  arr JI = ~J; //pseudoInverse(J);
  arr d = JI * S;
#endif

  if(length(d)<1e-10) HALT("???");

  return d;
}

arr QueryResult::getForwardStep(){
  arr goal_JI = pseudoInverse(goal_J);
  arr d = goal_JI * (-goal_y);
  return d;
}

arr QueryResult::getBackwardStep(double relativeStepLength, double margin, const arr& nullStep){
//  CHECK(!isFeasible, "");
  CHECK(coll_y.N>0, "");

  arr y,J;
  getViolatedContacts(y, J, margin);
  y -= margin;

  arr Jinv = pseudoInverse(J, NoArr, 1e-4);
  arr d = Jinv * (-relativeStepLength * y);

  if(!!nullStep) d += (eye(J.d1) - Jinv * J) * nullStep;

  return d;
}

void QueryResult::write(std::ostream& os) const{
  os <<"query: h_goal: " <<sumOfAbs(goal_y)
    <<" g_coll: " <<sum(elemWiseHinge(-coll_y))
   <<" isGoal: " <<isGoal
  <<" isFeasible: " <<isFeasible;
}

void QueryResult::writeDetails(std::ostream& os, const ConfigurationProblem& P, double margin) const{
  write(os);
  for(uint i=0;i<coll_y.N;i++){
    if(coll_y.elem(i)<margin && computed(i)){ // make sure that we do not print 'uncomputed' things
      os <<"\ncoll " <<i <<':' <<collisions[i]
           <<':' <<P.C.frames(collisions(i,0))->name <<'-' <<P.C.frames(collisions(i,1))->name
          <<" y:" <<coll_y.elem(i) <<" normal:" <<normal_y[i];
    }
  }
  os <<endl;
}

bool makePoseFeasible(arr& x, ConfigurationProblem& P, double IKstepSize, double maxQStepSize, uint trials){
  shared_ptr<QueryResult> qr = P.query(x);
  for(uint k=0;k<trials;k++){
    if(qr->isFeasible){
      break;
    }else{
    }
    arr delta = qr->getBackwardStep(IKstepSize);
    double l = length(delta);
    if(maxQStepSize>0. && l>maxQStepSize) delta *= maxQStepSize/l;
    x += delta;
    qr = P.query(x);
  }
  return qr->isFeasible;
}


#include <Kin/proxy.h>
#include <Optim/constrained.h>

#include "ConfigurationProblem.h"

#include <Geo/fclInterface.h>

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
        // LOG(-1) <<"QUERY OUT OF LIMIT: joint " <<i <<": " <<x.elem(i) <<' ' <<limits[i];
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

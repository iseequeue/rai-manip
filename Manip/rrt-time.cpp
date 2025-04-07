#include "rrt-time.h"

#include <chrono>
#include <memory>

#include <thread>
#include <mutex>

#define PERIODIC

arr PathFinder_RRT_Time::getDelta(const arr &p1, const arr &p2){
#ifdef PERIODIC
  arr tmp(delta_buffer.N);
  // delta_buffer = p1 - p2;
  // for (uint i=0; i<p1.d0; ++i){delta_buffer(i) = p1(i) - p2(i);}
  for (uint i=0; i<delta_buffer.N; ++i){
    if (periodicDimensions[i]){
      // this is an angular joint -> we need to check the other direction
      const double start = p1.p[i];
      const double end = p2.p[i];
      tmp.p[i] = std::fmod(start - end + 3.*RAI_PI, 2*RAI_PI) - RAI_PI;
    }
    else{
      tmp.p[i] = p1.p[i] - p2.p[i];
    }
  }
#else
  const arr delta = p1 - p2;
#endif

  return tmp;
}

double PathFinder_RRT_Time::q_metric(const arr& d) const{
#if 1
  double dist = 0;
  for (auto *j: TP.C.activeJoints){
    //double tmp = length(d({j->qIndex, j->qIndex+j->dim-1}));
    // const double tmp = absMax(d({j->qIndex, j->qIndex+j->dim-1}));
    double tmp = 0;
    for (uint i=j->qIndex; i < j->qIndex+j->dim; ++i){
      tmp = std::max({tmp, abs(d(i))});
    }
    // std::cout << tmp << " " << absMax(d({j->qIndex, j->qIndex+j->dim-1})) << std::endl;
    dist = std::max({tmp, dist});
    // if (tmp > dist){dist = tmp;}
  }
  /*for (uint i=0; i<d.N; ++i){ 
    if(std::fabs(d(i)) > dist) {dist = std::fabs(d(i));}
  }*/
  return dist;
#else
  return length(d);
#endif
}

double PathFinder_RRT_Time::distance(const Node &n1, const Node &n2){
  // gives back the distance to go from n1 to n2
  const double dt = n2.t - n1.t;
  if (dt < 0){
    return std::numeric_limits<double>::infinity();
  }
  
  const arr tmp = getDelta(n2.q, n1.q);
  const double dist = q_metric(tmp);

  if (false){
    std::cout << dist << " " << dt << " " << dist / dt << std::endl;
  }

  if (dt >= 0 && dt < 1e-6 && dist < 1e-6){
    return 0;
  }

  if (dist / dt > vmax){
    return std::numeric_limits<double>::infinity();
  }

  return  lambda * dist + (1-lambda) * dt;
}

Node* PathFinder_RRT_Time::steer(const Node &start, const Node &goal, bool reverse){
  const arr delta = getDelta(goal.q, start.q);
  const double dist = q_metric(delta);
  double dt = std::fabs(goal.t - start.t);

  // std::cout << "goal.t: " << goal.t << std::endl;
  // std::cout << "start.t: " << start.t << std::endl;
  // std::cout << "dt: " << dt << std::endl;

  //std::cout << start.q << std::endl << goal.q << std::endl << dir << std::endl << goal.q - start.q << std::endl;

  const double u = dist / dt;

  if (u >= vmax){
    std::cout << "reversed: " << reverse << std::endl;
    std::cout << "t start: " << start.t << std::endl;
    std::cout << "t goal: " << goal.t << std::endl;
    std::cout << "dt: " << dt << std::endl;
    std::cout << "dist: " << dist << std::endl;

    CHECK(false, "Violated vmax");
  }

  if (dt > step_time){
    dt = step_time;
  }

  double t;
  if (!reverse){
    t = start.t + dt;
  }
  else{
    t = start.t - dt;
  }
  //std::cout << "t " << t << std::endl;

  arr q = start.q + delta / dist * u * dt;

  // go over all pts to make sure they are periodic
#ifdef PERIODIC
  for (uint i=0; i<q.N; ++i){
    if (periodicDimensions[i]){
      if (q(i) >= RAI_PI){ q(i) -= 2*RAI_PI; }
      if (q(i) <= -RAI_PI){ q(i) += 2*RAI_PI; }
    }
  }
#endif

  Node* n = new Node(q, t);

  if (prePlannedFrames.N != 0){
    arr tmp = projectToManifold(q, t);
    n->q = tmp * 1.;
  }

  return n;
}

arr PathFinder_RRT_Time::projectToManifold(const arr &q, const double t){
  // need to know what frames are active
  // set them to the previously planned path
  // if there is something planned for them

  if (prePlannedFrames.N == 0){
    return q;
  }

  // handle the mode switch-crap here
  //TP.A.setToTime(TP.C, tPrePlanned);

  //std::cout << t << std::endl;

  if (t > tPrePlanned){
    return q;
  }

  // TP.query sets the state to the constrained one automatically
  TP.query(q, t);
  arr qNew = TP.C.getJointState();
  if (false){
    //Ccpy.setJointState(qPre, _joints);
    //Ccpy.watch(true);

    TP.C.setJointState(q);
    TP.C.watch(true);

    TP.C.setJointState(qNew);
    TP.C.watch(true);
  }
  return qNew * 1.;
}

Node* PathFinder_RRT_Time::extend(Tree* tree, const Node &goal, const bool connect){
  const double resolution = 0.05;
  
  if(!connect){
    Node* close = tree->getNearest(goal);
    if (!close){return nullptr;}

    Node *n = steer(*close, goal, tree->reverse);
    close = tree->getNearest(*n);
    
    const auto state_res = TP.query(n->q, n->t);

    if (state_res->isFeasible){
      const auto start_time = std::chrono::high_resolution_clock::now();

      // TODO: this way of checking an edge is bad, the resolution should not always be 20 pts
      const uint discretization = std::ceil(length(close->q - n->q) / resolution);
      // std::cout << discretization << std::endl;
      const bool edge_valid = TP.checkEdge(close->q, close->t, n->q, n->t, discretization);

      const auto end_time = std::chrono::high_resolution_clock::now();
      const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time -start_time).count();

      edge_checking_time_us += duration;

      if (edge_valid){
        n->parent = close;
        close->addChild(n);

        // double d;
        // if (!tree->reverse) {d = distance(*close, *n);}
        // else {d = distance(*n, *close);}
        //n->cost = close->cost + d;
        n->cost = close->cost + abs(close->t - n->t);
        n->x = state_res->disp3d;

        tree->addNode(n);
        return n;
      }
    }
    delete n;
    return nullptr;
  }
  else{
    Node* prev = nullptr;
    //std::cout << "Start extending" << std::endl;
    double prevDist = -1.;

    uint cnt = 0;
    while(true){
      Node * close = tree->getNearest(goal);
      if (!close){return nullptr;}

      // early exit
      if (cnt > 50){
        return prev;
      }

      Node *n = steer(*close, goal, tree->reverse);
      close = tree->getNearest(*n);

      // if (prev) std::cout << "prev" << prev->q << std::endl;
      // if (close) std::cout << close->q << std::endl;

      if (!close){
        // std:cout << "AAAAAAAAAAAA" << std::endl;
        delete n;
        return prev;
      }

      // TODO: make sure that this "close" is fine with the normal planning
      auto state_res = TP.query(n->q, n->t);
      // TP.C.watch(true);

      if (!state_res->isFeasible){
        // std::cout << "not feas" << std::endl;

        delete n;
        return prev;
      }

      const auto start_time = std::chrono::high_resolution_clock::now();

      const uint discretization = std::ceil(length(close->q - n->q) / resolution);
      bool edge_valid = TP.checkEdge(close->q, close->t, n->q, n->t, discretization);

      const auto end_time = std::chrono::high_resolution_clock::now();
      const auto duration =
          std::chrono::duration_cast<std::chrono::microseconds>(end_time -
                                                                start_time)
              .count();

      edge_checking_time_us += duration;

      if (edge_valid){
        n->parent = close;
        close->addChild(n);

        double d;
        if (!tree->reverse) {d = distance(*close, *n);}
        else {d = distance(*n, *close);}

        if (d > prevDist && prevDist > 0){
          // TODO: finish this implementation
          //break;
        }

        //n->cost = close->cost + d;
        n->cost = close->cost + abs(close->t - n->t);
        n->x = state_res->disp3d;

        tree->addNode(n);
        prev = n;

        // std::cout << "dist: " <<  distance(*n, goal) << " " <<  distance(goal, *n) << std::endl;
        if (distance(*n, goal) < tol || distance(goal, *n) < tol){
          return n;
        }

        prevDist = d;
      }
      else{
        // std::cout << "BBBBBBBBBBBB" << std::endl;
        if(true || verbose){
          // res->writeDetails(cout, TP.C);
          // std::cout << "collisiton at time: " << n->t << std::endl;
          // TP.A.setToTime(TP.C, n->t);
          // TP.C.setJointState(n->q);
          // TP.C.watch(true);
        }
        delete n;
        return prev;
      }
      ++cnt;
    }
  }
}

TimedPath PathFinder_RRT_Time::plan(const arr &q0, const double t0, const TimedGoalSampler gs, double tGoalLowerBound, double tGoalUpperBound){
  const bool fixedTime = rai::getParameter<bool>("assembly/fixedTime", false); 

  //auto start_rrt_time = std::chrono::high_resolution_clock::now();

  // check initial and final nodes
  if (!TP.query(q0, t0, 0)->isFeasible){
    spdlog::error("Initial point not feasible");
    TP.C.watch(true);
    return TimedPath({}, {});
  }
  
  // TODO: check final node without dynamic things

  // scale with the DoF of the object we are steering
  const double time_multiple = std::sqrt(q0.N);

  if (verbose){
    std::stringstream ss;
    ss << q0;
    spdlog::info("q0: {}", ss.str());
    //std::cout << "qT: " << qT << std::endl;
    //std::cout << "mindt: " << min_dt << std::endl;
  }

  // TODO: deal with this properly
  arr qT;
  for (uint i=0; i<3; ++i){
    arr rnd(1);
    rndUniform(rnd, 0, 1);
    int t = ((t0 > tGoalLowerBound)?t0:tGoalLowerBound) + 100. * rnd(0);

    if (fixedTime){
      t = tGoalLowerBound;
    }

    gs(t, qT);

    if (qT.N != 0){break;}
  }

  if (qT.N == 0){
    spdlog::error("Failed bc. no initial goal was found");
    return TimedPath({}, {});
  }

  // TODO: approximate min_l by sampling without dynamic obstacles:
  // this gives an estimate of the lower bound
  const double min_l = q_metric(getDelta(qT, q0));
  const double min_dt = min_l / vmax;

  double min_time = t0 + min_dt;
  //double min_time = t0 + min_dt * 1.5;
  if (tGoalLowerBound > min_time){
    spdlog::info("Replacing min time {} (from max speed) with given lower bound {}", min_time, tGoalLowerBound);
    min_time = tGoalLowerBound;
  }

  double max_time = t0 + time_multiple * min_dt;
  if (max_time < min_time){
    max_time = min_time + time_multiple * min_dt;
  }

  if (max_time - min_time < 10){
    max_time = min_time + 10;
  }

  if (tGoalUpperBound < max_time && tGoalUpperBound > 0){
    spdlog::info("Replacing computed max time {} with upper bound {}", max_time, tGoalUpperBound);
    max_time = tGoalUpperBound;
  }

  //std::cout << "q0: " << q0 << std::endl;
  //std::cout << "qT: " << qT << std::endl;
  spdlog::info("t0: {}", t0);
  spdlog::info("q_metric: {}", min_l);
  spdlog::info("min_dt: {}", min_dt);
  spdlog::info("min_time: {}", min_time);
  spdlog::info("max_time: {}", max_time);

  // setup trees
  auto distFun = std::bind(&PathFinder_RRT_Time::distance, this, std::placeholders::_1, std::placeholders::_2);

  Tree startTree(distFun);
  {
    // sample start
    // const arr q = q0;

    const arr q = projectToManifold(q0, t0);

    Node *start_node = new Node(q, t0);
    start_node->cost = 0;
    start_node->x = TP.query(q0, t0)->disp3d;
    startTree.addNode(start_node);
  }

  rai::Array<std::pair<arr, double>> goals;

  // setup goal tree
  int max_goal_time = 0;

  spdlog::info("Setting up goal tree");
  Tree goalTree(distFun, true);
  for(uint i=0;;++i){
    if (i > maxInitialSamples || (fixedTime && i > 10)){
      spdlog::error("No final nodes found");
      //HALT("Q");
      return TimedPath({}, {});
    }

    // increase goal-time-range
    max_time += time_multiple;
    if (max_time > tGoalUpperBound && tGoalUpperBound > 0){ max_time = tGoalUpperBound;}

    // sample goal time
    arr rnd(1);
    rndUniform(rnd, 0, 1);

    // TODO: figure out what to do here
    int t = min_time + 1. * rnd(0) * (max_time - min_time) * i / 10.;
    //double t = min_time + 1. * rnd(0) * (max_time - min_time);

    if (fixedTime){
      t = tGoalLowerBound;
      std::cout << "A" << std::endl;
      std::cout << "A" << std::endl;
      std::cout << "switchng to fixed time" << std::endl;
      std::cout << "A" << std::endl;
      std::cout << "A" << std::endl;
    }
    //std::cout << t << " " << rnd(0) << " " << min_time << " " << max_time << " " << max_time - min_time << std::endl;

    // sample goal at time t
    arr q; 
    gs(t, q);
    
    if(q.N == 0){continue;} // infeasible sample

    //q = projectToManifold(q, t);
    
    //std::cout << q << std::endl;

    //TP.C.setJointState(q);
    //TP.C.watch(true);

    auto res = TP.query(q, t, t0);
    if (q_metric(getDelta(q, q0))/vmax > t - t0){
      spdlog::info("Final point not reachable due to v_max. Dt: {}, dist/vmax: {}", t-t0, q_metric(getDelta(q, q0))/vmax);
      continue;
    }
    if (!res->isFeasible){
      //res->writeDetails(cout, TP.C);
      spdlog::error("Final point not feasible at time {}", t);
      res->writeDetails(cout, TP.C);
      // TP.C.watch(true);
      
      continue;
    }

    spdlog::info("Adding goal at time {}", t);
    // res->writeDetails(cout, TP.C);
    // ConfigurationProblem cp(TP.C);
    // {
    //   auto res = cp.query(q);

    //   if (res->isFeasible == false) {
    //     res->writeDetails(cout, cp.C);
    //     cp.C.watch(true);

    //     TP.C.setJointState(q);
    //     TP.C.watch(true);
    //   }
    // }
    Node *goal_node = new Node(q, t);
    goal_node->x = res->disp3d;
    goal_node->cost = 0;
    goalTree.addNode(goal_node);

    goals.append({q, t});
    max_goal_time = std::max({max_goal_time, t});

    //TP.C.setJointState(q);
    //TP.C.watch(true);

    break;
  }
  
  Tree *ta = &startTree;
  Tree *tb = &goalTree;

  rai::Configuration DSP;
  DSP.copy(TP.C);
  DSP.gl()->add(startTree);
  DSP.gl()->add(goalTree);
  if (disp){
    DSP.copy(TP.C);
    DSP.gl()->add(startTree);
    DSP.gl()->add(goalTree);

    DSP.watch(false);
  }

  TimedPath tp({}, {});

  bool success = false;
  bool converged = false;

  std::vector<double> costs;
  std::vector<uint> iters;

  double minDist = length(getDelta(q0, qT));

  std::array<Node*, 2> finalNodes = {nullptr, nullptr};

  uint cnt = 0;
  int cnt2 = 0;
  while(true){
    if (disp){
      DSP.watch(false);
    }
    // cnt2++;
    // if (cnt2 %1000 == 0) { std::cout << cnt2 << std::endl; }

    if (cnt%100==0){
      spdlog::info("iteration: {}", cnt);
      spdlog::info("nodes in start tree: {} nodes in goal tree: {}, min dist between trees {}", startTree.nodes.size(), goalTree.nodes.size(),  minDist);
    }

    if (cnt > maxIter){
      // Display the direct interpolation
      /*const arr delta = getDelta(qT, q0);
      for (uint i=0; i<20; ++i){
        const arr p = q0 + 1. * i / 20 * delta;
        TP.A.setToTime(TP.C, t0 + i);
        TP.C.setJointState(p);
        TP.C.watch(true);

        auto res = TP.query(p, t0 + i);
        if (!res->isFeasible){std::cout << "Direct interpolation not possible" << std::endl;}
      }*/

      /*{
        std::ofstream f;
        f.open("./out/not_feasible.txt", std::ios_base::app);
        f << "q0: " << q0 << std::endl;
        f << "qT: " << qT << std::endl;
      }*/

      if (!success){
        // TP.C.setJointState(qT);
        // TP.C.watch(true);

        // TP.A.setToTime(DSP, tb_new->t);
        // DSP.setJointState(tb_new->q);
        // DSP.watch(true);

        spdlog::debug("Didnt find a path in the maximum number of iterations.");
        return TimedPath({}, {});
      }
      else{
        break;
      }
    }

    /*if (tb->nodes.size() + ta->nodes.size() > 200){
      if ((tb->nodes.size() > ta->nodes.size() && tb->nodes.size() > 10*ta->nodes.size()) || 
          (ta->nodes.size() > tb->nodes.size() && ta->nodes.size() > 10*tb->nodes.size())){
        std::cout << "Too many nodes without arriving at goal" << std::endl;
        return TimedPath({}, {});
      }
    }*/

    // swap trees
    Tree *tmp = ta;
    ta = tb;
    tb = tmp;

    ++cnt;
    
    arr r(1);  //dim0 = 1 - [x]?
    rndUniform(r, 0, 1); // [random x between 0 and 1]?
    if (r(0) > goalSampleProbability){// && goals.N < 10){  // indexation
      // sampleGoal();

      // increase goal-time-range
      //adjustTimeRange(min_time, max_time, tGoalLowerBound, tGoalUpperBound);
      max_time += time_multiple;
      //min_time -= time_multiple; // TODO: think about how this should be done

      if (max_time > tGoalUpperBound && tGoalUpperBound > 0){ max_time = tGoalUpperBound;}
      if (min_time < tGoalLowerBound) {min_time = tGoalLowerBound;}
      if (costs.size() > 0){max_time = (costs.back() - lambda*min_l) / (1-lambda) ;}

      // sample goal time
      arr rnd(1);
      rndUniform(rnd, 0, 1);
      int t = min_time + rnd(0) * (max_time - min_time);

      if (t < tGoalLowerBound){t = tGoalLowerBound;}

      if (fixedTime){
        t = tGoalLowerBound;
      }

      // sample goal at time t
      arr q;
      gs(t, q);
      if(q.N == 0){continue;} // infeasible sample
      //q = projectToManifold(q, t);
      //std::cout << q << std::endl;

      //TP.C.setJointState(q);
      //TP.C.watch(true);

      // TODO: first check if this is feasible in the static env
      auto res = TP.query(q, t, t0);
      if (!res->isFeasible){
        //res->writeDetails(cout, TP.C);
        std::cout << "Final point not feasible at time " << t << " (in loop)" << std::endl;
      }
      else if (q_metric(getDelta(q0, q)) / vmax > t - t0){
        std::cout << "Final pt. not reachable due to speed" << std::endl;
      }
      else{
        //std::cout << "adding in loop at time " << t << std::endl;
        {
          Node *goal_node = new Node(q, t);
          goal_node->cost = 0;
          goal_node->x = res->disp3d;
          goalTree.addNode(goal_node);
          goals.append({q, t});

          max_goal_time = std::max({max_goal_time, t});

          // spdlog::info("Adding additional goal at time {}", t);
        }

        // if a goal is found, we sample it additionally at various times!
        if (sampleAdditionalNodes){
          for (uint i=0; i<10; ++i){
            uint t_ = min_time + (max_time - min_time)/10. * i;
            auto res_ = TP.query(q, t_, t0);
            if (!res_->isFeasible){
              Node *goal_node = new Node(q, t);
              goal_node->x = res->disp3d;
              goal_node->cost = 0;
              goalTree.addNode(goal_node);
              goals.append({q, t});

              max_goal_time = std::max({max_goal_time, t});
            }
          }
        }
      }
    }

    arr qs;
    double ts;
    if (sampling_type == SamplingType::BOX_CONSTRAINED_CONDITIONAL_SAMPLING){
      if (informed_sampling) {
        qs = TP.sample(q0, qT, (max_goal_time - t0) * vmax, min_l);
      } else {
        qs = TP.sample();
      }

      // sample t
      const double min_dt_from_start = q_metric(getDelta(qs, q0)) / vmax;
      
      double max_t_sample = 0;
      for (const auto &g: goals){
        const double goal_time = g.second;
        const arr &goal_pose = g.first;
        
        const double tmp = goal_time - q_metric(getDelta(qs, goal_pose)) / vmax;
        //std::cout << "time to goal: " << tmp << std::endl;
        if (tmp > max_t_sample) {max_t_sample = tmp;}
      }

      //std::cout << max_goal_time << std::endl;
      //std::cout << min_dt_from_goal << std::endl;

      const double min_t_sample = t0 + min_dt_from_start;

      arr ts_arr(1);
      rndUniform(ts_arr, min_t_sample, max_t_sample);
      ts = ts_arr(0);

      //std::cout << "sampling between " << min_t_sample << " " << max_t_sample << std::endl;

      if (min_t_sample > max_t_sample){
        //std::cout << min_t_sample << " " << max_t_sample << std::endl;
        // std::cout << "Rejected sample" << std::endl;
        //HALT("A");
        continue;
      }
    }
    else{
      spdlog::info("Rejection sampling");

      bool found_valid_sample = false;
      while (true){
        qs = TP.sample(q0, qT, (max_goal_time - t0) * vmax, min_l);
        const double min_dt_from_start = q_metric(getDelta(qs, q0)) / vmax;

        arr ts_arr(1);
        rndUniform(ts_arr, t0, max_goal_time);
        ts = ts_arr(0);

        if (ts >= t0 + min_dt_from_start){
          for (const auto &g: goals){
            const double goal_time = g.second;
            const arr &goal_pose = g.first;
            
            const double tmp = goal_time - q_metric(getDelta(qs, goal_pose)) / vmax;

            if (ts <= tmp){
              found_valid_sample = true;
              break;
            }
          }

          if (found_valid_sample){
            // spdlog::info("Found sample at time {}", ts);
            break;
          }
        }
      }
    }

    if (prePlannedFrames.N != 0 && ts <= tPrePlanned){
      // std::cout << ts << std::endl;
      // std::cout << qs << std::endl;
      qs = projectToManifold(qs, ts);
      // std::cout << qs << std::endl << std::endl;
    }

    if (verbose){
      std::cout << "Sampled time: " << ts << std::endl;
      std::cout << "Sampled pos: " << qs << std::endl;
    }

    // try to connect ta to the sample
    Node sampled(qs, ts);    
    // steer from close towards sampled
    spdlog::trace("Added sample at time {}", ts);
    Node* ta_new = extend(ta, sampled, connect);

    if (ta_new){
      spdlog::trace("Successfully extended");
      // TP.C.watch(true);
      if (disp){
        TP.A.setToTime(DSP, ta_new->t);
        DSP.setJointState(ta_new->q);
        DSP.watch(true);
      }
      if (optimize){ ta->rewire(ta_new, TP);}

      // attempt to connect to the node from the other tree
      Node* tb_new = extend(tb, *ta_new, connect);

      if(tb_new){ 
        spdlog::trace("Successfully connected");
        if (disp){
          TP.A.setToTime(DSP, tb_new->t);
          DSP.setJointState(tb_new->q);
          DSP.watch(true);
        }

        if (optimize){ tb->rewire(tb_new, TP);}

        const double dist = length(getDelta(ta_new->q, tb_new->q));
        if (minDist > dist) {minDist = dist;}

        //std::cout << distance(*tb_new, *ta_new) << " " <<  distance(*ta_new, *tb_new) << std::endl;
        const double new_cost = tb_new->cost + ta_new->cost;
        if ((costs.size() == 0 || finalNodes[0]->cost + finalNodes[1]->cost > new_cost) &&
            (distance(*tb_new, *ta_new) < tol || distance(*ta_new, *tb_new) < tol)){
          success = true;
          tp = extractPath(ta, tb, ta_new, tb_new);

          costs.push_back(new_cost);
          iters.push_back(cnt);

          finalNodes[0] = ta_new;
          finalNodes[1] = tb_new;

          spdlog::info("Path found, cost: {}", new_cost);
        }
      }
    }

    if (finalNodes[0] && finalNodes[0]->cost + finalNodes[1]->cost < costs.back()){
      spdlog::info("Found cheaper path: {}", finalNodes[0]->cost + finalNodes[1]->cost);

      costs.push_back(finalNodes[0]->cost + finalNodes[1]->cost);
      iters.push_back(cnt);
    }

    // check convergence criterion
    if (iters.size() > 0 && cnt - iters.back() > conv_iter){
      spdlog::info("Stopped due to iterations");
      converged = true;
    }

    if (costs.size() > 0 && costs.back() < 1.1 * (tGoalLowerBound - t0)){
      spdlog::info("Converged due to cost" );
      converged = true;
    }

    if (success && (!optimize || converged)){
      break;
    }
  }

  //auto end_rrt_time = std::chrono::high_resolution_clock::now();
  //auto duration = std::chrono::duration_cast<std::chrono::microseconds>( end_rrt_time - start_rrt_time ).count();

  //std::cout << time << std::endl;
  //std::cout << path << std::endl;
  //std::cout << "done in " << tp.time(tp.time.N-1) - tp.time(0) << "steps" << std::endl;
  //std::cout << "planning time " << duration / 1e6 << "s" << std::endl;

  if (false){
    for (uint i=0; i<tp.path.d0; i+=2){
      std::cout << tp.time(i) << std::endl;
      TP.query(tp.path[i], tp.time(i));
      TP.C.watch(true);
    }
  }
  nn_time_us = startTree.comp_time_us + goalTree.comp_time_us;
  return tp;
}

TimedPath PathFinder_RRT_Time::plan(const arr &q0, const double t0, const arr &qT, const double tGoalLowerBound, const double tGoalUpperBound){
  // construct 'discrete' sampler
  TimedGoalSampler gs = [&](const double t, arr &q) { q=qT;};
  return plan(q0, t0, gs, tGoalLowerBound, tGoalUpperBound);
}


std::vector<std::pair<int, int>> PathFinder_SIRRT_Time::get_safe_intervals(const Eigen::VectorXd &qq) {
  std::vector<std::pair<int, int>> safeIntervals;
  bool isSafe = false;
  int intervalStartFrame = 0;
  arr q(qq.size());
  for (int i=0; i<qq.size(); i++)
  {
    q(i) = qq(i);
  }

  std::vector<std::pair<double, double>> raw_safe_intervals =  this->TP.get_safe_intervals(q);
  for (int safe_int_id=0;safe_int_id<raw_safe_intervals.size();safe_int_id++){
     safeIntervals.emplace_back((raw_safe_intervals[safe_int_id].first-this->t_start)/this->dt,(raw_safe_intervals[safe_int_id].second-this->t_start)/this->dt);
  }
  std::cout << "get safe interval: " ;
  for (int i=0; i<safeIntervals.size(); i++)
    std::cout << safeIntervals[i].first << ' ' << safeIntervals[i].second << ", ";
  std::cout << std::endl;
  return safeIntervals;
}

std::vector<std::pair<int, int>> PathFinder_SIRRT_Time::get_safe_intervals_naive(const Eigen::VectorXd &qq) {
  std::vector<std::pair<int, int>> safeIntervals;
  bool isSafe = false;
  int intervalStartFrame = 0;
  arr q(qq.size());
  for (int i=0; i<qq.size(); i++)
  {
    q(i) = qq(i);
  }

  for (int frame = 0; frame < n_frames; ++frame) {
      double t = this->t_start + frame * dt;
      bool currentState = TP.query(q, t)->isFeasible;

      if (currentState && !isSafe) {
          intervalStartFrame = frame;
          isSafe = true;
      } else if (!currentState && isSafe) {
          safeIntervals.push_back({intervalStartFrame, frame - 1});
          isSafe = false;
      }
  }
  if (isSafe) {
      safeIntervals.push_back({intervalStartFrame, n_frames -1 });
  }
  std::cout << "get safe interval: " ;
  for (int i=0; i<safeIntervals.size(); i++)
    std::cout << safeIntervals[i].first << ' ' << safeIntervals[i].second << ", ";
  std::cout << std::endl;
  return safeIntervals;
}

arr PathFinder_SIRRT_Time::getDelta(const arr &p1, const arr &p2){
  #ifdef PERIODIC
    arr tmp(delta_buffer.N);
    for (uint i=0; i<delta_buffer.N; ++i){
      if (periodicDimensions[i]){
        const double start = p1.p[i];
        const double end = p2.p[i];
        tmp.p[i] = std::fmod(start - end + 3.*RAI_PI, 2*RAI_PI) - RAI_PI;
      }
      else{
        tmp.p[i] = p1.p[i] - p2.p[i];
      }
    }
  #else
    const arr delta = p1 - p2;
  #endif
  
    return tmp;
  }
  
double PathFinder_SIRRT_Time::q_metric(const arr& d) const{
  #if 1
    double dist = 0;
    for (auto *j: TP.C.activeJoints){
      double tmp = 0;
      for (uint i=j->qIndex; i < j->qIndex+j->dim; ++i){
        tmp = std::max({tmp, abs(d(i))});
      }
      dist = std::max({tmp, dist});
    }
    return dist;
  #else
    return length(d);
  #endif
  }

bool PathFinder_SIRRT_Time::extend(Eigen::VectorXd &coords_of_new)
{
  Vertex *q_nearest = this->get_nearest_node(coords_of_new); 
  Eigen::VectorXd delta_vector = coords_of_new - q_nearest->coords;
  double delta = delta_vector.norm();
  if (delta <= (1 / 100000000)) {return false;}

  if (delta < this->planner_range) { return true; }
  coords_of_new = q_nearest->coords + delta_vector.normalized() * this->planner_range;
  return true;
}

std::vector<Vertex *> PathFinder_SIRRT_Time::set_parent(Eigen::VectorXd &coord_rand, std::vector<std::pair<int, int>> &safe_intervals_of_coord_rand)
{   
    std::vector<std::pair<Vertex *, int>> nearest_nodes = this->get_nearest_node_by_radius(coord_rand, 9 * this->planner_range * this->planner_range, (this->current_tree));
    std::vector<Vertex *> added_vertices;
    int safe_interval_ind = 0;
    if (this->current_tree == this->start_tree)
    {
      // std::cout << "start tree\n";
      double time_to_goal = (coord_rand - this->goal_coords).norm() / this->dt / this->vmax; //changed fps here
        std::sort(nearest_nodes.begin(), nearest_nodes.end(), [](const std::pair<Vertex *, int> &a, std::pair<Vertex *, int> &b)
                  { return a.first->arrival_time < b.first->arrival_time; });
        int interval_cnt=0;
        std::cout << "start tree " << safe_intervals_of_coord_rand.size()<<std::endl;

        for (std::pair<int, int> safe_int : safe_intervals_of_coord_rand) // For each interval
        {
          interval_cnt++;
            // if we can't reach goal from that interval - skip
            if ((safe_int.first + time_to_goal) >= this->n_frames)
            {   std::cout<<"if we can't reach goal from that interval - skip" <<std::endl;
                safe_interval_ind++;
                continue;
            }
            bool found_parent = false;

            // for each parent
            int vertex_cnt = 0;
          
            for (std::pair<Vertex *, int> candidate_node : nearest_nodes)
            {
              vertex_cnt++;
              std::cout << "vertex " << vertex_cnt << std::endl;

                double time_to_node = (coord_rand - candidate_node.first->coords).norm() / this->dt / this->vmax;

                // candidate nodes are sorted by ascending arrival time. If arrival time > safe int bound + time to node - > break;
                if (candidate_node.first->arrival_time + time_to_node > safe_int.second)
                {
                  std::cout << candidate_node.first->coords  << std::endl;
                  std::cout <<  candidate_node.first->arrival_time <<" "<< time_to_node << " " << safe_int.second  << std::endl;
                  std::cout << "arrival_time + time_to_node > safe_int "  << std::endl;

                    break;
                }

                // if intervals don't overlap - skip
                if ((candidate_node.first->safe_interval.second + time_to_node < safe_int.first))
                {
                  std::cout << "don't overlap "  << std::endl;

                    continue;
                }

                Eigen::VectorXd start_coords = candidate_node.first->coords;

                // for each departure time
                std::cout << std::max(candidate_node.first->arrival_time, (double)safe_int.first - time_to_node) << ' ' << std::min((double)candidate_node.first->safe_interval.second, (double)safe_int.second - time_to_node)  << std::endl;
                for (double departure_time = std::max(candidate_node.first->arrival_time, (double)safe_int.first - time_to_node); departure_time <= std::min((double)candidate_node.first->safe_interval.second, (double)safe_int.second - time_to_node); departure_time += 1)
                {
                    double arrival_time = departure_time + time_to_node;
                    // fix rounding errors
                    if (arrival_time > safe_intervals_of_coord_rand[safe_interval_ind].second)
                    {
                        //TODO: Check, if this breaks everything
                        departure_time -= arrival_time - safe_intervals_of_coord_rand[safe_interval_ind].second;
                        if(departure_time <candidate_node.first->arrival_time){
                            std::cout << "departure_time <candidate_node.first->arrival_time "  << std::endl;

                            continue;
                        }
                        arrival_time =  safe_intervals_of_coord_rand[safe_interval_ind].second;
                    }
                    if (arrival_time < safe_intervals_of_coord_rand[safe_interval_ind].first)
                    {
                        // departure_time += safe_intervals_of_coord_rand[safe_interval_ind].first - arrival_time;
                        arrival_time = safe_intervals_of_coord_rand[safe_interval_ind].first;

                    }

                    double lv = departure_time*this->dt+this->t_start;
                    double lvv = arrival_time*this->dt+this->t_start;
                    
                    if (!is_collision_motion(start_coords, coord_rand, lv, lvv))
                    {
                      std::cout << "found parent"  << std::endl;

                        this->current_tree->add_vertex(Vertex(coord_rand, safe_intervals_of_coord_rand[safe_interval_ind]), candidate_node.first, departure_time, arrival_time);
                        added_vertices.push_back(&(this->current_tree->array_of_vertices.back()));
                        found_parent = true;
                        break;
                    }
                }
                if (found_parent)
                {
                    break;
                }
            }

            safe_interval_ind++;
        }
    }
    
    else if (this->current_tree == this->goal_tree)
    {
      std::cout << "goal tree \n";
        double time_to_start = (coord_rand - this->root_node->coords).norm() * (1.0/dt) / this->vmax;
        std::sort(nearest_nodes.begin(), nearest_nodes.end(), [](const std::pair<Vertex *, int> &a, std::pair<Vertex *, int> &b)
                  { return a.first->arrival_time > b.first->arrival_time; });
        int interval_cnt1 = 0;
        for (std::pair<int, int> safe_int : safe_intervals_of_coord_rand) // For each interval
        {
          interval_cnt1++;
            // if we can't reach start from that interval - skip
            if ((safe_int.second - time_to_start) < 0)
            {
                safe_interval_ind++;
                continue;
            }
            bool found_parent = false;

            // for each parent
            int near_node_cnt=0;
            for (std::pair<Vertex *, int> candidate_node : nearest_nodes)
            {
              near_node_cnt++;
              // std::cout << "near_node_cnt " << near_node_cnt << " / "  << nearest_nodes.size() << std::endl;

                double time_to_node = (coord_rand - candidate_node.first->coords).norm() / this->dt / this->vmax;

                // candidate nodes are sorted by descending arrival time. If arrival time of parent < safe int bound + time to node - > break;
                if (candidate_node.first->arrival_time - time_to_node < safe_int.first)
                {
                    break;
                }

                // if intervals don't overlap - skip
                if ((candidate_node.first->safe_interval.first - time_to_node > safe_int.second))
                {
                    continue;
                }
                // assert(candidate_node.first);

                Eigen::VectorXd start_coords = candidate_node.first->coords;

                // for each departure time
                std::cout << std::min(candidate_node.first->arrival_time, (double)safe_int.second + time_to_node) << ' ' << std::max((double)candidate_node.first->safe_interval.first, (double)safe_int.first + time_to_node) << std::endl;
                for (double departure_time = std::min(candidate_node.first->arrival_time, (double)safe_int.second + time_to_node); departure_time >= std::max((double)candidate_node.first->safe_interval.first, (double)safe_int.first + time_to_node); departure_time -= 1)
                {
                    // assert(!is_collision_motion(candidate_node.first->coords, candidate_node.first->coords, candidate_node.first->arrival_time, departure_time));
                    double arrival_time = departure_time - time_to_node; // we travel to the past
                    // fix rounding errors
                    
                    if (arrival_time > safe_intervals_of_coord_rand[safe_interval_ind].second)
                    {
                        // departure_time -= arrival_time - safe_intervals_of_coord_rand[safe_interval_ind].second;
                        arrival_time = safe_intervals_of_coord_rand[safe_interval_ind].second;
                    }
                    if (arrival_time < safe_intervals_of_coord_rand[safe_interval_ind].first)
                    {
                        departure_time += safe_intervals_of_coord_rand[safe_interval_ind].first - arrival_time;
                        if(departure_time >candidate_node.first->arrival_time){
                            continue;
                        }
                        arrival_time = safe_intervals_of_coord_rand[safe_interval_ind].first;
                    }

                    double lv = arrival_time*this->dt+this->t_start;
                    double lvv = departure_time*this->dt+this->t_start;

                    if (!is_collision_motion(coord_rand, start_coords, lv, lvv))
                    {
                        this->current_tree->add_vertex(Vertex(coord_rand, safe_intervals_of_coord_rand[safe_interval_ind]), candidate_node.first, departure_time, arrival_time);
                        added_vertices.push_back(&(this->current_tree->array_of_vertices.back()));
                        found_parent = true;
                        break;
                    }
                }
                if (found_parent)
                {
                    break;
                }
            }
            safe_interval_ind++;
        }
    }
    return added_vertices;
}


void PathFinder_SIRRT_Time::swap_trees()
{
    std::swap(this->current_tree, this->other_tree);
}

std::vector<Vertex *> PathFinder_SIRRT_Time::grow_tree(Eigen::VectorXd &coord_rand, std::vector<std::pair<int, int>> &safe_intervals_of_coord_rand)
{
    std::vector<Vertex *> result;
    std::vector<Vertex *> array_of_new_nodes = this->set_parent(coord_rand, safe_intervals_of_coord_rand);

    if (array_of_new_nodes.size() == 0)
    {
        return result;
    }
    result.reserve(result.size() + std::distance(array_of_new_nodes.begin(), array_of_new_nodes.end()));
    result.insert(result.end(), array_of_new_nodes.begin(), array_of_new_nodes.end());

    // while (array_of_new_nodes.size() != 0)
    // {
    //     Vertex *node = array_of_new_nodes.back();
    //     array_of_new_nodes.pop_back();
    // }
    return result;
}


bool PathFinder_SIRRT_Time::connect_trees(Eigen::VectorXd &coord_rand, std::vector<std::pair<int, int>> &safe_intervals_of_coord_rand, std::vector<Vertex *> another_tree_new_nodes)
{
  
    const Eigen::VectorXd original_coord_rand = coord_rand;
    while (true)
    {
        coord_rand = original_coord_rand;
        this->swap_trees();
        bool is_ok = this->extend(coord_rand);
        this->swap_trees();
        if (!is_ok)
        {
            break;
        }


        // for (int i = 0; i < this->dimensionality; i++) {
        //   robot_q(i) = coord_rand(i);
        // }


        std::vector<std::pair<int, int>> safe_intervals_of_coord_rand;

        safe_intervals_of_coord_rand = this->get_safe_intervals(coord_rand);
        this->swap_trees();
        std::vector<Vertex *> new_nodes = this->grow_tree(coord_rand, safe_intervals_of_coord_rand);
        this->swap_trees();
        if (new_nodes.size() == 0)
        {
            break;
        }

        if (original_coord_rand == coord_rand) // if we reached original nodes of first tree
        {
            for (Vertex *node : new_nodes) // for new node in second tree
            {
                for (Vertex *another_tree_node : another_tree_new_nodes) // for new node in first tree
                {
                    if (node->coords == another_tree_node->coords)
                    {
                        if (node->safe_interval == another_tree_node->safe_interval)
                        {
                            if (node->tree_id == 1) // TODO: change tree_id to enums. 1==goal_Tree
                            {
                                if (node->arrival_time >= another_tree_node->arrival_time)
                                {
                                    this->goal_reached = true;
                                    this->goal_nodes = std::pair<Vertex *, Vertex *>(another_tree_node, node);
                                    this->prune_goal_tree();
                                }
                            }
                            else
                            {
                                if (node->arrival_time <= another_tree_node->arrival_time)
                                {
                                    this->goal_reached = true;
                                    this->goal_nodes = std::pair<Vertex *, Vertex *>(node, another_tree_node);
                                    this->prune_goal_tree();
                                }
                            }
                            return true;
                        }
                    }
                }
                if (this->goal_reached)
                {
                    break;
                }
            }
            break;
        }
        coord_rand = original_coord_rand;
    }

    return false;
}


void PathFinder_SIRRT_Time::prune_goal_tree()
{
    assert(this->goal_reached);

    Vertex *start_tree_node = this->goal_nodes.first;
    Vertex *goal_tree_node_child = this->goal_nodes.second;

    assert(goal_tree_node_child->parent); // Must have parent, because RRTConnect can't just connect to root node
    Vertex *goal_tree_node = goal_tree_node_child->parent; 
    
    while (goal_tree_node)
    {
        double time_to_node = (start_tree_node->coords - goal_tree_node->coords).norm() * (1.0/dt) / this->vmax;

        auto start_coords = start_tree_node->coords;
        bool found_dep_time = false;
        // for each departure time
        for (double departure_time = std::max(start_tree_node->arrival_time, (double)goal_tree_node->safe_interval.first - time_to_node); departure_time <= goal_tree_node_child->arrival_time; departure_time += 1)
        {
            
            // assert(!is_collision_motion(start_coords, start_coords, start_tree_node->arrival_time, departure_time));
            double arrival_time = departure_time + time_to_node;
            double lv = departure_time*this->dt+this->t_start;
            double lvv = arrival_time*this->dt+this->t_start;
            if (!is_collision_motion(start_coords, goal_tree_node->coords, lv, lvv))
            {
                this->start_tree->add_vertex(Vertex(goal_tree_node->coords, goal_tree_node->safe_interval), start_tree_node, departure_time, arrival_time);
                start_tree_node = &(this->start_tree->array_of_vertices.back());
                found_dep_time = true;
                break;
            }
        }
        if(!found_dep_time)
        {
            this->start_tree->add_vertex(Vertex(goal_tree_node->coords, goal_tree_node->safe_interval), start_tree_node, goal_tree_node_child->arrival_time, goal_tree_node_child->departure_from_parent_time);
            start_tree_node = &(this->start_tree->array_of_vertices.back());
        }
        goal_tree_node_child = goal_tree_node;
        goal_tree_node = goal_tree_node->parent;

    }
    this->finish_node = start_tree_node;
}

bool PathFinder_SIRRT_Time::check_planner_termination_condition() const
{
    if (this->stop_when_path_found)
    {
        return !this->goal_reached && std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - this->solver_start_time).count() < this->max_planning_time;
    }
    return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - this->solver_start_time).count() < this->max_planning_time;
}

TimedPath PathFinder_SIRRT_Time::plan(const arr &q0, const double &t0, const arr &q_goal, const double &t_up){
  this->TP.min_time = t0;
  std::cout << "t0 " << t0 << std::endl;
  std::cout << "t_up " << t_up << std::endl;

  this->TP.max_time = t_up;  
  this->TP.init_safe_interval_collisison_check(q0);
  
  std::cout << "РАЗМЕРНОСТЬ " << q0.N << std::endl;
  const bool fixedTime = rai::getParameter<bool>("assembly/fixedTime", false); 
    this->dimensionality = q0.N;
    this->t_start = t0;
    this->t_max = t_up;
    this->n_frames = std::ceil((t_up - t0) / dt);
    this->goal_reached = false;
    this->max_planning_time = 15;
    goal_nodes = std::pair<Vertex *, Vertex *>(nullptr,nullptr);

    this->start_tree = new Tree_nf("start_tree", 0, this->dimensionality);
    this->goal_tree = new  Tree_nf("goal_tree", 1, this->dimensionality);

    this->start_tree->array_of_vertices.reserve(1000000);  // DONT REMOVE IT!
    this->goal_tree->array_of_vertices.reserve(1000000);   // DONT REMOVE IT!
    std::cout << "start:" <<  q0    << " time: " << t0 <<std::endl;
    std::cout << "goal:"  << q_goal << " time: " << t_up <<std::endl;
    Eigen::VectorXd q0_coords(this->dimensionality);
    Eigen::VectorXd q_goal_coords(this->dimensionality);
    for (int i = 0; i < q0.N; ++i)
    {
      q0_coords(i) = q0(i);
      q_goal_coords(i) = q_goal(i);

    }

    std::vector<std::pair<int, int>> start_safe_intervals = this->get_safe_intervals(q0_coords);

    this->start_tree->add_vertex(Vertex(q0_coords, start_safe_intervals[0]), nullptr, -1, 0);

    this->root_node = &this->start_tree->array_of_vertices[0];
    // this->root_node->arrival_time = t0;
    
    goal_coords.resize(this->dimensionality);
    for (int i = 0; i < this->dimensionality; i++) {
      this->goal_coords(i) = q_goal(i);
    }

    this->goal_safe_intervals = this->get_safe_intervals(q_goal_coords);

    this->goal_tree->add_vertex(Vertex(q_goal_coords, this->goal_safe_intervals.back()), nullptr, -1, this->goal_safe_intervals.back().second);
    this->current_tree = this->goal_tree;
    this->other_tree = this->start_tree;
    
  
    //auto start_rrt_time = std::chrono::high_resolution_clock::now();
  
    if (!TP.query(q0, t0, 0)->isFeasible){
      spdlog::error("Initial point not feasible");
      TP.C.watch(true);
      return TimedPath({}, {});
    }
 
    rai::Configuration DSP;
    DSP.copy(TP.C);
    // DSP.gl()->add(this->start_tree);
    // DSP.gl()->add(goalTree);
    // if (disp){
    //   DSP.copy(TP.C);
    //   DSP.gl()->add(startTree);
    //   DSP.gl()->add(goalTree);
  
    //   DSP.watch(false);
    // }

  
    TimedPath tp({}, {});

    this->solver_start_time = std::chrono::steady_clock::now();
    int iter = 0;

    // spdlog::info("МЫ ЗДЕСЬ 0");
    // int t = min_time + 1. * rnd(0) * (max_time - min_time) * i / 10.;
    // int max_goal_time 0;
    //max_goal_time = std::max({max_goal_time, t});

    int v_count = 0;
    const double min_l = q_metric(getDelta(q_goal, q0));


    while (this->check_planner_termination_condition() && !this->goal_reached)
    {        
        v_count++;
        std::cout<< v_count << " " << this->start_tree->array_of_vertices.size() << " " << this->goal_tree->array_of_vertices.size() <<  std::endl;
        Eigen::VectorXd coord_rand(this->dimensionality);
  
        // arr qs = TP.sample(q0, q_goal, (t_up - t0) * vmax, min_l);
        arr qs = TP.sample(); //q0, q_goal, (max_goal_time - t0) * vmax, min_l);
        for (int i = 0; i < this->dimensionality; i++) {
          coord_rand(i) = qs(i);
        }
        std::cout << qs << std::endl;

        bool is_ok = this->extend(coord_rand);
        if (!is_ok)
        {
            continue;
        }
        std::vector<std::pair<int, int>> safe_intervals_of_coord_rand;

        safe_intervals_of_coord_rand = this->get_safe_intervals(coord_rand);

        std::vector<Vertex *> new_nodes = this->grow_tree(coord_rand, safe_intervals_of_coord_rand);

        if (new_nodes.size() == 0)
        {
            this->swap_trees();
            continue;
        }
        bool connected = connect_trees(coord_rand, safe_intervals_of_coord_rand, new_nodes);
        this->swap_trees();
    }

    // spdlog::info("МЫ ЗДЕСЬ 1");
    if (this->goal_reached) std::cout << "МЫ НАШЛИ РЕШЕНИЕ, МЫ НАШЛИ ЕГО\n";

    if (!this->goal_reached) {delete this->start_tree; delete this->goal_tree;      start_tree = nullptr;
      goal_tree = nullptr;
      current_tree = nullptr;
      other_tree = nullptr;
      root_node = nullptr; return TimedPath({}, {});}
    else
    {
      arr path;
      arr time;
      path.resize(0, this->dimensionality);

      std::vector<Vertex *> result;
      Vertex *current = this->finish_node;  
      while (current)
      {
          result.push_back(current);
          current = current->parent;
      }
      std::reverse(result.begin(), result.end());
      std::cout << "NODES IN SOLUTION: " << result.size() << "\n";
      
      // for (Vertex* v:result)
      // {
      //   std::cout << v->q << std::endl;
      //   std::cout << "eigen = ";
      //   for (auto angle:v->coords)
      //     std::cout << angle << ' ';

      //     std::cout << std::endl;
      // }
       std::cout << "====================================================\n"; 
      int v_count = 0;
      for (Vertex* v:result)
      {
        double t = 0.0;
        arr v_q(v->coords.size());
        for (int i = 0; i<v->coords.size(); i++)
          v_q(i) = v->coords(i);
        path.append(v_q);
        if (v_count == 0)
        {
          time.append(this->t_start);
          t = this->t_start;
        }
        else
        {
          time.append(v->departure_from_parent_time + (v->coords - v->parent->coords).norm() / this->vmax +this->t_start);
          t = v->departure_from_parent_time + (v->coords - v->parent->coords).norm() / this->vmax +this->t_start;
        }
        std::cout << "time: " << t << ", " << v_q << std::endl;
        v_count++;
      }
      delete this->start_tree;
      delete this->goal_tree;
      start_tree = nullptr;
      goal_tree = nullptr;
      current_tree = nullptr;
      other_tree = nullptr;
      root_node = nullptr;
      // std::cout << time << std::endl;
      // std::cout << path.d <<' ' << path.d0 << ' ' << path.d1 << ' ' << path.d2 << std::endl;
      // std::cout << path << std::endl;
      return TimedPath(path, time);
    }
  }
  
// TimedPath PathFinder_SIRRT_Time::plan(const arr &q0, const double t0, const arr &qT, const double tGoalLowerBound, const double tGoalUpperBound){
//     // construct 'discrete' sampler
//     TimedGoalSampler gs = [&](const double t, arr &q) { q=qT;};
//     return plan(q0, t0, gs, tGoalLowerBound, tGoalUpperBound);
//   }
  

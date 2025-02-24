#include "rrt-time.h"

#include <chrono>

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
      const auto duration =
          std::chrono::duration_cast<std::chrono::microseconds>(end_time -
                                                                start_time)
              .count();

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
  while(true){
    if (disp){
      DSP.watch(false);
    }

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
    if (r(0) > goalSampleProbability){// && goals.N < 10){
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


//=========================================================================================================================================================================================
//=========================================================================================================================================================================================


arr PathFinder_SIRRT_Time::getDelta(const arr &p1, const arr &p2){
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
  
double PathFinder_SIRRT_Time::q_metric(const arr& d) const{
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
  
double PathFinder_SIRRT_Time::distance(const Node &n1, const Node &n2){
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
  
Vertex* PathFinder_SIRRT_Time::steer(const Node &start, const Node &goal, bool reverse){
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
  
arr PathFinder_SIRRT_Time::projectToManifold(const arr &q, const double t){
    if (prePlannedFrames.N == 0){
      return q;
    }
  
    if (t > tPrePlanned){
      return q;
    }
      TP.query(q, t);
    arr qNew = TP.C.getJointState();
    if (false){
  
      TP.C.setJointState(q);
      TP.C.watch(true);
  
      TP.C.setJointState(qNew);
      TP.C.watch(true);
    }
    return qNew * 1.;
  }
  
Vertex* PathFinder_RRT_Time::extend(Tree* tree, const Node &goal, const bool connect){
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
        const auto duration =
            std::chrono::duration_cast<std::chrono::microseconds>(end_time -
                                                                  start_time)
                .count();
  
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
  
TimedPath PathFinder_SIRRT_Time::plan(const arr &q0, const double t0, const TimedGoalSampler gs, double tGoalLowerBound, double tGoalUpperBound){
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
    while(true){
      if (disp){
        DSP.watch(false);
      }
  
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
      if (r(0) > goalSampleProbability){// && goals.N < 10){
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
  
TimedPath PathFinder_SIRRT_Time::plan(const arr &q0, const double t0, const arr &qT, const double tGoalLowerBound, const double tGoalUpperBound){
    // construct 'discrete' sampler
    TimedGoalSampler gs = [&](const double t, arr &q) { q=qT;};
    return plan(q0, t0, gs, tGoalLowerBound, tGoalUpperBound);
  }
  

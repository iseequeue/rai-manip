#pragma once

#include "spdlog/spdlog.h"

#include <algorithm>

#include <Gui/opengl.h>
#include <GL/gl.h>

#include <PlanningSubroutines/ConfigurationProblem.h>
#include "timedPath.h"

#include <nanoflann.hpp>
#include <Eigen/Core>

struct Node{
  arr q;
  arr x;
  double t;
  double cost = -1;

  Node *parent = nullptr;
  std::vector<Node*> children;

  Node(arr _q, double _t): q(_q), t(_t) {};

  void addChild(Node* node){
    children.push_back(node);
  }

  void print() const{
    std::cout << "t: " << t << std::endl;
    std::cout << "q: " << q << std::endl;
  }

  void updateCost(const double newCost){
    const double diff = cost - newCost;

    for (auto c: children){
      c->updateCost(c->cost - diff);
    }
    cost = newCost;
  }
};

typedef Eigen::Matrix<double, 6, 1> VertexCoordType;

struct Vertex{
  arr q;
  
  VertexCoordType coords; 
  std::pair<double, double> safe_interval;
  double cost = -1;

  Vertex *parent = nullptr;
  std::vector<Vertex*> children{};

  double arrival_time = -1;
  double before_rewiring_arrival_time = -1;
  double departure_from_parent_time = -1;
  int tree_id = -1;

  Vertex(arr _q, const std::pair<double, double> &_safe_interval): 
  q(_q), safe_interval(_safe_interval) {
    for (int i = 0; i < q.N; i++) {
      coords(i) = q(i);
  }
  };

  void addChild(Vertex* v){
    children.push_back(v);
  }

  void print() const{
    std::cout << "t: [" << safe_interval.first << " " << safe_interval.second <<"] \n";
    std::cout << "q: " << q << std::endl;
  }

  void updateCost(const double newCost){
    const double diff = cost - newCost;

    for (auto c: children){
      c->updateCost(c->cost - diff);
    }
    cost = newCost;
  }
};

#define TREE_DIMENSIONALITY 6
typedef nanoflann::KDTreeSingleIndexDynamicAdaptor<nanoflann::L2_Simple_Adaptor<double, Tree_nf>, Tree_nf, TREE_DIMENSIONALITY> KdTree;

struct Tree_nf: GLDrawer
{
    Tree_nf(): tree_name(""), tree_idx(-1), kd_tree(TREE_DIMENSIONALITY, *this, nanoflann::KDTreeSingleIndexAdaptorParams(25)) {};;
    Tree_nf(const std::string &tree_name_, size_t tree_idx_) : tree_name(tree_name_), tree_idx(tree_idx_),
    kd_tree(TREE_DIMENSIONALITY, *this, nanoflann::KDTreeSingleIndexAdaptorParams(25)) {};;
    ~Tree_nf(){};                                     // destructor
    Tree_nf(const Tree_nf &other) = delete;            // copy constructor
    Tree_nf(Tree_nf &&other) = default;                // move constructor
    Tree_nf &operator=(const Tree_nf &other) = delete; // copy assignment
    Tree_nf &operator=(Tree_nf &&other) = default;     // move assignment

    std::string tree_name;
    size_t tree_idx;

    std::vector<Vertex> array_of_vertices;
    KdTree kd_tree;

    Vertex &getNearestState(const Vertex q){
      const size_t num_results{1};
      size_t q_near_idx{0};
      double out_dist_sqr{0};
      nanoflann::KNNResultSet<double> result_set(num_results);
      result_set.init(&q_near_idx, &out_dist_sqr);
      this->kd_tree.findNeighbors(result_set, q.coords.data(), {0});
      return array_of_vertices[q_near_idx];
  }

    void add_vertex(Vertex q_new, Vertex *q_parent, double departure_time, double arrival_time)
    {
      q_new.tree_id = this->tree_idx;
      q_new.parent = q_parent;
      q_new.arrival_time = arrival_time;
      q_new.departure_from_parent_time = departure_time;
      this->array_of_vertices.push_back(q_new);
  
      size_t N{array_of_vertices.size() - 1};
        this->kd_tree.addPoints(N, N);
  
      if (q_parent != nullptr)
      {
          q_parent->children.push_back(&(this->array_of_vertices.back()));
      }
  }
    void delete_vertex(int vertex_id){   this->kd_tree.removePoint(vertex_id); }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX & /* bb */) const { return false; } 
    inline size_t kdtree_get_point_count() const { return array_of_vertices.size(); }
    inline double kdtree_get_pt(const size_t idx, const size_t dim) const { return array_of_vertices[idx].q(dim); }  // у них индексация через круглые скобки
};


struct Tree: GLDrawer {
  double comp_time_us{0.};

  // TODO: replace with multiset, use orderedness to discard 
  //       nodes that can not be close
  // OR:   implement kd-tree that workd with our weird distance function
  std::vector<Node*> nodes;
  bool reverse;

  std::function<double (const Node&, const Node&)> distance;

  Tree(std::function<double (const Node&, const Node&)> _distance, bool _reverse=false):reverse(_reverse), distance(_distance){};
  ~Tree(){
    //std::cout << "deleting stuff" << std::endl;
    for (const auto &n: nodes){
      delete n;
    }
  }

  Node * getNearest(const Node &target, const bool print = false){
    const auto start_time = std::chrono::high_resolution_clock::now();

    Node *close = nullptr;
    double min_dist = std::numeric_limits<double>::max();
    double d;
    for (Node* n: nodes){
      if (reverse) {d = distance(target, *n);}
      else {d = distance(*n, target);}

      if (print){
        target.print();
        n->print();
        std::cout << "dist: " << d << std::endl << std::endl;
      }

      if (d < min_dist){
        close = n;
        min_dist = d;
      }
    }

    if (print){
      std::cout << std::endl;
    }

    const auto end_time = std::chrono::high_resolution_clock::now();
    const auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(end_time -
                                                              start_time)
            .count();

    comp_time_us += duration;

    return close;
  };

  std::vector<Node*> getNear(const Node &start, const double radius){
    const auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<Node*> closeNodes;
    double d;
    for (Node* n: nodes){
      if (n == &start){
        continue;
      }

      if (!reverse) {d = distance(start, *n);}
      else {d = distance(*n, start);}

      if (d < radius){
        closeNodes.push_back(n);
      }
    }

    const auto end_time = std::chrono::high_resolution_clock::now();
    const auto duration =
        std::chrono::duration_cast<std::chrono::microseconds>(end_time -
                                                              start_time)
            .count();

    comp_time_us += duration;

    return closeNodes;
  }

  bool rewire(Node *nNode, TimedConfigurationProblem &TP){
    bool rewired = false;
    auto closeNodes = getNear(*nNode, 1);
    for (Node *n: closeNodes){
      //std::cout << reverse << std::endl;
      //std::cout << nNode->q << " " << nNode->t << std::endl;
      //std::cout << n->q << " " << n->t << std::endl;

      if (n->parent == nNode){
        continue;
      }

      double d;
      if (!reverse) {d = distance(*nNode, *n);}
      else {d = distance(*n, *nNode);}

      //std::cout << d << std::endl;

      if (n->cost > nNode->cost + d){
        //std::cout << n->cost << " " << nNode->cost + d << std::endl;
        if (!TP.checkEdge(nNode->q, n->t, n->q, n->t)) {continue;}

        n->updateCost(nNode->cost + d);

        std::vector<Node*> &childVec = n->parent->children;
        childVec.erase(std::remove(childVec.begin(), childVec.end(), n), childVec.end());
        n->parent = nNode;
        nNode->addChild(n);

        rewired = true;
      }
    }
    return rewired;
  }

  void addNode(Node* n){
    nodes.push_back(n);
  }

  void glDraw(OpenGL &gl){
    if (reverse){
      glColor(0., 0., 0.);
    }
    else{
      glColor(1., 1., 1.);
    }
    glLineWidth(2.f);
    glBegin(GL_LINES);
    for(uint i=0;i<nodes.size();i++){
      if (nodes[i]->parent){
        glVertex3dv(&(nodes[i]->parent->x(0)));
        glVertex3dv(&(nodes[i]->x(0)));
      }
    }
    glEnd();
    glLineWidth(1.f);
  }
};

enum class SamplingType{
  BOX_CONSTRAINED_CONDITIONAL_SAMPLING,
  BOX_CONSTRAINED_REJECTION_SAMPLING,

  CONDITIONAL_SAMPLING, // TODO, not implemented
  REJECTION_SAMPLING, // TODO, not implemented
};

using TimedGoalSampler = std::function<void (const double, arr&)>;

struct PathFinder_RRT_Time{
  TimedConfigurationProblem &TP;

  double vmax = .1;
  double lambda = .9;
  double step_time = 1;

  double goalSampleProbability = 0.9;//0.9

  uint maxInitialSamples = 50;

  uint maxIter = 1000;
  bool verbose = false;
  bool disp = false;

  double tol = 1e-3;
  bool connect = true;
  bool optimize = false;
  uint conv_iter = 500;
  bool sampleAdditionalNodes = false;

  SamplingType sampling_type = SamplingType::BOX_CONSTRAINED_CONDITIONAL_SAMPLING;
  // SamplingType sampling_type = SamplingType::REJECTION_SAMPLING;

  bool informed_sampling = false;

  std::vector<bool> periodicDimensions;

  double edge_checking_time_us{0.};
  double nn_time_us{0.};

  FrameL prePlannedFrames;
  uint tPrePlanned;

  PathFinder_RRT_Time(TimedConfigurationProblem &_TP) :TP(_TP){
    delta_buffer = arr(TP.C.getJointState().N);
    periodicDimensions = std::vector<bool>(TP.C.getJointState().N, false);

    for (auto *j: TP.C.activeJoints){
      if (j->type == rai::JT_hingeX || j->type == rai::JT_hingeY|| j->type == rai::JT_hingeZ){
        periodicDimensions[j->qIndex] = true;
      }
    }
  };

  arr projectToManifold(const arr &q, const double t);

  TimedPath plan(const arr &q0, const double t0, const arr &qT, double tGoalLowerBound=0, double tGoalUpperBound=-1);
  TimedPath plan(const arr &q0, const double t0, const TimedGoalSampler gs, double tGoalLowerBound=0, double tGoalUpperBound=-1);

  arr delta_buffer;
  double q_metric(const arr& d) const;
  double distance(const Node &n1, const Node &n2);

  Node* extend(Tree* tree, const Node &goal, const bool connect = false);

  arr getDelta(const arr &p1, const arr &p2);
  Node* steer(const Node &start, const Node &goal, bool reverse);

  TimedPath extractPath(Tree* t1, Tree* t2, Node* leafNode1, Node* leafNode2){
    //std::cout << "cost: " << leafNode1->cost + leafNode2->cost << std::endl;

    // traverse over both sides and export result to time, path
    std::vector<Node*> np;
    {
      Node* n = leafNode1;
      while(!!n->parent){
        np.push_back(n);
        n = n->parent;
      }
      np.push_back(n);

      std::reverse(np.begin(), np.end());
    }

    {
      Node* n = leafNode2;
      while(n->parent){
        np.push_back(n);
        n = n->parent;
      }
      np.push_back(n);
    }

    if (t1->reverse){
      std::reverse(np.begin(), np.end());
    }

    arr path;
    arr time;

    path.resize(0, leafNode1->q.N);
    for (Node* n: np){
      time.append(n->t);
      path.append(n->q);
    }

    return TimedPath(path, time);
  }
};


struct PathFinder_SIRRT_Time{
  // Ихний стафф
  TimedConfigurationProblem &TP;

  FrameL prePlannedFrames;
  uint tPrePlanned;
  arr delta_buffer;
  std::vector<bool> periodicDimensions;

  PathFinder_SIRRT_Time(TimedConfigurationProblem &_TP) :TP(_TP){
    this->n_frames = std::ceil((t_max - t_start) / dt);
    delta_buffer = arr(TP.C.getJointState().N);
    periodicDimensions = std::vector<bool>(TP.C.getJointState().N, false);

    for (auto *j: TP.C.activeJoints){
      if (j->type == rai::JT_hingeX || j->type == rai::JT_hingeY|| j->type == rai::JT_hingeZ){
        periodicDimensions[j->qIndex] = true;
      }
    }
  };


  arr projectToManifold(const arr &q, const double t){
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
  };

  bool verbose = false;
  bool disp = false;

  // Наш стафф
  Tree_nf *start_tree;
  Tree_nf *goal_tree;
  Tree_nf *current_tree;
  Tree_nf *other_tree;
  Vertex *root_node = nullptr;

  int n_frames;
  int dimensionality;

  double t_start;
  double t_max;
  double dt = 1.0/15.0; // HARDCODE!, fps analog

  std::vector<std::pair<int, int>> get_safe_intervals(const arr &q);

  double vmax = .1;
  double goal_bias = 0.4;
  double planner_range = 1.0;


  bool is_collision_motion(const arr &start_coords, const arr &end_coords, double &start_time, double &end_time)
  {
    return TP.checkEdge(start_coords, start_time, end_coords, end_time, 20);
  }
  bool is_collision_state(const arr &q, int &time) { return TP.query(q, time * dt)->isFeasible; }

  bool goal_reached = false;
  Vertex * finish_node;
  
  VertexCoordType goal_coords;

  Vertex *get_nearest_node(const VertexCoordType &coords)
  {
    const size_t num_results = 1;
    size_t ret_index;
    double out_dist_sqr;
    // nanoflann::SearchParams search_params;
    nanoflann::KNNResultSet<double> resultSet(num_results);
    resultSet.init(&ret_index, &out_dist_sqr);
    this->current_tree->kd_tree.findNeighbors(resultSet, coords.data(), {0});
    return &(this->current_tree->array_of_vertices[ret_index]);
  }

  std::vector<std::pair<Vertex *, int>> get_nearest_node_by_radius(VertexCoordType &coords, double raduis, Tree_nf *tree)
  {
    std::vector<std::pair<Vertex *, int>> result;
    std::vector<std::pair<size_t, double>> indices_dists;
    nanoflann::RadiusResultSet<double, size_t> resultSet(raduis, indices_dists);

    tree->kd_tree.findNeighbors(resultSet, coords.data(), {0});

    result.reserve(resultSet.m_indices_dists.size());
    for (auto node_ind_dist_pair : resultSet.m_indices_dists)
    {
        result.emplace_back(&(tree->array_of_vertices[node_ind_dist_pair.first]), node_ind_dist_pair.first);
    }
    return result;
  }

  bool extend(VertexCoordType &coords_of_new);

  std::vector<Vertex *> set_parent(VertexCoordType &coord_rand, std::vector<std::pair<int, int>> &safe_intervals_of_coord_rand);

  
  std::pair<Vertex *, Vertex *> goal_nodes = std::pair<Vertex *, Vertex *>(nullptr,nullptr);
  bool connect_trees(VertexCoordType& coord_rand, std::vector<std::pair<int, int>>& safe_intervals_of_coord_rand,std::vector<Vertex* > new_nodes);
  void swap_trees();
  std::vector<Vertex*> grow_tree(VertexCoordType &coord_rand, std::vector<std::pair<int, int>> &safe_intervals_of_coord_rand);
  bool connect_trees(VertexCoordType &coord_rand, std::vector<std::pair<int, int>> &safe_intervals_of_coord_rand, std::vector<Vertex *> another_tree_new_nodes);
  void prune_goal_tree();

  bool check_planner_termination_condition() const;
  bool stop_when_path_found = true;
  std::chrono::time_point<std::chrono::steady_clock> solver_start_time;
  float max_planning_time = 180;

  //Финальное творение
  TimedPath plan(const arr &q0, const double t0, const arr &q_goal);
  // TimedPath plan(const arr &q0, const double t0, const TimedGoalSampler gs, double tGoalLowerBound=0, double tGoalUpperBound=-1);

};
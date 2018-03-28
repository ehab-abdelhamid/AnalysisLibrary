#include <iostream>
#include <vector>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <cfloat>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <ctime>
#include <mpi.h>
#include <fstream>
#include <sstream>
using namespace std;

#define NDEBUG

int world_rank;
int world_size;

int num_points = 100;
int num_dimensions = 2;
bool print_canopies = true;

float T1 = 200; // loose distance
float T2 = 100; // tight distance

string fileName;

class DataPoint {
private:
  int point_id;
  int cluster_id = -1; // from 0..(k-1), or -1 if unassigned
protected:
  float* vals;
public:
  static int point_id_counter;
  DataPoint() {
    vals = new float[num_dimensions];
    memset(vals, 0, sizeof(vals));
    cluster_id = -1;
    point_id = point_id_counter++;
  }

  DataPoint(float *vals) {
    this->vals = new float[num_dimensions];
    memcpy(this->vals, vals, sizeof(float)*num_dimensions);
    cluster_id = -1;
    point_id = point_id_counter++;
  }

  DataPoint(vector<float>& vals) {
    assert(vals.size() == num_dimensions);
    this->vals = new float[num_dimensions];
    memcpy(this->vals, vals.data(), sizeof(this->vals));
    cluster_id = -1;
    point_id = point_id_counter++;
  }

  int get_point_id() const {
    return point_id;
  }

  int get_cluster_id() const {
    return cluster_id;
  }

  void set_cluster_id(int cluster_id) {
    this->cluster_id = cluster_id;
  }

  float get_val(int i) const {
    //assert(i >= 0 && i < k);
    return vals[i];
  }

  float get_squared_dist(DataPoint& other_point) const {
    float result = 0.0;
    for (int i=0; i<num_dimensions; i++) {
      result += pow(vals[i] - other_point.vals[i], 2);
    }
    return result;
  }

  void print() const {
    for (int i=0; i<num_dimensions; i++) {
      cout << vals[i] << ' ';
    }
    cout << endl;
  }

  bool operator==(const DataPoint& other_point) const {
    return other_point.get_point_id() == point_id;
  }
};

class CSVRow
{
    public:
        float const& operator[](std::size_t index) const
        {
            return m_data[index];
        }
        std::size_t size() const
	{
		return m_data.size();
	}
        void readNextRow(std::istream& str)
	{
		std::string         line;
		std::getline(str, line);

		std::stringstream   lineStream(line);
		std::string         cell;

		m_data.clear();
		while(std::getline(lineStream, cell, ','))
		{
			float cellF = atof(cell.c_str());
			m_data.push_back(cellF);
		}
		// This checks for a trailing comma with no data after it.
		if (!lineStream && cell.empty())
		{
			// If there was a trailing comma then add an empty element.
			m_data.push_back(-1);
		}
	}

    private:
        std::vector<float> m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}

class Canopy {
private:
  vector<DataPoint*> data_points;
  DataPoint* centre;
public:
  Canopy(DataPoint* centre) {
    assert(centre);
    this->centre = centre;
    data_points.push_back(centre);
  }

  void add_point(DataPoint* point) {
    data_points.push_back(point);
  }

  vector<DataPoint*> get_data_points() {
    return data_points;
  }

  void print() {
    centre->print();
  }

  void printElements() const {
    for (int i=0; i<data_points.size(); i++) {
      cout << "\t";
      (*data_points[i]).print();
    }
  }
};

void generate_points(vector<DataPoint*>& points, int num_points) {
  for (int i=0; i<num_points; i++) {
    float vals[num_dimensions];
    for (int j=0; j<num_dimensions; j++) {
      vals[j] = ((rand() / float(RAND_MAX) ) * 1000);
    }
    DataPoint* point = new DataPoint(vals);
    points.push_back(point);
  }
}

void read_points(vector<DataPoint*>& points, string file_name, bool ignore) {
	std::ifstream file(file_name.c_str());

	CSVRow row;
	while(file >> row)
	{
		if(ignore) {
			ignore = false;
			continue;
		}
		num_dimensions=row.size();
		float vals[row.size()];
		for(unsigned int i=0;i<row.size();i++) {
			vals[i]=(float)row[i];
		}
		DataPoint* point = new DataPoint(vals);
		points.push_back(point);
	}
}

vector<Canopy> canopy_mpi(vector<DataPoint*>& all_points) {
  //printf("point_id_counter %d %d\n", DataPoint::point_id_counter, world_rank);

  // the main process will generate points and scatter them among
  // other processes
  int rec_buff_cnt =
    (num_points / world_size + (world_rank < num_points % world_size ? 1 : 0)) * num_dimensions;
  float* rec_buff = new float[rec_buff_cnt];

  vector<DataPoint*> points;

  int* scatter_counts = new int[world_size];
  int* scatter_displs = new int[world_size];
  int sum = 0;
  for (int i=0; i<world_size; i++) {
    scatter_counts[i] =  num_points / world_size * num_dimensions;
    if (i < num_points % world_size) {
      scatter_counts[i] += num_dimensions;
    }
    scatter_displs[i] = sum;
    sum += scatter_counts[i];
  }

  float* data = NULL;
  if (world_rank == 0) {
    data = new float[num_points * num_dimensions];
    int i = 0;
    for (const DataPoint* p : all_points) {
      for (int j=0; j<num_dimensions; j++) {
        data[i++] = p->get_val(j);
      }
    }
  }
  MPI_Scatterv(data, scatter_counts, scatter_displs, MPI_FLOAT,
               rec_buff, rec_buff_cnt, MPI_FLOAT,
               0, MPI_COMM_WORLD);

  // reconstruct points from rec_points_data
  DataPoint::point_id_counter = 0;
  for (int i=0; i<rec_buff_cnt; i += num_dimensions) {
    DataPoint* data_point = new DataPoint(rec_buff + i);
    points.push_back(data_point);
  }

  // each process has a chunk of points

#ifdef DEBUG
  // TODO : remove this
  printf("process %d\n", world_rank);
  for (DataPoint* p : points) {
    printf("%d   ", world_rank);
    p->print();
  }
#endif

  unordered_set<DataPoint*> point_set;
  for (DataPoint* p : points) {
    point_set.insert(p);
  }

  vector<Canopy> canopies;
  
  while (!point_set.empty()) {
    DataPoint* new_canopy_centre = *(point_set.begin());
    Canopy new_canopy(new_canopy_centre);
    point_set.erase(new_canopy_centre);
    vector<DataPoint*> points_to_erase;
    for (DataPoint* p : point_set) {
      float squared_dist = new_canopy_centre->get_squared_dist(*p);
      if (squared_dist < T1 * T1) {
        new_canopy.add_point(p);
      }
      if (squared_dist < T2 * T2) {
        points_to_erase.push_back(p);
      }
    }
    for (DataPoint* p : points_to_erase) {
      point_set.erase(p);
    }
    canopies.push_back(new_canopy);
  }

  // process that picks the new canopy centre
  int root_process = 0;
  int canopy_id = 0;
  int* canopy_id_send_data = new int[scatter_counts[world_rank] / 2];
  memset(canopy_id_send_data, -1, sizeof(int) * scatter_counts[world_rank] / 2);

  while (true) {
    MPI_Barrier(MPI_COMM_WORLD);

    // try to update the root_process
    int prev_root_process = root_process;
    int temp_root_process = root_process;
    if (world_rank == root_process && point_set.empty()) {
      temp_root_process++;
    }
    MPI_Bcast(&temp_root_process, 1, MPI_INT, root_process, MPI_COMM_WORLD);
    root_process = temp_root_process;

    if (root_process >= world_size) {
      // we got through all the processes and can break
      //printf("breaking %d\n", world_rank);
      break;
    }

    if (root_process != prev_root_process) {
      //printf("continue\n");
      continue;
    }

    DataPoint* new_canopy_centre = NULL;
    float* new_canopy_centre_data = new float[num_dimensions];
    if (world_rank == root_process) {
      new_canopy_centre = *(point_set.begin());
      point_set.erase(new_canopy_centre);
      for (int i=0; i<num_dimensions; i++) {
        new_canopy_centre_data[i] = new_canopy_centre->get_val(i);
      }
    }
    //broadcast new canopy centre
    MPI_Bcast(new_canopy_centre_data, num_dimensions, MPI_FLOAT, root_process, MPI_COMM_WORLD);
    if (world_rank != root_process) {
      new_canopy_centre = new DataPoint(new_canopy_centre_data);
    }
    if (world_rank == 0) {
      canopies.push_back(Canopy(new_canopy_centre));
    }
 #ifdef DEBUG
    printf("centre %d ", world_rank);
    new_canopy_centre->print();
 #endif

    vector<DataPoint*> points_to_erase;
    for (DataPoint* p : point_set) {
      float squared_dist = new_canopy_centre->get_squared_dist(*p);
      if (squared_dist < T1 * T1) {
 #ifdef DEBUG
        if (p->get_point_id() >= scatter_counts[world_rank]) {
          printf("ERROR ERROR ERROR %d\n", p->get_point_id());
          continue;
        }
 #endif
        canopy_id_send_data[p->get_point_id()] = canopy_id;
      }
      if (squared_dist < T2 * T2) {
        points_to_erase.push_back(p);
      }
    }
    for (DataPoint* p : points_to_erase) {
      point_set.erase(p);
    }
    canopy_id++;
  }

  // gather the information about canopy points in process 0
  int* canopy_id_data = NULL;
  if (world_rank == 0) {
    canopy_id_data = new int[all_points.size()];
  }

  int* gather_counts = new int[world_size];
  int* gather_displs = new int[world_size];
  sum = 0;
  for (int i=0; i<world_size; i++) {
    gather_counts[i] =  num_points / world_size;
    if (i < num_points % world_size) {
      gather_counts[i]++;
    }
    gather_displs[i] = sum;
    sum += gather_counts[i];
  }

  MPI_Gatherv(canopy_id_send_data, gather_counts[world_rank], MPI_INT,
             canopy_id_data, gather_counts, gather_displs, MPI_INT,
             0, MPI_COMM_WORLD);

  // create new canopy in process 0
  if (world_rank == 0) {
    for (int i=0; i<all_points.size(); i++) {
      if (canopy_id_data[i] != -1) {
        canopies[canopy_id_data[i]].add_point(all_points[i]);
      }
    }
  }
  return canopies;
}

void parse_args(int argc, char** argv) {
  for (int i = 1; i<argc;) {
    /*if (strcmp(argv[i], "-d") == 0) {
      num_dimensions = atoi(argv[i+1]);
      i += 2;
    } else if (strcmp(argv[i], "-n") == 0) {
      num_points = atoi(argv[i+1]);
      i += 2;
    }*/
    if (strcmp(argv[i], "-file") == 0) {
      fileName = string(argv[i+1]);
      i += 2;
    } else if (strcmp(argv[i], "-T1") == 0) {
      T1 = atoi(argv[i+1]);
      i += 2;
    } else if (strcmp(argv[i], "-T2") == 0) {
      T2 = atoi(argv[i+1]);
      i += 2;
    } else if (strcmp(argv[i], "--no-print") == 0) {
      print_canopies = false;
      i++;
    }
  }
}

int DataPoint::point_id_counter = 0;

int main(int argc, char** argv) {
  MPI_Init(NULL, NULL);

  srand(time(NULL));
  parse_args(argc, argv);

  assert(T1 > T2);

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

#ifdef DEBUG
  printf("Rank : %d\n", world_rank);
#endif

  vector<DataPoint*> all_points;
  if (world_rank == 0) {
    //generate_points(all_points, num_points);
    read_points(all_points, fileName, false);
    num_points = all_points.size();
#ifdef DEBUG
    for (const DataPoint* p : all_points) {
      p->print();
    }
#endif
  }

  // call canopy
  vector<Canopy> canopies = canopy_mpi(all_points);

  // print points
  if (world_rank == 0 && print_canopies) {
    /*for (DataPoint* p : all_points) {
      p->print();
    }
    printf("\n");
    */

    //cout << T1 << ' ' << T2 << endl << endl;
    //printf("%f %f\n", T1, T2);

    // print canopies
    cout<<"Resulted clusters: "<<endl;
    for (Canopy& c : canopies) {
      cout<<"Cluster center: ";
      c.print();
      c.printElements();
    }
  }

  MPI_Finalize();
  return 0;
}

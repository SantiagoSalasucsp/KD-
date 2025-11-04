#include <iostream>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <sstream>
#include <string>
#include <functional>
#include <map>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <set>
#include <limits>
#include <numeric>
#include <iomanip>

using namespace std;

class Point
{
public:
    std::vector<double> coords;

    Point() {}
    Point(int k) : coords(k) {}
    Point(const std::vector<double>& coords_) : coords(coords_) {}
    Point(std::vector<double>&& coords_) : coords(std::move(coords_)) {}

    bool operator==(const Point& other) const { return coords == other.coords; }
    bool operator<(const Point& other) const { return coords < other.coords; }
    Point operator+(const Point& other) const {
        std::vector<double> result(coords.size());
        std::transform(coords.begin(), coords.end(), other.coords.begin(), result.begin(), std::plus<>());
        return result;
    }
    Point operator/(double k) const {
        std::vector<double> result(coords.size());
        std::transform(coords.begin(), coords.end(), result.begin(), [k](double c) { return c / k; });
        return result;
    }
    
    double distanceSq(const Point& other) const {
        double dist = 0.0;
        for (size_t i = 0; i < coords.size(); ++i) {
            dist += pow(coords[i] - other.coords[i], 2);
        }
        return dist;
    }
};

namespace std {
    template <>
    struct hash<Point> {
        size_t operator()(const Point& p) const {
            size_t seed = 0;
            for (double val : p.coords) {
                seed ^= (std::hash<double>()(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2));
            }
            return seed;
        }
    };
}

class KDNode
{
public:
    Point point;
    KDNode *left, *right;
    
    KDNode(const Point& p) : point(p), left(NULL), right(NULL) {}
};

class KDTree
{
private:
    KDNode* root;
    int k_dim;

    KDNode* insert(std::vector<Point>& points, int start, int end, int depth)
    {
        if (start >= end) return nullptr;

        int dim = depth % k_dim, median = (start + end) / 2;
        
        std::nth_element(points.begin() + start, points.begin() + median,
                         points.begin() + end,
                         [dim](const Point& a, const Point& b) {
                             return a.coords[dim] < b.coords[dim];
                         });

        KDNode* node = new KDNode(points[median]);
        node->left = insert(points, start, median, depth + 1);
        node->right = insert(points, median + 1, end, depth + 1);

        return node;
    }
    
    Point* nearestPointRecursive(KDNode* root_, const Point& target, int depth, Point* best, double& best_dist_sq)
    {
        if (!root_) return best;

        double dist_sq = root_->point.distanceSq(target);
        if (dist_sq < best_dist_sq) {
            best_dist_sq = dist_sq;
            best = &root_->point;
        }

        int dim = depth % k_dim;
        bool isLeft = target.coords[dim] < root_->point.coords[dim];
        KDNode* next = isLeft ? root_->left : root_->right;
        KDNode* other = isLeft ? root_->right : root_->left;

        best = nearestPointRecursive(next, target, depth + 1, best, best_dist_sq);
        
        double diff = target.coords[dim] - root_->point.coords[dim];
        if (diff * diff < best_dist_sq) {
            best = nearestPointRecursive(other, target, depth + 1, best, best_dist_sq);
        }
        return best;
    }

public:
    KDTree(const std::vector<Point>& points, int k_dim_val) : k_dim(k_dim_val)
    {
        std::vector<Point> temp = points;
        root = insert(temp, 0, temp.size(), 0);
    }
    
    Point* nearestPoint(const Point& point)
    {
        double best_dist_sq = std::numeric_limits<double>::max();
        return nearestPointRecursive(root, point, 0, nullptr, best_dist_sq);
    }
};

class Kmeans {
private:
    vector<pair<double, double>> bounds;
    vector<Point> points;

    set<Point> chooseCentroids(int k) {
        set<Point> centroids;
        while (centroids.size() < k) {
            int randNum = rand() % points.size();
            centroids.emplace(points[randNum]);
        }
        return centroids;
    }

    Point calcMean(const vector<Point>& clusterPoints) {
        if (clusterPoints.empty()) return Point(points[0].coords.size());
        
        Point sum = clusterPoints[0];
        for (size_t i = 1; i < clusterPoints.size(); ++i) {
            sum = sum + clusterPoints[i];
        }
        return sum / clusterPoints.size();
    }
    
    unordered_map<Point, vector<Point>> clustering(const set<Point>& centroids) {
        unordered_map<Point, vector<Point>> clusters;
        
        for(const auto& c : centroids) clusters[c] = {};
        
        for (const auto& p : points) {
            double min_dist = std::numeric_limits<double>::max();
            Point closest_centroid;
            
            for (const auto& c : centroids) {
                double dist = p.distanceSq(c);
                if (dist < min_dist) {
                    min_dist = dist;
                    closest_centroid = c;
                }
            }
            if (clusters.count(closest_centroid)) {
                clusters[closest_centroid].push_back(p);
            }
        }
        return clusters;
    }
    
    unordered_map<Point, vector<Point>> clustering_kd(const set<Point>& centroids) {
        unordered_map<Point, vector<Point>> clusters;
        
        for(const auto& c : centroids) clusters[c] = {};

        if (centroids.empty()) return clusters;
        int k_dim = centroids.begin()->coords.size();
        KDTree tree(std::vector<Point>(centroids.begin(), centroids.end()), k_dim);

        for (const auto& p : points) {
            Point* closest_centroid_ptr = tree.nearestPoint(p);
            if (closest_centroid_ptr) {
                clusters[*closest_centroid_ptr].push_back(p);
            }
        }
        return clusters;
    }
    
    set<Point> recalculateCentroids(const unordered_map<Point, vector<Point>>& clusters) {
        set<Point> newCentroids;
        for (const auto& pair : clusters) {
            newCentroids.insert(calcMean(pair.second));
        }
        return newCentroids;
    }

public:
    Kmeans(const std::string& filename, int lim) {
        ifstream file(filename);
        if (!file) throw runtime_error("No se pudo abrir el archivo csv: " + filename);

        string line;
        getline(file,line);

        for (int lineCount = 0; lineCount < lim && getline(file, line); ++lineCount) {
            istringstream iss(line); string svalue;
            vector<double> coords;
            
            while (getline(iss, svalue, ',')) {
                coords.push_back(std::stod(svalue));
            }
            if (!coords.empty()) {
                points.push_back(std::move(coords));
            }
        }
        if (points.empty()) throw runtime_error("No se cargaron puntos.");
    }

    unordered_map<Point, vector<Point>> kmeans(int k, int iterations) {
        unordered_map<Point, vector<Point>> clusters;
        set<Point> centroids = chooseCentroids(k);

        for (size_t i = 0; i < iterations; i++) {
            clusters = clustering(centroids);
            set<Point> newCentroids = recalculateCentroids(clusters);
            if (centroids == newCentroids) break;
            centroids = std::move(newCentroids);
        }
        return clusters;
    }

    unordered_map<Point, vector<Point>> kmeans_kd(int k, int iterations) {
        unordered_map<Point, vector<Point>> clusters;
        set<Point> centroids = chooseCentroids(k);

        for (size_t i = 0; i < iterations; i++) {
            clusters = clustering_kd(centroids);
            set<Point> newCentroids = recalculateCentroids(clusters);
            if (centroids == newCentroids) break;
            centroids = std::move(newCentroids);
        }
        return clusters;
    }
    
    void ExportClustersToJson(const unordered_map<Point, vector<Point>>& clusters, const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "ERROR: No se pudo crear el archivo JSON de clusters en: " << filename << std::endl;
            return;
        }
        
        file << "{\n  \"clusters\": {\n";
        int cluster_idx = 0;
        
        for (const auto& entry : clusters) {
            file << "    \"" << cluster_idx << "\": [\n      [";
            for (size_t i = 0; i < entry.first.coords.size(); ++i) {
                file << std::fixed << std::setprecision(4) << entry.first.coords[i] << (i < entry.first.coords.size() - 1 ? ", " : "");
            }
            
            file << "],\n      [\n";
            for (size_t j = 0; j < entry.second.size(); ++j) {
                file << "        [";
                for (size_t i = 0; i < entry.second[j].coords.size(); ++i) {
                    file << std::fixed << std::setprecision(4) << entry.second[j].coords[i] << (i < entry.second[j].coords.size() - 1 ? ", " : "");
                }
                file << "]" << (j < entry.second.size() - 1 ? ",\n" : "\n");
            }
            file << "      ]\n    ]" << (cluster_idx < clusters.size() - 1 ? ",\n" : "\n");
            cluster_idx++;
        }
        file << "  }\n}\n";
    }
};

double executionTime(const std::function<void()>& func) {
    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(end - start).count();
}


struct AnalysisTest
{
    int constant;
    map<int, vector<pair<double, double>>> times_Kmeans;
};

void exportAnalysisToJsonManual(const std::vector<AnalysisTest>& analysis, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "ERROR: No se pudo crear el archivo JSON en: " << filename << std::endl;
        return;
    }
    
    file << "[\n";
    for (size_t i = 0; i < analysis.size(); ++i) {
        const auto& test = analysis[i];
        file << "  {\n";
        file << "    \"constant\": " << test.constant << ",\n";
        file << "    \"times_Kmeans\": {\n";

        size_t x_count = 0;
        for (const auto& [x_val, times] : test.times_Kmeans) {
            file << "      \"" << x_val << "\": [\n";
            for (size_t j = 0; j < times.size(); ++j) {
                file << "        { \"Brute Force\": " << std::fixed << std::setprecision(4) << times[j].first
                     << ", \"KD-Tree\": " << std::fixed << std::setprecision(4) << times[j].second << " }"
                     << (j < times.size() - 1 ? ",\n" : "\n");
            }
            file << "      ]" << (x_count < test.times_Kmeans.size() - 1 ? ",\n" : "\n");
            x_count++;
        }
        
        file << "    }\n";
        file << "  }" << (i < analysis.size() - 1 ? ",\n" : "\n");
    }
    file << "]\n";

    file.close();
    std::cout << "Datos de analisis guardados en: " << filename << std::endl;
}

std::vector<AnalysisTest> runAnalysisTests(const std::string& data_path,
                                         const std::vector<int>& constants,
                                         const std::vector<int>& variables,
                                         int iteraciones,
                                         bool is_N_variable,
                                         int repeticiones)
{
    std::vector<AnalysisTest> analysis;
    
    for (auto& constV : constants) {
        AnalysisTest times{ constV };

        for (auto& var : variables) {
            
            for (int r = 0; r < repeticiones; ++r) {
                
                int N = is_N_variable ? var : 2400;
                int K = is_N_variable ? constV : var;
                
                srand(r + constV * 100);
                
                Kmeans kc(data_path, N);
                double timeBF = executionTime([&]() { kc.kmeans(K, iteraciones); });
                
                Kmeans kc_kd(data_path, N);
                srand(r + constV * 100);
                double timeKD = executionTime([&]() { kc_kd.kmeans_kd(K, iteraciones); });
                
                times.times_Kmeans[var].emplace_back(timeBF, timeKD);
            }
            std::cout << (is_N_variable ? "N=" : "K=") << var << " completado. " << std::endl;
        }
        analysis.push_back(times);
    }
    return analysis;
}

int main() {
    const string DATA_FILE_PATH = "/Users/santiagosalas/Desktop/Punto/data2k.csv";
    const string OUTPUT_PATH = "/Users/santiagosalas/Desktop/resultados/";
    const int MAX_ITERATIONS = 100;
    const int NUM_REPETITIONS = 10;
    
    try {
        vector<int> K_CONSTANTES = { 5, 15 };
        vector<int> N_VARIABLES = { 500, 800, 1100, 1400, 1700, 2000, 2400 };

        cout << "--- Ejecutando Analisis: N Variable (Tiempo vs N) ---" << endl;
        auto analysis_N = runAnalysisTests(DATA_FILE_PATH, K_CONSTANTES, N_VARIABLES, MAX_ITERATIONS, true, NUM_REPETITIONS);
        exportAnalysisToJsonManual(analysis_N, OUTPUT_PATH + "tiempos_N_variable.json");

        vector<int> N_CONSTANTES = { 1100, 1500 };
        vector<int> K_VARIABLES = { 25, 50, 75, 100, 125, 150, 175, 200 };

        cout << "\n--- Ejecutando Analisis: K Variable (Tiempo vs K) ---" << endl;
        auto analysis_K = runAnalysisTests(DATA_FILE_PATH, N_CONSTANTES, K_VARIABLES, MAX_ITERATIONS, false, NUM_REPETITIONS);
        exportAnalysisToJsonManual(analysis_K, OUTPUT_PATH + "tiempos_K_variable.json");

        cout << "\n--- Ejecutando para Visualizacion (K=18, N=2400) ---" << endl;
        srand(1);
        Kmeans km_viz(DATA_FILE_PATH, 2400);
        auto clusters_final = km_viz.kmeans_kd(18, MAX_ITERATIONS);
        km_viz.ExportClustersToJson(clusters_final, OUTPUT_PATH + "clusters_final_kd.json");
        cout << "Clusters finales guardados en: " << OUTPUT_PATH << "clusters_final_kd.json" << endl;
    } catch (const std::exception& e) {
        cerr << "Excepcion fatal en la ejecucion: " << e.what() << endl;
        cerr << "Asegurate de que la ruta de salida exista: " << OUTPUT_PATH << endl;
    } catch (...) {
        cerr << "Excepcion desconocida en la ejecucion." << endl;
    }
    
    return 0;
}

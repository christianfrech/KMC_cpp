diff --git a/hpp_files/math_func.hpp b/hpp_files/math_func.hpp
index 85b69d6..eebceb4 100644
--- a/hpp_files/math_func.hpp
+++ b/hpp_files/math_func.hpp
@@ -8,46 +8,67 @@
 #include <cmath>
 #include <math.h>
 
+#include <iostream>
+#include <sstream>
+#include <stdio.h>
+#include <string>
+#include <vector>
+#include <cstdlib>
+#include <fstream>
 
+#include "classes_strucs.hpp"
+
+/*! \brief Calculates the power of an integer base raised to an exponent.
+ * \param [in] base The base integer.
+ * \param [in] exponent The exponent integer.
+ * \return The result of base raised to the power of exponent.
+ */
 int exp_int(int base, int exponent) {
     int res = 0;
-    
-    if (exponent == 0) {return 1;}
-    else {
-        for (int i=0; i<(exponent); i++) {
+
+    if (exponent == 0) {
+        return 1;
+    } else {
+        for (int i = 0; i < exponent; i++) {
             if (i == 0) {
                 res = base;
+            } else {
+                res = res * base;
             }
-            else {res = res*base;}
         }
     }
     return res;
 }
 
+/*! \brief Sums the elements in a vector of doubles.
+ * \param [in] vec The vector of double values to be summed.
+ * \return The sum of all elements in the vector.
+ */
 double sum_float(std::vector<double> vec) {
     double sum = 0;
-    for (int i=0; i<(int)vec.size(); i++) {
-        sum += vec[i]; 
-        
-        }
-
+    for (int i = 0; i < (int)vec.size(); i++) {
+        sum += vec[i];
+    }
     return sum;
 }
 
+/*! \brief Calculates the binomial coefficient, "n choose r."
+ * \param [in] n The total number of items.
+ * \param [in] r The number of items to choose.
+ * \return The binomial coefficient (n choose r).
+ */
 int NCR(int n, int r) {
     if (r == 0) return 1;
 
-    /*
-     Extra computation saving for large R,
-     using property:
-     N choose R = N choose (N-R)
-    */
-    if (r > n / 2) return NCR(n, n - r); 
+    /*! 
+     * \note Extra computation saving for large r, using property:
+     * N choose R = N choose (N - R)
+     */
+    if (r > n / 2) return NCR(n, n - r);
 
-    long res = 1; 
+    long res = 1;
 
-    for (int k = 1; k <= r; ++k)
-    {
+    for (int k = 1; k <= r; ++k) {
         res *= n - k + 1;
         res /= k;
     }
@@ -55,8 +76,14 @@ int NCR(int n, int r) {
     return res;
 }
 
+/*! \brief Finds the greatest common factor (GCF) of two integers.
+ * \param [in] a The first integer.
+ * \param [in] b The second integer.
+ * \return The greatest common factor of a and b.
+ */
 int find_gcf(int a, int b) {
     std::cout << "a: " << a << " b: " << b << "\n";
+
     // Everything divides 0 
     if (a == 0) 
         return b; 
@@ -73,14 +100,25 @@ int find_gcf(int a, int b) {
     return find_gcf(a, b - a);
 }
 
+/*! \brief Creates a partition by finding the GCF of two integers.
+ * \param [in] a The first integer.
+ * \param [in] b The second integer.
+ * \return The GCF of a and b, with conditional adjustments.
+ */
 int create_partition(int a, int b) {
-    if ((a%2 == 0) && (b%2 == 0)) {
-        if (a == b) { return find_gcf(a, (int)(b/2));}
+    if ((a % 2 == 0) && (b % 2 == 0)) {
+        if (a == b) { return find_gcf(a, (int)(b / 2)); }
         else return find_gcf(a, b);
-    }
-    else return find_gcf(a, b);
+    } else return find_gcf(a, b);
 }
 
+/*! \brief reduces two integers modulo given bounds, ensuring they fall within specified bounds.
+ * \param [in] i The first integer.
+ * \param [in] j The second integer.
+ * \param [in] bound_x The bound for the first integer.
+ * \param [in] bound_y The bound for the second integer.
+ * \return A vector of size_t values containing the modulo results for i and j.
+ */
 std::vector<size_t> mod_with_bounds(int i, int j, int bound_x, int bound_y) {
     std::vector<size_t> idxs(2);
 
@@ -90,6 +128,11 @@ std::vector<size_t> mod_with_bounds(int i, int j, int bound_x, int bound_y) {
     return idxs;
 }
 
+/*! \brief Reduces an integer modulo a given bound, ensuring it falls within the bound.
+ * \param [in] i The integer to modify.
+ * \param [in] bound_x The bound.
+ * \return A size_t value containing the modulo result of i.
+ */
 size_t mod_with_bounds(int i, int bound_x) {
     size_t idx = (size_t)((i % bound_x + bound_x) % bound_x);
     return idx;
diff --git a/hpp_files/str_func.hpp b/hpp_files/str_func.hpp
index 157dee2..805ee3c 100644
--- a/hpp_files/str_func.hpp
+++ b/hpp_files/str_func.hpp
@@ -8,158 +8,174 @@
 
 #include "classes_strucs.hpp"
 
-
+/*! \brief Structure for portability of information of a site read in from custom region input file . */
 struct reg_line_struct {
-    public:
-        /*! \brief Constructor
-        */
-
-        reg_line_struct(std::string lattice_type, int x, int y, int z):
-            lattice_type_(lattice_type),
-            x_(x),
-            y_(y),
-            z_(z)
-            {}
-
-        /*! \return get lattice_type field of struc */
-        std::string get_latice_type() const {
-            return lattice_type_;
-        }
-
-        /*! \return get x field of struc */
-        int get_x() const {
-            return x_;
-        }
-
-        /*! \return get y field of struc */
-        int get_y() const {
-            return y_;
-        }
-
-        /*! \return get y field of struc */
-        int get_z() const {
-            return z_;
-        }
-        
-    private:
-        std::string lattice_type_; int x_; int y_; int z_;
+public:
+    /*! \brief Constructor
+     * \param [in] lattice_type Type of lattice.
+     * \param [in] x X-coordinate.
+     * \param [in] y Y-coordinate.
+     * \param [in] z Z-coordinate.
+     */
+    reg_line_struct(std::string lattice_type, int x, int y, int z)
+        : lattice_type_(lattice_type), x_(x), y_(y), z_(z) {}
+
+    /*! \brief Get the lattice type.
+     * \return String representing the lattice type.
+     */
+    std::string get_latice_type() const { return lattice_type_; }
+
+    /*! \brief Get the X coordinate.
+     * \return Integer X coordinate.
+     */
+    int get_x() const { return x_; }
+
+    /*! \brief Get the Y coordinate.
+     * \return Integer Y coordinate.
+     */
+    int get_y() const { return y_; }
+
+    /*! \brief Get the Z coordinate.
+     * \return Integer Z coordinate.
+     */
+    int get_z() const { return z_; }
+
+private:
+    std::string lattice_type_; ///< Lattice type.
+    int x_; ///< X coordinate.
+    int y_; ///< Y coordinate.
+    int z_; ///< Z coordinate.
 };
- 
-typedef struct reg_line_struct reg_line_struct;
 
+typedef struct reg_line_struct reg_line_struct;
 
+/*! \brief Structure for storing information of location of atom from geometry input file. */
 struct line_struct {
-    public:
-        /*! \brief Constructor
-        */
-
-        line_struct(std::string lattice_pos, double x, double y, double z, int atype):
-            lattice_pos_(lattice_pos),
-            x_(x),
-            y_(y),
-            z_(z),
-            atype_(atype)
-            {}
-
-        /*! \return get lattice_type field of struc */
-        std::string get_latice_pos() const {
-            return lattice_pos_;
-        }
-
-        /*! \return get x field of struc */
-        double get_x() const {
-            return x_;
-        }
-
-        /*! \return get y field of struc */
-        double get_y() const {
-            return y_;
-        }
-
-        /*! \return get z field of struc */
-        double get_z() const {
-            return z_;
-        }
-
-        /*! \return get z field of struc */
-        int get_atype() const {
-            return atype_;
-        }
-        
-    private:
-        std::string lattice_pos_; double x_; double y_; double z_; int atype_;
+public:
+    /*! \brief Constructor
+     * \param [in] lattice_pos Lattice position identifier.
+     * \param [in] x X-coordinate.
+     * \param [in] y Y-coordinate.
+     * \param [in] z Z-coordinate.
+     * \param [in] atype Atom type.
+     */
+    line_struct(std::string lattice_pos, double x, double y, double z, int atype)
+        : lattice_pos_(lattice_pos), x_(x), y_(y), z_(z), atype_(atype) {}
+
+    /*! \brief Get the lattice position identifier.
+     * \return String representing the lattice position.
+     */
+    std::string get_latice_pos() const { return lattice_pos_; }
+
+    /*! \brief Get the X coordinate.
+     * \return Double X coordinate.
+     */
+    double get_x() const { return x_; }
+
+    /*! \brief Get the Y coordinate.
+     * \return Double Y coordinate.
+     */
+    double get_y() const { return y_; }
+
+    /*! \brief Get the Z coordinate.
+     * \return Double Z coordinate.
+     */
+    double get_z() const { return z_; }
+
+    /*! \brief Get the atom type.
+     * \return Integer representing the atom type.
+     */
+    int get_atype() const { return atype_; }
+
+private:
+    std::string lattice_pos_; ///< Lattice position.
+    double x_; ///< X coordinate.
+    double y_; ///< Y coordinate.
+    double z_; ///< Z coordinate.
+    int atype_; ///< Atom type.
 };
- 
-typedef struct line_struct line_struct;
 
+typedef struct line_struct line_struct;
 
+/*! \brief Structure for cataloging rate data for regions. */
 struct ratecatalog_struct {
-    public:
-        /*! \brief Constructor
-        */
-
-        ratecatalog_struct(std::vector< std::vector<int> > configs, 
-        std::vector< std::vector< std::vector<double> > > energies, 
-        std::vector< std::vector< std::vector< std::vector<double> > > > regions_cat, 
-        int region_num):
-            configs_(configs),
-            energies_(energies),
-            regions_cat_(regions_cat),
-            region_num_(region_num)
-            {}
-
-        /*! \return get x field of struc */
-        std::vector< std::vector<int> > get_configs() const {
-            return configs_;
-        }
-
-        /*! \return get y field of struc */
-        std::vector< std::vector< std::vector<double> > > get_energies() const {
-            return energies_;
-        }
-
-        /*! \return get z field of struc */
-        std::vector< std::vector< std::vector< std::vector<double> > > > get_regions_cat() const {
-            return regions_cat_;
-        }
-
-        /*! \return get z field of struc */
-        int get_region_num() const {
-            return region_num_;
-        }
-        
-    private:
-        std::vector< std::vector<int> > configs_; 
-        std::vector< std::vector< std::vector<double> > > energies_; 
-        std::vector< std::vector< std::vector< std::vector<double> > > > regions_cat_; 
-        int region_num_;
+public:
+    /*! \brief Constructor
+     * \param [in] configs Vector of binary strings representing configurations.
+     * \param [in] energies Vector of energy data.
+     * \param [in] regions_cat Vector of region catalog data.
+     * \param [in] region_num Number of regions.
+     */
+    ratecatalog_struct(std::vector<std::vector<int>> configs,
+                       std::vector<std::vector<std::vector<double>>> energies,
+                       std::vector<std::vector<std::vector<std::vector<double>>>> regions_cat,
+                       int region_num)
+        : configs_(configs), energies_(energies), regions_cat_(regions_cat), region_num_(region_num) {}
+
+    /*! \brief Get the vector of binary strings representing configurations.
+     * \return Vector of configuration data.
+     */
+    std::vector<std::vector<int>> get_configs() const { return configs_; }
+
+    /*! \brief Get the vector of energies corresponding to each NN-configuration.
+     * \return Vector of energy data.
+     */
+    std::vector<std::vector<std::vector<double>>> get_energies() const { return energies_; }
+
+    /*! \brief Get the rate catalog corresponding to regions.
+     * \return Vector of rate catalogs (tripley-nested vectors) corresponding to regions.
+     */
+    std::vector<std::vector<std::vector<std::vector<double>>>> get_regions_cat() const { return regions_cat_; }
+
+    /*! \brief Get the number of regions.
+     * \return Integer number of regions.
+     */
+    int get_region_num() const { return region_num_; }
+
+private:
+    std::vector<std::vector<int>> configs_; ///< Vector of configurations.
+    std::vector<std::vector<std::vector<double>>> energies_; ///< Vector of energy data.
+    std::vector<std::vector<std::vector<std::vector<double>>>> regions_cat_; ///< Region catalog data.
+    int region_num_; ///< Number of regions.
 };
- 
+
 typedef struct ratecatalog_struct ratecatalog_struct;
 
+/*! \brief Tokenizes a string based on a given delimiter.
+ * \param [in] s The string to tokenize.
+ * \param [in] delimiter The delimiter string.
+ * \return Vector of tokens.
+ */
 std::vector<std::string> tokenizer(std::string s, std::string delimiter) {
     size_t pos_start = 0, pos_end, delim_len = delimiter.length();
     std::string token;
     std::vector<std::string> res;
 
     while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
-        token = s.substr (pos_start, pos_end - pos_start);
+        token = s.substr(pos_start, pos_end - pos_start);
         pos_start = pos_end + delim_len;
-        res.push_back (token);
+        res.push_back(token);
     }
-    if (s.substr (pos_start, pos_end - pos_start) != "") {
+    if (s.substr(pos_start, pos_end - pos_start) != "") {
         res.push_back(s.substr(pos_start));
     }
-    else {}
-    
+
     return res;
 }
 
+/*! \brief Splits a string into a vector of characters.
+ * \param [in] s The string to split.
+ * \return Vector of characters.
+ */
 std::vector<char> split_by_char(std::string s) {
     std::vector<char> v(s.begin(), s.end());
     return v;
 }
 
+/*! \brief Checks if a string contains only alphabetic characters.
+ * \param [in] str The string to check.
+ * \return True if the string is alphabetic, false otherwise.
+ */
 bool str_isalpha(std::string str) {
     auto it = std::find_if(str.begin(), str.end(), [](char const &c) {
         return !std::isalpha(c);
@@ -167,8 +183,11 @@ bool str_isalpha(std::string str) {
     return it == str.end();
 }
 
+/*! \brief Parses a line to create a line_struct object.
+ * \param [in] line The input line to parse.
+ * \return A line_struct object populated with parsed data.
+ */
 line_struct parse_line(std::string line) {
-
     std::vector<std::string> toks = tokenizer(line, " ");
     std::string lattice_pos = (toks[0]);
 
@@ -177,57 +196,63 @@ line_struct parse_line(std::string line) {
     double z = std::stof(toks[3]);
     int atype = std::stoi(toks[4]);
 
-    line_struct output_vals(lattice_pos,x,y,z,atype);
+    line_struct output_vals(lattice_pos, x, y, z, atype);
     return output_vals;
 }
 
+/*! \brief Parses a line from custom region input file to create a reg_line_struct object.
+ * \param [in] line The input line to parse.
+ * \return A reg_line_struct object populated with parsed data.
+ */
 reg_line_struct parse_reg_line(std::string line) {
-
     std::vector<std::string> toks = tokenizer(line, " ");
     std::string lattice_pos = (toks[0]);
     double x = std::stof(toks[1]);
     double y = std::stof(toks[2]);
     double z = std::stof(toks[3]);
 
-    reg_line_struct output_vals(lattice_pos,x,y,z);
+    reg_line_struct output_vals(lattice_pos, x, y, z);
     return output_vals;
 }
 
-
+/*! \brief Writes matrix data to a file.
+ * \param [in] filename The name of the output file.
+ * \param [in] values Matrix data to write to the file.
+ */
 void write_to_file(std::string filename, Matrix<int> values) {
     std::ofstream out_file;
     out_file.open(filename);
-    std::string s;
 
     std::cout << "writing filename: " << filename << "\n";
 
     if (out_file.is_open()) {
-        for (int i=0; i<(int)values.rows(); i++) {
-            for (int j=0; j<(int)values.cols(); j++) {     
+        for (int i = 0; i < (int)values.rows(); i++) {
+            for (int j = 0; j < (int)values.cols(); j++) {
                 out_file << values[i][j] << " ";
             }
-            out_file << "\n";  
+            out_file << "\n";
         }
     }
 
     out_file.close();
 }
 
-
-void write_to_file(std::string filename, std::vector< std::vector<int> > values) {
+/*! \brief Writes vector data to a file.
+ * \param [in] filename The name of the output file.
+ * \param [in] values Vector of integer vectors to write to the file.
+ */
+void write_to_file(std::string filename, std::vector<std::vector<int>> values) {
     std::ofstream out_file;
     out_file.open(filename);
-    std::string s;
 
     if (out_file.is_open()) {
-        for (int i=0; i<(int)values.size(); i++) {
-            for (int j=0; j<(int)values[i].size(); j++) { 
+        for (int i = 0; i < (int)values.size(); i++) {
+            for (int j = 0; j < (int)values[i].size(); j++) {
                 out_file << values[i][j] << " ";
             }
-        out_file << "\n";
+            out_file << "\n";
         }
     }
 
     out_file.close();
 }
-
diff --git a/hpp_files/vec_func.hpp b/hpp_files/vec_func.hpp
index 76dbc9b..d146282 100644
--- a/hpp_files/vec_func.hpp
+++ b/hpp_files/vec_func.hpp
@@ -7,27 +7,8 @@
 #include <fstream>
 #include "print_func.hpp"
 
-/*
-template <typename A, typename B>
-A vect_create_2D(size_t N, size_t M) {
-    A vec_out(N);
 
-    for(int i = 0; i < (int)N; i++) { 
-        vec_out[i] = B(M);
-    }
-
-    return vec_out;
-}
-
-template <typename A, typename B, typename C>
-A vect_create_3D_float(int L, int N, int M) {
-
-    std::vector< std::vector< std::vector<double> > > vec_out(L, std::vector< std::vector<double> >(N, std::vector<double>(M, 0)));
-    return vec_out;
-}
-*/
-
-std::vector< std::vector<int> > vect_create_2D(size_t N, size_t M) {
+std::vector< std::vector<int> > vect_create_2D(size_t N, size_t M, double val) {
     std::vector< std::vector<int> > vec_out(N);
 
     for(int i = 0; i < (int)N; i++) { 
@@ -84,13 +65,12 @@ std::vector<int> vec_where_1D(std::vector<int> vec1, int value) {
     int len1 = vec1.size();
 
     for (int i=0;i<len1;i++) { 
-            curr = vec1[i];
-            if (curr == value) {
-                coordinates.push_back(i);
-            } 
-        }
-    //std::cout << "coordinates: \n";
-    //print_1Dvector(coordinates);
+        curr = vec1[i];
+        if (curr == value) {
+            coordinates.push_back(i);
+        } 
+    }
+        
     return coordinates;
 }
 
@@ -251,7 +231,6 @@ std::vector< std::vector<int> > FourD_idxs(std::vector<int> start_vec, std::vect
         for (int j=start_vec[1]; j<(int)(end_vec[1]); j++) {
             for (int k=start_vec[2]; k<((int)end_vec[2]); k++) {
                 for (int l=start_vec[3]; l<(int)(end_vec[3]); l++) { 
-                    //std::cout << "[ " << i << " " << j << " " << k << " " << l << " ] \n";
                     idx = {i, j, k, l};
                     idxs.push_back(idx);
                 }
@@ -320,7 +299,7 @@ std::vector<double> slice_1Dvec_float(std::vector<double> vec, int i1_start, int
     }
     return vec;
 }
-
+/*
 std::vector<std::string> slice_1Dvec_str(std::vector<std::string> vec, int i1_start, int i1_end) {
     int size1 = i1_end - i1_start;
     std::vector<std::string> slice(size1);
@@ -331,6 +310,7 @@ std::vector<std::string> slice_1Dvec_str(std::vector<std::string> vec, int i1_st
     }
     return slice;
 }
+*/
 
 int find_max_element(std::vector<double> * max_rates) {
     int max_idx;
@@ -410,12 +390,20 @@ Matrix<int> comparison(Matrix<int> &mat1, Matrix<int> &mat2) {
         } 
     }
     
-    //std::cout << "reshaping: \n";
     mat_out.reshape((idx), 4);
-    //std::cout << "nonzero elem: " << elem << "\n";
     return mat_out;
 }
 
+template <typename T, , typename U>
+bool is_in(T v, U sub_v) {
+    for (int i=0; i<(int)v.size(); i++) {
+        if (v[i] == sub_v) {return true;}
+    }
+    
+    return false;
+}
+
+
 bool is_in(std::vector<std::vector<size_t>> v, std::vector<size_t> sub_v) {
     for (int i=0; i<(int)v.size(); i++) {
         if (v[i] == sub_v) {return true;}
@@ -447,4 +435,21 @@ bool is_in(std::vector< int > check_vec, int elem) {
         if (check_vec[i] == elem) return true;
     }
     return false;
-}
\ No newline at end of file
+}
+
+/*
+template<typename Container>
+void is_in(const Container &c)
+{
+    for (const auto &i:c)
+    {
+        read_vectors(i);
+    }
+}
+
+template<>
+void is_in(const vector<int> &container){
+    for(auto i:container)
+        cout<<i<<endl;
+}
+*/
\ No newline at end of file
diff --git a/kmc_cpp/kmc_stripping_parallel.hpp b/kmc_cpp/kmc_stripping_parallel.hpp
index 62f5b4b..cff6c56 100644
--- a/kmc_cpp/kmc_stripping_parallel.hpp
+++ b/kmc_cpp/kmc_stripping_parallel.hpp
@@ -29,9 +29,20 @@
  lattice */
 class Lattice {
 
+    std::vector< std::vector<int> > diag_directions;
+    std::vector< std::vector<int> > edge_directions;
+
     std::map<int, double> rate_typedict;
 
+    std::vector<bool> periodic;
+
+    std::vector<int> ver_tuples_cumsum;
+    std::vector<int> bc_tuples_cumsum;
+    std::vector<int> ver_edge_tuples_cumsum; 
+    std::vector<int> bc_edge_tuples_cumsum;
+
     public:
+
         std::vector<Region*> regions;
         std::vector<int> sublattice_dim;
         std::vector<int> total_dims;
@@ -53,9 +64,6 @@ class Lattice {
         std::vector<double> rates;
         std::vector<double> rate_cumsum;
         std::vector< std::vector< std::vector< std::vector<double> > > > region_energies;
-
-        Matrix<int> diag_directions;
-        Matrix<int> edge_directions;
         Matrix<int> moves_coords; 
         Matrix<int> moves_shifts;
         Matrix<int> moves_lattice;
@@ -93,10 +101,9 @@ class Lattice {
         int conflict_done_flag;
         int void_threshold;
         double void_barrier;
-        int solo_vacs;
             
 
-        Lattice(int xdim, int ydim, int zdim, int num_vacancies, int num_regions, int number_procs, int x_neigh, int y_neigh, int xtot, int ytot, int ztot, int rank_in):
+        Lattice(int xdim, int ydim, int zdim, int num_vacancies, int num_regions, int number_procs, int x_neigh, int y_neigh, int xtot, int ytot, int ztot):
             vertex_sites((size_t)1, (size_t)xdim, (size_t)ydim, (size_t)zdim),
             vacancies((size_t)2, (size_t)xdim, (size_t)ydim, (size_t)zdim),
             bc_sites((size_t)1, (size_t)xdim, (size_t)ydim, (size_t)zdim),
@@ -105,11 +112,7 @@ class Lattice {
             proc_neg_x_neighbors((size_t)1, (size_t)2, (size_t)y_neigh, (size_t)zdim), 
             proc_neg_y_neighbors((size_t)1, (size_t)2, (size_t)x_neigh, (size_t)zdim), 
             proc_pos_x_neighbors((size_t)1, (size_t)2, (size_t)y_neigh, (size_t)zdim), 
-            proc_pos_y_neighbors((size_t)1, (size_t)2, (size_t)x_neigh, (size_t)zdim), 
-
-            diag_directions((size_t)8, (size_t)3),
-            edge_directions((size_t)8, (size_t)3),
-
+            proc_pos_y_neighbors((size_t)1, (size_t)2, (size_t)x_neigh, (size_t)zdim),  
             proc_neighbors((size_t)number_procs, (size_t)8),
             moves_coords((size_t)(14 * num_vacancies), 4),
             moves_shifts((size_t)(14 * num_vacancies), 3),
@@ -124,12 +127,14 @@ class Lattice {
             configs_111(1, (size_t)exp_int(2,8)),
             configs_100(1, (size_t)exp_int(2,14)),
             vacancies_pos((size_t)num_vacancies, 4),
-            mt_obj((unsigned int)(std::chrono::high_resolution_clock::now().time_since_epoch().count() + rank_in)),
+            mt_obj((unsigned int)(std::chrono::high_resolution_clock::now().time_since_epoch().count())),
             //mt_obj((unsigned int)725863834569),
             x_rand(0, (size_t)xdim),
             y_rand(0, (size_t)ydim)
 
             {
+                diag_directions = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}, {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}};
+                edge_directions = {{0,0,1}, {-1,0,0}, {0,-1,0}, {0,1,0}, {1,0,0}, {0,0,-1}};
                 sublattice_dim = {xdim,ydim,zdim};
                 total_dims = {xtot,ytot,ztot};
                 t = 0;
@@ -138,11 +143,25 @@ class Lattice {
                 ghost_done_tag = 10;
                 par_done_tag = 11;
                 conflict_done_flag = 12;
-                rank = rank_in;
             }
 
-        /*
-        subroutine for adding move information to matrices stored as attributes of Lattice struc
+        /**
+        * @brief Subroutine for adding move information to matrices stored as attributes of Lattice structure.
+        * 
+        * This function calculates the new coordinates, shifts, and lattice for a move, then stores
+        * this information in the relevant arrays. It also computes the rate for the move using a rate 
+        * constant function and updates the cumulative rate.
+        * 
+        * @param i            Initial i coordinate.
+        * @param j            Initial j coordinate.
+        * @param k            Initial k coordinate.
+        * @param l            Initial l coordinate.
+        * @param curr_move_num Current move number to be processed.
+        * @param direc_sign   Direction sign for the move.
+        * @param s            Index of direction in diag_directions/edge_directions.
+        * @param idx          Index of vacancy in vacancies_pos Matrix.
+        * @param lattice      Integer corresponding to move type .
+        * @return int         Updated move number after processing.
         */
         int add_move(int i, int j, int k, int l, int curr_move_num, int direc_sign, int s, int idx, int lattice, int NN_count) {
             double rate = 0; 
@@ -188,8 +207,20 @@ class Lattice {
             return curr_move_num;
         }
 
-        /*
-        routine to check if a adjacent site is unoccupied for move
+        /**
+        * @brief Routine to check if an adjacent site is unoccupied for a move.
+        * 
+        * This function checks whether the target position in the lattice is free by verifying
+        * if the corresponding site in the vacancy matrix is occupied.
+        * 
+        * @param i           Initial i coordinate.
+        * @param j           Initial j coordinate.
+        * @param k           Initial k coordinate.
+        * @param l           Initial l coordinate.
+        * @param direc_sign  Direction sign for the move.
+        * @param s           Index of direction in diag_directions/edge_directions.
+        * @param lattice     Integer corresponding to move type (0, 1, 2, or 3).
+        * @return bool       True if the site is has no vacancy, false if occupied with vacancy.
         */
         bool check_move_free(int i, int j, int k, int l, int direc_sign, int s, int lattice) {
             int i1; int i2; int i3; int i4;
@@ -211,9 +242,18 @@ class Lattice {
 
             return true;
         }
+
         
-        /*
-        routine to count the number of vacancies in NN shell
+        /** 
+        * @brief Routine to count the number of vacancies in the nearest neighbor (NN) shell.
+        * 
+        * This function calculates and returns the number of neighboring vacancies for a given
+        * lattice site, taking into account the periodic boundary conditions along all directions.
+        * It checks various diagonal and edge directions and communicates with neighboring processors
+        * if the vacancy is located at the boundary.
+        * 
+        * @param idx The index of the vacancy whose neighbors are to be counted.
+        * @return int The number of nearest neighbor vacancies.
         */
         int get_NN_count(int idx) {
             int i = vacancies_pos(idx,0); int j = vacancies_pos(idx,1); int k = vacancies_pos(idx,2); int l = vacancies_pos(idx,3);
@@ -223,7 +263,7 @@ class Lattice {
             std::vector<size_t> x_dims = proc_pos_x_neighbors.size_vec;
             std::vector<size_t> y_dims = proc_pos_y_neighbors.size_vec;
 
-            for (int s=0; s < (int)diag_directions.rows(); s++) {
+            for (int s=0; s < (int)diag_directions.size(); s++) {
 
                 if (i == 0) {
                     new_j = (((j - diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
@@ -236,7 +276,11 @@ class Lattice {
                     new_l = (((l + diag_directions[s][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]);
                 }
 
-                if ((i == 0) && (j == 0) && (diag_directions[s][0] == 1)) {/*communicate with proc to -x direction*/
+                if ((l == 0) && (i == 0) && (diag_directions[s][2] == 1)) {/* checking for leftmost non-periodic boundary along z-axis*/}
+                
+                else if ((l == (int)(sublattice_dim[2]-1)) && (i == 1) && (diag_directions[s][2] == 1)) {/* checking for rightmost non-periodic boundary along z-axis*/}
+                
+                else if ((i == 0) && (j == 0) && (diag_directions[s][0] == 1)) {/*communicate with proc to -x direction*/
                     
                     if ((k == 0) && (diag_directions[s][1] == 1)) {
                         /*check neighbor -x,-y array */
@@ -253,12 +297,12 @@ class Lattice {
                 else if ((i == 1) && (j == (sublattice_dim[0] - 1)) && (diag_directions[s][0] == 1)) {/*communicate with proc to +x direction*/ 
                     
                     if ((k == ((sublattice_dim[1] - 1))) && (diag_directions[s][1] == 1)) {
-                        /*check neighbor (+x,+y) array */
+                        /*check neighbor +x,+y array */
                         if (((proc_neighbors(rank,1) == rank)) && (!check_move_free(i,j,k,l,1,s,1))) {NN_count ++;}
                         else if (((proc_neighbors(rank,1) != rank)) && (proc_pos_x_neighbors(0, 0, (size_t)(x_dims[2]-1), (size_t)(new_l)))) {NN_count ++;}
                     }
                     else {
-                        /*check neighbor (+x) array */
+                        /*check neighbor +x array */
                         if (((proc_neighbors(rank,0) == rank)) && (!check_move_free(i,j,k,l,1,s,1))) {NN_count ++;}
                         else if (((proc_neighbors(rank,0) != rank)) && (proc_pos_x_neighbors(0, 0, (size_t)(new_k), (size_t)(new_l)))) {NN_count ++;}
                     }
@@ -267,18 +311,16 @@ class Lattice {
                 else if ((i == 0) && (k == 0) && (diag_directions[s][1] == 1)) {/*communicate with proc to -y direction*/
                     
                     if ((j == 0) && (diag_directions[s][0] == 1)) {
-                        /*check neighbor (-x,-y) array */
+                        /*check neighbor -x,-y array */
                         if ((proc_neighbors(rank,5) == rank) && (!check_move_free(i,j,k,l,-1,s,0))) {NN_count ++;}
                         else if ((proc_neighbors(rank,5) != rank) && (proc_neg_y_neighbors(0, 1, 0, (size_t)(new_l)))) {NN_count ++;}
                     }
-                    else {
-                        if (((proc_neighbors(rank,6) == rank)) && (!check_move_free(i,j,k,l,-1,s,0))) {NN_count ++;}
-                        else if (((proc_neighbors(rank,6) != rank)) && (proc_neg_y_neighbors(0, 1, (size_t)(new_j+1), (size_t)(new_l)))) {NN_count ++;}
-                    }
+                    else if (((proc_neighbors(rank,6) == rank)) && (!check_move_free(i,j,k,l,-1,s,0))) {NN_count ++;}
+                    else if (((proc_neighbors(rank,6) != rank)) && (proc_neg_y_neighbors(0, 1, (size_t)(new_j+1), (size_t)(new_l)))) {NN_count ++;}
                 }
 
                 else if ((i == 1) && (k == (sublattice_dim[1] - 1)) && (diag_directions[s][1] == 1)) {/*communicate with proc to +y direction*/
-
+                    
                     if ((j == ((sublattice_dim[0] - 1))) && (diag_directions[s][0] == 1)) {
                         /*check neighbor +x,+y array */
                         if ((proc_neighbors(rank,1) == rank) && (!check_move_free(i,j,k,l,1,s,1))) {NN_count ++;}
@@ -309,8 +351,16 @@ class Lattice {
             return NN_count;
         }
 
-        /*
-        find actions in given subdomain assigned to processor
+
+        /** 
+        * @brief Find actions in the subdomain assigned to the current processor.
+        * 
+        * This function loops over all vacancies in the system and identifies possible moves
+        * based on the nearest neighbor and edge directions. It performs communication between
+        * processors when necessary, handles boundary conditions, and resizes data structures 
+        * when required to store new moves.
+        * 
+        * The function also updates the cumulative sum of move rates.
         */
         void parallel_get_actions() {
             //std::cout << "rank: " << rank << " enter parallel get actions \n";
@@ -355,7 +405,7 @@ class Lattice {
                 NN_count = get_NN_count(idx);
 
                 // finding all moves along the {111} family of vectors
-                for (int s=0; s < (int)diag_directions.rows(); s++) {
+                for (int s=0; s < (int)diag_directions.size(); s++) {
                     if (i == 0) {
                         new_j = (((j - diag_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                         new_k = (((k - diag_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
@@ -438,7 +488,7 @@ class Lattice {
                 }
 
                 // finding all moves along the {100} family of vectors
-                for (int s=0; s < (int)edge_directions.rows(); s++) {
+                for (int s=0; s < (int)edge_directions.size(); s++) {
 
                     new_j = (((j + edge_directions[s][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]);
                     new_k = (((k + edge_directions[s][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]);
@@ -483,8 +533,29 @@ class Lattice {
 
         }
 
-        /*
-        function for finding rate corresponding to NN-encoding and move type
+
+        /**
+        * @brief Determines the rate constant corresponding to the nearest neighbor (NN) encoding 
+        *        and move type of a vacancy in the lattice.
+        * 
+        * This function computes the rate constant for a vacancy move based on the configuration 
+        * of the lattice, the coordinates of the vacancy, and the direction of the move. It handles 
+        * both vertex and boundary condition (bc) sites and checks for special regions with 
+        * predefined rate constants.
+        * 
+        * @param coord A pointer to an array of integers representing the coordinates of the vacancy 
+        *              in the lattice. The array has 4 elements: coord[0] to coord[3] are the lattice coordinates.
+        * @param shift A pointer to an array of integers representing the direction of the move.
+        *              The array has 3 elements corresponding to the movement along each axis.
+        * @param lattice An integer representing the type of move:
+        *                - 0: Moving vacancy from vertex site to bc site
+        *                - 1: Moving vacancy from bc site to vertex site
+        *                - 2 or 3: Moving vacancy from vertex site to vertex site, or bc site to bc site
+        * 
+        * @return A double representing the rate constant for the vacancy move. If the move is not allowed, 
+        *         the function returns -1.
+        * 
+        * @throws std::exception If the lattice type is invalid, an exception is thrown.
         */
         double new_get_rateconstants(int* coord, int* shift, int lattice) {  
             
@@ -549,8 +620,24 @@ class Lattice {
             return rate;
         }
 
-        /*
-        generating encoding corresponding to configuration of nearest neighbors for a vacancy in bulk of processor domain
+
+        /**
+        * @brief Generates the encoding corresponding to the configuration of nearest neighbors
+        *        for a vacancy in the bulk of the processor domain.
+        * 
+        * This function computes an integer corresponding to a binary encoding for configuration of 
+        * nearest neighbors around the vacancy's position
+        * 
+        * @param vac A pointer to an array of integers representing the vacancy's position in the lattice. 
+        *            The array has 4 elements: vac[0] is the type of site (0 or 1), and vac[1], vac[2], vac[3]
+        *            represent the coordinates of the vacancy in the lattice.
+        * @param lattice An integer representing the type of move:
+        *                - 3: Moving vacancy from bc site to bc site
+        *                - 2: Moving vacancy from vertex site to vertex site
+        *                - 1: Moving vacancy from bc site to vertex site
+        *                - 0: Moving vacancy from vertex site to bc site
+        * 
+        * @return An integer representing the sum of the encoding for nearest neighbor interactions.
         */
         int new_get_neighbors(int* vac, int lattice) {
             //
@@ -565,34 +652,34 @@ class Lattice {
 
             // moving vacancy from bc site to bc site 
             if (lattice == 3) {
-                for (int i=0; i<(int)diag_directions.rows(); i++) {
+                for (int i=0; i<(int)diag_directions.size(); i++) {
                     sum += exp_int(m,i) * vertex_sites(0, (((vac[0] - diag_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), (((vac[1] - diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), (((vac[2] - diag_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]));
                 }
-                for (int i=0; i<(int)edge_directions.rows(); i++) {
+                for (int i=0; i<(int)edge_directions.size(); i++) {
                     sum += exp_int(m,i) *  bc_sites(0, (((vac[0] + edge_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), (((vac[1] + edge_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), (((vac[2] + edge_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]));
                 }
             }
 
             // moving vacancy from vertex site to vertex site 
             else if (lattice == 2) {
-                for (int i=0; i<(int)diag_directions.rows(); i++) {
+                for (int i=0; i<(int)diag_directions.size(); i++) {
                     sum +=  exp_int(m,i) * bc_sites(0, (((vac[0] + diag_directions[i][0])  % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), (((vac[1] + diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), (((vac[2] + diag_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]));
                 }
-                for (int i=0; i<(int)edge_directions.rows(); i++) {
+                for (int i=0; i<(int)edge_directions.size(); i++) {
                     sum += exp_int(m,i) *  vertex_sites(0, (((vac[0] + edge_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), (((vac[1] + edge_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), (((vac[2] + edge_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]));
                 }
             }
 
             // moving vacancy from bc site to vertex site
             else if (lattice == 1) {
-                for (int i=0; i<(int)diag_directions.rows(); i++) {
+                for (int i=0; i<(int)diag_directions.size(); i++) {
                     sum += exp_int(m,i) * vertex_sites(0, (((vac[0] - diag_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), (((vac[1] - diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), (((vac[2] - diag_directions[i][2]) % sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]));
                 }
             }
 
             // moving vacancy from vertex site to bc site     
             else if (lattice == 0) {
-                for (int i=0; i<(int)diag_directions.rows(); i++) {
+                for (int i=0; i<(int)diag_directions.size(); i++) {
                     sum += exp_int(m,i) * bc_sites(0, (((vac[0] + diag_directions[i][0]) % sublattice_dim[0] + sublattice_dim[0]) % sublattice_dim[0]), (((vac[1] + diag_directions[i][1]) % sublattice_dim[1] + sublattice_dim[1]) % sublattice_dim[1]), (((vac[2] + diag_directions[i][2])% sublattice_dim[2] + sublattice_dim[2]) % sublattice_dim[2]));
                 }
             }
@@ -600,9 +687,29 @@ class Lattice {
             return sum;
         }
 
-        /* 
-        subroutine to check if move exceeds bounds of processor domain and send
-        information to adjacent processor
+
+        /**
+        * @brief Checks if a move exceeds the bounds of the processor domain and communicates the 
+        *        information to adjacent processors if necessary.
+        *
+        * This function verifies if a vacancy move crosses the boundary of the current processor domain. 
+        * If the move exceeds the domain, the function sends the updated information to adjacent processors 
+        * using MPI communication. The function handles different boundary conditions and directions of the move.
+        *
+        * @param i_old The old lattice layer index of the vacancy before the move.
+        * @param j_old The old x-coordinate of the vacancy before the move.
+        * @param k_old The old y-coordinate of the vacancy before the move.
+        * @param l_old The old z-coordinate of the vacancy before the move.
+        * @param move_idx The index of the move information .
+        * @param new_loc A constant reference to a vector of integers representing the new coordinates 
+        *                of the vacancy after the move. The vector contains the lattice layer indices, x, y, and z.
+        * 
+        * @return True if the move crosses the boundary of the processor domain and information is sent 
+        *         to adjacent processors, otherwise False.
+        * 
+        * @note The function uses MPI non-blocking communication (MPI_Isend) to send data to neighboring 
+        *       processors. It handles multiple ghost regions and different data structures for various 
+        *       lattice configurations.
         */
         bool parallel_processes_check(int i_old, int j_old, int k_old, int l_old, int move_idx, const std::vector<int>& new_loc) {
             //
@@ -634,7 +741,7 @@ class Lattice {
             // (111) moves
             if ((moves_lattice(move_idx,0) == 0) || (moves_lattice(move_idx,0) == 1)) {
                 if ((j_old == 0) && (moves_shifts(move_idx,0) == -1)) {/*communicate with proc to -x direction*/
-                    if ((k_old == 0) && (moves_shifts(move_idx,1) == -1)) { /*communicate with proc to (-x,-y) direction*/
+                    if ((k_old == 0) && (moves_shifts(move_idx,1) == -1)) {
                         
                         if ((proc_neighbors(rank,6) == proc_neighbors(rank,5)) && (proc_neighbors(rank,4) == proc_neighbors(rank,5))) {}
                         
@@ -654,8 +761,7 @@ class Lattice {
                                 1,                  //Message Tag
                                 MPI_COMM_WORLD,      //MPI Communicator
                                 &request1 ); 
-
-                            MPI_Wait(&request1, MPI_STATUS_IGNORE);
+                            
                         }
                         else if ((proc_neighbors(rank,4) == proc_neighbors(rank,5)) && (proc_neighbors(rank,6) == rank)) {
                             new_proc = proc_neighbors(rank,4);
@@ -674,8 +780,6 @@ class Lattice {
                                 MPI_COMM_WORLD,      //MPI Communicator
                                 &request1 ); 
 
-                            MPI_Wait(&request1, MPI_STATUS_IGNORE);
-
                         }
                         else {
                             new_proc = proc_neighbors(rank,5);
@@ -686,14 +790,13 @@ class Lattice {
                             loc_buffer2[7] = k_old;
                             loc_buffer2[8] = l_old;
                             MPI_Isend(loc_buffer2.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request2);
-
-                            MPI_Wait(&request1, MPI_STATUS_IGNORE);
                         }
                         proc_neg_y_neighbors(0, i, (size_t)0, l) = 1;
                         proc_neg_x_neighbors(0, i, (size_t)0, l) = 1;
                     }
-                    else { /*communicate with proc to (-x) direction*/
-                        
+                    else {
+                        /*communicate with proc to (-x) direction*/
+                        //std::cout << "rank: " << rank << " communicate with proc to (-x) direction \n\n";
                         new_proc = proc_neighbors(rank,4); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y)
 
                         if (new_proc == rank) { return false; }
@@ -707,8 +810,7 @@ class Lattice {
                             loc_buffer1[8] = l_old; 
                             ////std::cout << "rank: " << rank << " sending to proc: " << new_proc <<  "\n";
                             MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1); 
-                            
-                            MPI_Wait(&request1, MPI_STATUS_IGNORE);
+                            //std::cout << "rank: " << rank << " sent to proc: " << new_proc <<  "\n";
                         }
                         proc_neg_x_neighbors(0, i, k+1, l) = 1;
                     }
@@ -716,7 +818,8 @@ class Lattice {
                     return true;
                 }
                 else if ((j_old == (sublattice_dim[0] - 1)) && (moves_shifts(move_idx,0) == 1)) {/*communicate with proc to +x direction*/
-                    if ((k_old == (sublattice_dim[1] - 1))  && (moves_shifts(move_idx,1) == 1)) { /*communicate with proc to (+x,+y) direction*/
+                    if ((k_old == (sublattice_dim[1] - 1))  && (moves_shifts(move_idx,1) == 1)) {
+                        /*communicate with proc to (+x,+y) direction*/
 
                         if ((proc_neighbors(rank,2) == proc_neighbors(rank,1)) && (proc_neighbors(rank,2) == proc_neighbors(rank,0))) {return false;}
                         
@@ -736,8 +839,6 @@ class Lattice {
                                 1,                  //Message Tag
                                 MPI_COMM_WORLD,      //MPI Communicator
                                 &request1 ); 
-
-                            MPI_Wait(&request1, MPI_STATUS_IGNORE);
                         }
                         else if ((proc_neighbors(rank,0) == proc_neighbors(rank,1)) && (proc_neighbors(rank,2) == rank)) {
                             new_proc = proc_neighbors(rank,0);
@@ -756,8 +857,7 @@ class Lattice {
                                 1,                  //Message Tag
                                 MPI_COMM_WORLD,      //MPI Communicator
                                 &request1 ); 
-                                
-                            MPI_Wait(&request1, MPI_STATUS_IGNORE);
+                            //std::cout << "rank: " << rank << " sent to proc: " << new_proc <<  "\n";
                         }
                         else {
                             new_proc = proc_neighbors(rank,1);
@@ -767,15 +867,15 @@ class Lattice {
                             loc_buffer2[6] = j_old;
                             loc_buffer2[7] = k_old;
                             loc_buffer2[8] = l_old;
-                            MPI_Isend(loc_buffer2.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);
-
-                            MPI_Wait(&request1, MPI_STATUS_IGNORE);                            
+                            MPI_Isend(loc_buffer2.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request2);
+                            
                         }
 
                         proc_pos_x_neighbors(0, i, (size_t)(x_dims[2]-1), (size_t)(l)) = 1;
                         proc_pos_y_neighbors(0, i, (size_t)(y_dims[2]-1), (size_t)(l)) = 1;
                     }
-                    else { /*communicate with proc to (+x) direction*/
+                    else {
+                        /*communicate with proc to (+x) direction*/
                         
                         new_proc = proc_neighbors(rank,0); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y) 
                         if (new_proc == rank) {return false;}
@@ -787,15 +887,15 @@ class Lattice {
                             loc_buffer1[7] = k_old;
                             loc_buffer1[8] = l_old;
                             MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);
-
-                            MPI_Wait(&request1, MPI_STATUS_IGNORE);
                         }
                         proc_pos_x_neighbors(0, i, (size_t)(k), (size_t)(l)) = 1;
                     }
 
                     return true;
                 }
-                else if ((k_old == 0) && (moves_shifts(move_idx,1) == -1)) {/*communicate with proc to -y direction*/ 
+                else if ((k_old == 0) && (moves_shifts(move_idx,1) == -1)) {
+                    //std::cout << "rank: " << rank << " communicate with proc to (-y) direction \n\n";
+                    /*communicate with proc to -y direction*/ 
                     new_proc = proc_neighbors(rank,6); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y) 
 
                     if (new_proc == rank) { return false; }
@@ -806,16 +906,15 @@ class Lattice {
                         loc_buffer1[6] = j_old;
                         loc_buffer1[7] = k_old;
                         loc_buffer1[8] = l_old; 
-                        
+                        //std::cout << "rank: " << rank << " sending to proc: " << new_proc <<  "\n";
                         MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);
-
-                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
-                        
+                        //std::cout << "rank: " << rank << " sent to proc: " << new_proc <<  "\n";
                     }
                     proc_neg_y_neighbors(0, i, (size_t)(j+1), (size_t)(l)) = 1;
                     return true;
                 }
-                else if ((k_old == (sublattice_dim[1] - 1)) && (moves_shifts(move_idx,1) == 1)) { /*communicate with proc to +y direction*/
+                else if ((k_old == (sublattice_dim[1] - 1)) && (moves_shifts(move_idx,1) == 1)) {
+                    /*communicate with proc to +y direction*/
                     new_proc = proc_neighbors(rank,2); // proc_neighbors indices follow: 0:+x, 1:(+x,+y), 2:+y, 3:(-x,+y), 4:-x, 5:(-x,-y), 6:-y, 7:(+x,-y) 
 
                     if (new_proc == rank) { return false; }
@@ -828,8 +927,6 @@ class Lattice {
                         loc_buffer1[8] = l_old;
                         
                         MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);
-
-                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
                     }
 
                     proc_pos_y_neighbors(0, i, (size_t)(j), (size_t)(l)) = 1;
@@ -853,7 +950,6 @@ class Lattice {
                     loc_buffer1[8] = l_old; 
                     MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);
 
-                    MPI_Wait(&request1, MPI_STATUS_IGNORE);
                     return true;
                 }
                 else if ((moves_shifts(move_idx,0) == -1) && (j_old == 0)) {
@@ -871,7 +967,6 @@ class Lattice {
                     loc_buffer1[8] = l_old; 
                     MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);
 
-                    MPI_Wait(&request1, MPI_STATUS_IGNORE);
                     return true;
                 }
                 else if ((moves_shifts(move_idx,1) == 1) && (k_old == (sublattice_dim[1] - 1))) {
@@ -889,7 +984,6 @@ class Lattice {
                     loc_buffer1[8] = l_old; 
                     MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);
 
-                    MPI_Wait(&request1, MPI_STATUS_IGNORE);
                     return true;
                 }
                 else if ((moves_shifts(move_idx,1) == -1) && (k_old == 0)) {
@@ -907,7 +1001,6 @@ class Lattice {
                     loc_buffer1[8] = l_old; 
                     MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, 1, MPI_COMM_WORLD, &request1);
 
-                    MPI_Wait(&request1, MPI_STATUS_IGNORE);
                     return true;
                 }
             }
@@ -915,8 +1008,23 @@ class Lattice {
             return false;
         }
 
-        /*
-        method to send information to other processor regarding ghost sites
+
+        /**
+        * @brief Sends information about ghost sites to other processors.
+        *
+        * This method communicates the status of ghost sites to neighboring processors.
+        * It performs checks to ensure that the old and new locations are within valid ranges.
+        * The data includes the indices of the ghost sites and the rank of the sending processor.
+        * Depending on the values of the input parameters, it sends the data to the appropriate neighboring processor
+        * using MPI_Isend for parallel transfer.
+        *
+        * @param i_old Old lattice coordinate (vertex/bc site) of site in ghost zone.
+        * @param j_old Old x-coordinate of the ghost site.
+        * @param k_old Old y-coordinate of the ghost site.
+        * @param l_old Old z-coordinate of the ghost site.
+        * @param idx Index of the ghost site (used for data management).
+        * @param new_loc Vector containing the new locations of the ghost sites.
+        * @param parallel_transfer Boolean flag indicating if the transfer is parallel.
         */
         void ghost_site_send(int i_old, int j_old, int k_old, int l_old, int idx, const std::vector<int>& new_loc, bool parallel_transfer) {         
             //
@@ -950,7 +1058,7 @@ class Lattice {
             if (j_old < 1) {/*communicate with proc to -x direction*/
                 if (k_old < 1) {
                     /*communicate with proc to (-x,-y) direction*/
-                    if ((proc_neighbors(rank,6) == proc_neighbors(rank,5)) && (proc_neighbors(rank,4) == proc_neighbors(rank,5))) {}//std::cout << "rank: " << rank << " pass (-x,-y)\n\n";}     
+                    if ((proc_neighbors(rank,6) == proc_neighbors(rank,5)) && (proc_neighbors(rank,4) == proc_neighbors(rank,5))) {std::cout << "rank: " << rank << " pass (-x,-y)\n\n";}     
                     else if ((proc_neighbors(rank,6) == rank) && (proc_neighbors(rank,4) != rank)) {
                         new_proc = proc_neighbors(rank,4);
                         for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
@@ -967,8 +1075,6 @@ class Lattice {
                             tag,                  //Message Tag
                             MPI_COMM_WORLD,      //MPI Communicator
                             &request1 ); 
-
-                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
                     }
                     else if ((proc_neighbors(rank,4) == rank) && (proc_neighbors(rank,6) != rank)) {
                         new_proc = proc_neighbors(rank,6);
@@ -986,8 +1092,6 @@ class Lattice {
                             tag,                  //Message Tag
                             MPI_COMM_WORLD,      //MPI Communicator
                             &request1 ); 
-
-                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
                     }
                     else if ((proc_neighbors(rank,6) != rank) && (proc_neighbors(rank,4) != rank)) {
                         new_proc = proc_neighbors(rank,6);
@@ -1004,8 +1108,7 @@ class Lattice {
                             new_proc,           //Rank of receiving process
                             tag,                  //Message Tag
                             MPI_COMM_WORLD,      //MPI Communicator
-                            &request1 );                 
-                        MPI_Wait(&request1, MPI_STATUS_IGNORE);            
+                            &request1 );                             
                         new_proc = proc_neighbors(rank,5);
                         for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer2[idx] = new_loc[idx];}
                         loc_buffer2[4] = rank;
@@ -1014,7 +1117,6 @@ class Lattice {
                         loc_buffer2[7] = k_old;
                         loc_buffer2[8] = l_old;
                         MPI_Isend(loc_buffer2.data(), bufferlen, MPI_INT, new_proc, tag, MPI_COMM_WORLD, &request2);
-                        MPI_Wait(&request2, MPI_STATUS_IGNORE);
                         new_proc = proc_neighbors(rank,4);
                         for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer3[idx] = new_loc[idx];}
                         loc_buffer3[4] = rank;
@@ -1023,7 +1125,6 @@ class Lattice {
                         loc_buffer3[7] = k_old;
                         loc_buffer3[8] = l_old;
                         MPI_Isend(loc_buffer3.data(), bufferlen, MPI_INT, new_proc, tag, MPI_COMM_WORLD, &request3);
-                        MPI_Wait(&request3, MPI_STATUS_IGNORE);
                     }
                 }
                 else {
@@ -1039,15 +1140,13 @@ class Lattice {
                         loc_buffer1[7] = k_old;
                         loc_buffer1[8] = l_old;
                         MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, tag, MPI_COMM_WORLD, &request1);
-
-                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
                     }
                 }
             }
             else if (j_old > (sublattice_dim[0] - 2)) {/*communicate with proc to +x direction*/
                 if (k_old > (sublattice_dim[1] - 2)) {
                     /*communicate with proc to (+x,+y) direction*/
-                    if ((proc_neighbors(rank,2) == proc_neighbors(rank,1)) && (proc_neighbors(rank,2) == proc_neighbors(rank,0))) {}//std::cout << "rank: " << rank << " pass (+x,+y)\n\n";}
+                    if ((proc_neighbors(rank,2) == proc_neighbors(rank,1)) && (proc_neighbors(rank,2) == proc_neighbors(rank,0))) {std::cout << "rank: " << rank << " pass (+x,+y)\n\n";}
                     else if ((proc_neighbors(rank,0) == rank) && (proc_neighbors(rank,2) != rank)) {
                         new_proc = proc_neighbors(rank,2);
                         for (int idx=0; idx<new_loc.size(); idx++) {loc_buffer1[idx] = new_loc[idx];}
@@ -1064,8 +1163,6 @@ class Lattice {
                             tag,                  //Message Tag
                             MPI_COMM_WORLD,      //MPI Communicator
                             &request1 ); 
-
-                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
                     }
                     else if ((proc_neighbors(rank,2) == rank) && (proc_neighbors(rank,0) != rank)) {
                         new_proc = proc_neighbors(rank,0);
@@ -1083,8 +1180,6 @@ class Lattice {
                             tag,                  //Message Tag
                             MPI_COMM_WORLD,      //MPI Communicator
                             &request1 ); 
-
-                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
                     }
                     else if ((proc_neighbors(rank,2) != rank) &&  (proc_neighbors(rank,0) != rank)) {
                         new_proc = proc_neighbors(rank,2);
@@ -1103,7 +1198,6 @@ class Lattice {
                             MPI_COMM_WORLD,      //MPI Communicator
                             &request1 ); 
                         new_proc = proc_neighbors(rank,1);
-                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
                         for (int idx=0; idx< new_loc.size(); idx++) {loc_buffer2[idx] = new_loc[idx];}
                         loc_buffer2[4] = rank;
                         loc_buffer2[5] = i_old;
@@ -1111,7 +1205,6 @@ class Lattice {
                         loc_buffer2[7] = k_old;
                         loc_buffer2[8] = l_old;
                         MPI_Isend(loc_buffer2.data(), bufferlen, MPI_INT, new_proc, tag, MPI_COMM_WORLD, &request2);
-                        MPI_Wait(&request2, MPI_STATUS_IGNORE);
                         new_proc = proc_neighbors(rank,0);
                         for (int idx=0; idx< new_loc.size(); idx++) {loc_buffer3[idx] = new_loc[idx];}
                         loc_buffer3[4] = rank;
@@ -1120,7 +1213,6 @@ class Lattice {
                         loc_buffer3[7] = k_old;
                         loc_buffer3[8] = l_old;
                         MPI_Isend(loc_buffer3.data(), bufferlen, MPI_INT, new_proc, tag, MPI_COMM_WORLD, &request3);
-                        MPI_Wait(&request3, MPI_STATUS_IGNORE);
                     }
                 }
                 else {
@@ -1135,7 +1227,6 @@ class Lattice {
                         loc_buffer1[7] = k_old;
                         loc_buffer1[8] = l_old;
                         MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, tag, MPI_COMM_WORLD, &request1);
-                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
                     }
                 }
             }
@@ -1152,7 +1243,6 @@ class Lattice {
                     loc_buffer1[7] = k_old;
                     loc_buffer1[8] = l_old; 
                     MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, tag, MPI_COMM_WORLD, &request1);
-                    MPI_Wait(&request1, MPI_STATUS_IGNORE);
                 }
             }
             else if ( (k_old > (sublattice_dim[1] - 2)) ) {/*communicate with proc to +y direction*/
@@ -1167,13 +1257,30 @@ class Lattice {
                     loc_buffer1[7] = k_old;
                     loc_buffer1[8] = l_old;
                     MPI_Isend(loc_buffer1.data(), bufferlen, MPI_INT, new_proc, tag, MPI_COMM_WORLD, &request1);
-                    MPI_Wait(&request1, MPI_STATUS_IGNORE);
                 }
             }
         }
 
-        /*
-        method to recieve information from other processor regarding ghost sites
+
+        /**
+        * @brief Handles the reception of ghost site information from neighboring processes.
+        *
+        * This function updates the neighbor process lists based on the new location of a ghost site 
+        * received in the `new_loc_buffer`. It accounts for both parallel and non-parallel transfers 
+        * and modifies the appropriate neighbors according to the ghost site's old and new positions.
+        *
+        * @param new_loc_buffer A vector containing the new location and other related data of the ghost site.
+        *                       Expected format:
+        *                       - new_loc_buffer[0]: New lattice coordinate (vertex/bc site)
+        *                       - new_loc_buffer[1]: New x-coordinate
+        *                       - new_loc_buffer[2]: New y-coordinate
+        *                       - new_loc_buffer[3]: New z-coordinate
+        *                       - new_loc_buffer[4]: Old processor ID
+        *                       - new_loc_buffer[5]: Old lattice coordinate (vertex/bc site)
+        *                       - new_loc_buffer[6]: Old x-coordinate
+        *                       - new_loc_buffer[7]: Old y-coordinate
+        *                       - new_loc_buffer[8]: Old z-coordinate
+        * @param parallel_transfer A boolean indicating whether the transfer is parallel.
         */
         void ghost_site_recieve(const std::vector<int>& new_loc_buffer, bool parallel_transfer) {
             //std::cout << "rank: " << rank << " ghost_site_recieve()\n\n";
@@ -1415,8 +1522,26 @@ class Lattice {
             
         }
 
-        /*
-        method to update ghost sites in local domain of processor
+
+        /**
+        * @brief Update ghost sites in the local domain of the processor.
+        * 
+        * This function checks the current location of ghost sites and updates 
+        * their references based on the new location provided. It distinguishes 
+        * between parallel and non-parallel transfer scenarios, adjusting 
+        * neighbor references accordingly. The function uses the global rank 
+        * of the processor to determine which ghost sites need to be updated 
+        * and applies the changes in their respective positions based on the 
+        * defined neighbor relationships.
+        * 
+        * @param i_old Old lattice coordinate (vertex/bc site) of site in ghost zone.
+        * @param j_old Old x-coordinate of the ghost site.
+        * @param k_old Old y-coordinate of the ghost site.
+        * @param l_old Old z-coordinate of the ghost site.
+        * @param new_loc A vector containing the new location coordinates: 
+        *                {lattice type, x, y, z}.
+        * @param parallel_transfer Boolean flag indicating whether the update 
+        *                         occurs in parallel.
         */
         void ghost_site_self_reference(int i_old,int j_old,int k_old,int l_old, const std::vector<int>& new_loc, bool parallel_transfer) {
             //std::cout << "rank: " << rank << " enter ghost_site_self_reference\n";
@@ -1515,17 +1640,25 @@ class Lattice {
             //std::cout << "rank: " << rank << " exit ghost_site_self_reference\n";
         }
 
-        /*
-        updating the positions of atoms on lattice according to selected move
+
+        /**
+        * @brief Updates the positions of atoms on the lattice according to the selected move.
+        * 
+        * This function updates the positions of vacancies and atoms in the lattice based on the 
+        * selected move for a given index. It handles both parallel and non-parallel transfer cases 
+        * while ensuring that the boundaries of the domain are respected.
+        * 
+        * @param idx The index of to access information about move being propagated from moves_vacs,
+        *               moves_lattice, and other Matrix data structures
+        * @param move_ticks The number of move that is being propagated
+        * 
+        * @note This function checks the boundaries of the domain against the total dimensions. 
+        *       If equal and the move is within ghost sites, a modified move is created to find 
+        *       a new location in ghost sites.
         */
-        //UPDATE TO CHECK BOUNDS OF DOMAIN VS TOTAL_DIMS: IF EQUAL & MOVE WITHIN GHOST SITES, THEN CREATE 
-        // MODDED MOVE TO FIND NEW LOCATION IN GHOST SITES
         void new_update_lattice(int idx, int move_ticks) {
 
             bool parallel_transfer = false;
-            //std::cout << "lattice rank: " << rank << " moves_lattice[idx][0]: " << moves_lattice(idx,0) << "\n";
-            //std::cout << "lattice rank: " << rank << " vacancies_pos.rows(): " << vacancies_pos.rows() << "\n";
-
             if (moves_lattice(idx,0) == 5) { 
                 prev_move_type.push_back(0);
                 prev_move_type_ticks.push_back(move_ticks);
@@ -1719,29 +1852,35 @@ class Lattice {
             }
         }
 
-        /*
-        calculating time elapsed for a move
+
+        /**
+        * @brief Calculates the time elapsed for a move.
+        * 
+        * This function generates a random double between 0 and 1, 
+        * then calculates the elapsed time based on the cumulative rate.
+        * 
+        * @return double The calculated time elapsed for a move.
         */
         double new_random_times() {
             // creating a random double between 0 and 1
-            double time = 0;
-            if (rank == 0) {
-                unsigned int random = mt_obj();
-                double random_double = ((random / (1.+ UINT32_MAX)) +  (1 / (1.+ UINT32_MAX)));
-                
-                //calculating the time elapsed
-                int last_idx = (int)rate_cumsum.size() - 1;
-                time = ((-1/ rate_cumsum[last_idx]) * log(random_double));
-            }
-
-            MPI_Bcast(&time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
+            unsigned int random = mt_obj();
+            double random_double = ((random / (1.+ UINT32_MAX)) +  (1 / (1.+ UINT32_MAX)));
+            
+            //calculating the time elapsed
+            int last_idx = (int)rate_cumsum.size() - 1;
+            double time = ((-1/ rate_cumsum[last_idx]) * log(random_double));
             
             return time;
         }
 
-        /*
-        communicating total rate in each processor domain at each timestep and creating null move
-        corresponding to difference between Rtot_i and Rmax
+
+        /**
+        * @brief Communicates the total rate in each processor domain at each timestep.
+        * 
+        * This function gathers the maximum rate across all processes, 
+        * broadcasts it to all processes, and creates a null move corresponding 
+        * to the difference between the total rate (Rtot_i) and the maximum rate (Rmax).
+        * 
         */
         void communicate_rates() {
             double max_i_rate = 0;
@@ -1752,19 +1891,28 @@ class Lattice {
             MPI_Request request;
             int end_idx;
             if (rate_cumsum.size() != 0) { max_i_rate = rate_cumsum[rate_cumsum.size()-1]; }
-            else { max_i_rate = 0; }  
+            else { max_i_rate = 0; }
+
+            MPI_Gather(&max_i_rate, 1, MPI_DOUBLE, max_rates.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
+                        
+            if (rank == 0) {
+                max_rate_idx = find_max_element(&max_rates);
+               //std::cout << "rank: " << rank << " max_rate_idx " << max_rate_idx << "\n"; 
+               //std::cout << "rank: " << rank << " max_rates \n";
+                //print_1Dvector(max_rates);
+                max_rate = max_rates[max_rate_idx];
+            }
             
-            //std::cout << "rank: " << rank << " max_i_rate: " << max_i_rate << "\n";
-            MPI_Reduce(&max_i_rate, &max_rate, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
             MPI_Bcast(&max_rate, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
             
             if (rate_cumsum.size() != 0 ) {
-                if (rate_cumsum[(rate_cumsum.size() - 1)] == max_rate) {//std::cout << "rank: " << rank << " has the max rate! \n"; 
+               //std::cout << "rank: " << rank << " rate_cumsum[(rate_cumsum.size() - 1)] " << rate_cumsum[(rate_cumsum.size() - 1)]  << "\n";
+                if (rate_cumsum[(rate_cumsum.size() - 1)] == max_rate) {//std::cout << "rank: " << rank << " has the max rate! \n";
                 }
                 else {
                     fflush(stdout);
                     end_idx = (int)moves_lattice.rows() + 1;
-                    //std::cout << "rank: " << rank << " doesn't have the max rate - end_idx: " << end_idx << "\n";
+                   //std::cout << "rank: " << rank << " doesn't have the max rate - end_idx: " << end_idx << "\n";
                     moves_lattice.reshape(end_idx, 1, rank);
                     moves_lattice((end_idx-1),0) = 5;
                     rate_cumsum.push_back(max_rate);
@@ -1787,12 +1935,20 @@ class Lattice {
                 moves_vacs(0,0) = -1;
                 moves_coords.reshape(1,4, rank);
                 for (int i=0; i<4; i++) moves_coords(0, i) = -1;
-            }                  
+            }
+            
         }
 
-        /*
-        selecting random move in vector of moves, with selection probability proportional
-        to rate constant corresponding to move
+
+        /**
+        * @brief Selects a random move from a vector of moves.
+        * 
+        * This function selects a random index in the vector of moves, 
+        * with selection probability proportional to the rate constant 
+        * corresponding to each move. It communicates the current rates 
+        * and generates a random number to determine the selected index.
+        * 
+        * @return int The index of the selected move.
         */
         int get_idx() {
             //
@@ -1816,31 +1972,33 @@ class Lattice {
             }
             //std::cout << "rank: " << rank << " random_idx: " << min_idx << " last_idx: " << last_idx << " rate_cumsum[last_idx]: " << rate_cumsum[last_idx] << "\n";
 
-            //MPI_Barrier(MPI_COMM_WORLD); 
-
-
             return min_idx;
         }
 
-        /* 
-        mpi routine to tell other processes to roll back two steps and restart run 
+
+        /**
+        * @brief Communicates a boundary conflict to other processes.
+        * 
+        * This MPI routine informs all other processes to roll back two steps 
+        * and restart the simulation run. It uses a non-blocking send to notify 
+        * other ranks of the rollback action.
         */
         void comm_boundary_conflict() {
             int tag = 6;
             MPI_Request request;
         
-            for (int new_proc=0; new_proc<num_procs; new_proc++) { 
-                if (new_proc != rank) {
-                    MPI_Isend(NULL, 0, MPI_CHAR, new_proc, tag, MPI_COMM_WORLD, &request ); 
-                    MPI_Wait(&request, MPI_STATUS_IGNORE);
-                }    
-            }
+            for (int new_proc=0; new_proc<num_procs; new_proc++) { if (new_proc != rank) MPI_Isend(NULL, 0, MPI_CHAR, new_proc, tag, MPI_COMM_WORLD, &request ); }
         }
 
 
-        /* 
-        truncating vector of previous parallel moves to only include those
-        from previous two iterative steps
+        /**
+        * @brief Truncates the vector of previous parallel moves.
+        * 
+        * This function removes old parallel moves from par_prev_moves, retaining 
+        * only those moves from the previous two iterative steps. It updates 
+        * par_prev associated vectors to ensure they reflect the current state of moves.
+        * 
+        * @param move_ticks The current number move in simulation
         */
         void remove_old_par_moves(int move_ticks) {
             
@@ -1877,10 +2035,14 @@ class Lattice {
             }
         }
 
-        /* 
-        rolling back move originating in native process when error encountered 
-        in simulation - one case for involving parallel transfer, one for no
-        parallel transfer
+
+        /** 
+        * @brief Rolls back a move originating in the native process when an error is encountered in the simulation. 
+        * 
+        * This function handles two scenarios: one involving parallel transfer and the other without. 
+        * It updates the relevant prev_moves vectors based on the previous state of the system. 
+        * 
+        * @param parallel_transfer A boolean indicating whether the move involved parallel transfer.
         */
         void reverse_move(bool parallel_transfer) {
             std::vector<int> new_vac;
@@ -1891,7 +2053,6 @@ class Lattice {
             int i1; int i2; int i3; int i4;
             int vacs_idx;
  
-            /*
             std::cout << "rank: " << rank << " prev_newlocs: \n";
             print_2Dvector(prev_newlocs);
             std::cout << "rank: " << rank << " prev_oldlocs \n";
@@ -1900,7 +2061,6 @@ class Lattice {
             print_2Dvector(prev_moves);
             std::cout << "rank: " << rank << " prev_lattice: \n";
             print_1Dvector(prev_lattice); 
-            */
                
             new_vac = prev_newlocs.at((prev_newlocs.size() - 1));
             old_vac = prev_oldlocs.at((prev_oldlocs.size() - 1));
@@ -2093,11 +2253,21 @@ class Lattice {
             prev_lattice.pop_back();
             prev_idxs.pop_back(); 
                     
-        }        
-
-        /* 
-        rolling back moves involving parallel communication when error 
-        encountered in simulation
+        }      
+
+
+        /**
+        * @brief Rolls back moves involving parallel communication when an error 
+        * encountered in simulation.
+        *
+        * This function reverses the most recent move by updating the occupancy of
+        * lattice sites and adjusting the status of vacancies. It handles restoring 
+        * the original occupancy in the native processor - the removal of vacancy
+        * from adjacent processor domain is handled by reverse_move() in adjacent 
+        * rank process 
+        *
+        * It also manages the rollback of associated data structures used for 
+        * tracking previous moves.
         */
         void reverse_move_parallel() {
             int vacs_idx;
@@ -2106,10 +2276,18 @@ class Lattice {
             
             vacs_idx = par_prev_idx.at(par_prev_idx.size()-1),0;
             
+            std::cout << "rank: " << rank << " parallel_prev_moves: \n";
+            print_2Dvector(par_prev_moves);
+
             i1 = (par_prev_moves.at((par_prev_moves.size()-1)).at(0));
             i2 = (par_prev_moves.at((par_prev_moves.size()-1)).at(1));
             i3 = (par_prev_moves.at((par_prev_moves.size()-1)).at(2));
             i4 = (par_prev_moves.at((par_prev_moves.size()-1)).at(3));
+            
+            std::cout << "rank: " << rank << " parallel site (set 0): " << "[ " << i1 << " " << i2 << " " << i3 << " " << i4 << " ]\n";
+            std::cout << "rank: " << rank << " parallel old site status: " << vacancies(i1,i2,i3,i4) << " \n"; 
+            std::cout << "rank: " << rank << " parallel vacancies_pos(vacs_idx): [ " << vacancies_pos(vacs_idx,0) << " " << vacancies_pos(vacs_idx,1) << " " << vacancies_pos(vacs_idx,2) << " " << vacancies_pos(vacs_idx,3) << " ]\n"; 
+
 
             
             if (i1 == 1) {
@@ -2131,21 +2309,39 @@ class Lattice {
             par_prev_idx.pop_back();
         }
 
-        /* 
-        rolling back clock associated with previous moves when error encountered
-        in simulation
+
+        /**
+        * @brief Decrements the simulation clock based on previously stored time values.
+        *
+        * This function rolls back the simulation clock by subtracting the recorded
+        * time intervals from the previous moves. It also clears the list of previous 
+        * time values after adjusting the clock.
         */
         void deincrement_time() {
-            //print_1Dvector(prev_times);
+            print_1Dvector(prev_times);
             for (int i=0; i < (int)prev_times.size(); i++) {
                 t = t - prev_times[i];}
             
             prev_times.clear();
         }
 
-        /*
-        store information for move originating in native process - two different storage
-        types depending on whether or not move involved interprocessor communication
+
+        /**
+        * @brief Stores information about a move originating in the native process.
+        *
+        * This function manages the storage of move information based on whether the
+        * move involved interprocessor communication. It updates the locations of 
+        * vacancies and tracks the shifts made during the move. Two different storage
+        * paths are utilized depending on whether the move is a parallel transfer or 
+        * a standard move.
+        *
+        * @param idx The index of the move being stored.
+        * @param parallel_transfer Indicates if the move involves parallel transfer.
+        * @param vac_idx The index of the vacancy being moved.
+        * @param lattice Optional parameter for lattice type, (default is -1).
+        * @param shift A vector representing the shift of the move (default is {0,0,0}).
+        * @param new_loc A vector representing the new location of the vacancy (default is {0,0,0,0}).
+        * @param old_loc A vector representing the old location of the vacancy (default is {0,0,0,0}).
         */
         void store_move_info(int idx,  bool parallel_transfer, int vac_idx, int lattice = -1, std::vector<int> shift = {0,0,0}, const std::vector<int>& new_loc = {0,0,0,0}, const std::vector<int>& old_loc = {0,0,0,0}) {
 
@@ -2202,9 +2398,17 @@ class Lattice {
                 prev_lattice.push_back(lattice);
             }
         }
+        
 
-        /*
-        store information for move involving parallel transfer from another process
+        /**
+        * @brief Stores information about a move involving parallel transfer from another process.
+        *
+        * This function captures details of a move that was transferred in parallel, 
+        * updating the necessary vacancy and move tracking structures.
+        *
+        * @param parallel_buffer A vector containing the details of the vacancy involved in the move.
+        * @param move_ticks The number of ticks associated with the move.
+        * @param vac_idx The index of the vacancy that was moved.
         */
         void store_parallel_info(const std::vector<int>& parallel_buffer, int move_ticks, int vac_idx) {
             std::vector<int> vac(4); 
@@ -2218,8 +2422,15 @@ class Lattice {
             par_prev_idx.push_back(vac_idx);
         }
 
-        /*
-        store time elapsed during previous iteration of KMC process
+
+        /**
+        * @brief Stores the time elapsed during the previous iteration of the KMC process.
+        *
+        * This function updates the time tracking structures by storing the time increment 
+        * from the last Kinetic Monte Carlo (KMC) process iteration. If there are 
+        * already two previous time entries, the oldest one is overwritten.
+        *
+        * @param timestep The time increment to be stored.
         */
         void store_time_incr(double timestep) {
             if ((int)prev_times.size() >= 2) {
@@ -2228,12 +2439,15 @@ class Lattice {
             }          
         }
 
-        /*
-        helper function to reverse previous two moves depending on whether 
-        move involves parallel transfer or not
+
+        /**
+        * @brief Helper function to reverse the previous two moves.
+        *
+        * This function checks the type of each previous move (in-lattice or parallel transfer)
+        * and calls the appropriate reverse function to undo the moves. After reversing,
+        * it clears the recorded previous move types and ticks.
         */
         void reverse_moves_wrapper() {
-            /*
             std::cout << "rank: " << rank << " reverse_moves_wrapper: \n";
             std::cout << "rank: " << rank << " prev_move_type: \n";
             print_1Dvector(prev_move_type);
@@ -2241,10 +2455,9 @@ class Lattice {
             print_1Dvector(prev_move_type_ticks);
             std::cout << "rank: " << rank << " prev_newlocs: \n";
             print_2Dvector(prev_newlocs);
-            */
 
             for (int i=(prev_move_type.size()-1); i>=0; i--) { 
-                //std::cout << "rank: " << rank << " prev_move_type[i]: " << prev_move_type[i] << "\n";
+                std::cout << "rank: " << rank << " prev_move_type[i]: " << prev_move_type[i] << "\n";
                 if (prev_move_type[i] == 0) reverse_move(false); // reversing in-lattice move with no parallel transfer
                 if (prev_move_type[i] == 1) reverse_move(true); // reversing move with sending parallel transfer
                 if (prev_move_type[i] == 2) reverse_move_parallel(); // reversing move with recieving parallel transfer
@@ -2252,22 +2465,21 @@ class Lattice {
 
             prev_move_type.clear();
             prev_move_type_ticks.clear();
-            /*
-            Matrix<int> only_vacancies = vacancies.nonzero();
-        
-            std::cout << "rank: " << rank << " vacancies_pos.rows: " << vacancies_pos.rows() << "\n";
-            std::cout << "rank: " << rank << " num of vacancies: " << only_vacancies.rows() << "\n";
-            std::cout << "rank: " << rank << " vacancies_pos: \n";
-            vacancies_pos.print();
-            std::cout << "rank: " << rank << " only_vacancies: \n";
-            only_vacancies.print();
-            */
 
         }
 
-        /*
-        method to receive move from other process - utilizes nonblocking receive to
-        update lattice prior to finding new moves
+
+        /**
+        * @brief Receives a move from another process using non-blocking communication.
+        *
+        * This function updates the lattice with the received move details while checking for 
+        * potential conflicts with existing vacancies. It extracts the new location and previous 
+        * location data from the input buffer, performs necessary checks, and stores relevant 
+        * information about the move.
+        *
+        * @param new_loc_buffer A vector containing the new location and metadata of the move.
+        * @param move_ticks The number of ticks associated with the received move.
+        * @return true if an interboundary conflict occurred, false otherwise.
         */
         bool recieve_move_parallel(const std::vector<int>& new_loc_buffer, int move_ticks) {
             //
@@ -2304,7 +2516,7 @@ class Lattice {
             
             
             if (vacancies((size_t)i, (size_t)j, (size_t)k, (size_t)l) == 1) {
-                //std::cout << "rank: " << rank << " ERROR: Vacancy overwriting already existing vacancy in interprocessor communication \n";
+                std::cout << "rank: " << rank << " ERROR: Vacancy overwriting already existing vacancy in interprocessor communication \n";
                 interboundary_conflict = true;
                 comm_boundary_conflict();
             }
@@ -2314,7 +2526,10 @@ class Lattice {
                 size_t rows = vacancies_pos.rows();
                 size_t cols = vacancies_pos.cols();
                 vacancies_pos.reshape(rows+1, cols, rank);
-                                
+                //std::cout << "rank: " << rank << " parallel recv site status: " << vacancies(i,j,k,l) << " \n"; 
+                //std::cout << "rank: " << rank << " parallel recv vacancies_pos(vacs_idx): [ " << vacancies_pos(rows,0) << " " << vacancies_pos(rows,1) << " " << vacancies_pos(rows,2) << " " << vacancies_pos(rows,3) << " ]\n"; 
+
+                
                 store_parallel_info(new_loc_buffer, move_ticks, rows);
                 
                 vacancies((size_t)i, (size_t)j, (size_t)k, (size_t)l) = 1;
@@ -2334,9 +2549,15 @@ class Lattice {
             return interboundary_conflict;
         }
 
-        /*
-        recieve communication (ghost site update, lattice site update, conflict) using
-        blocking communication
+
+
+        /**
+        * @brief Receive communication (ghost site update, lattice site update, conflict) using blocking communication.
+        *
+        * This function handles the reception of messages from other processes,
+        * including updates to ghost sites, lattice sites, and conflicts.
+        *
+        * @param move_ticks The current iteration for moves in the simulation.
         */
         void receive_parallel_comm_helper(int move_ticks) {
 
@@ -2408,16 +2629,11 @@ class Lattice {
             }
 
             MPI_Request request1;
-            for (int new_proc=0; new_proc<num_procs; new_proc++) { 
-                if (new_proc != rank) {
-                    MPI_Isend( NULL, 0, MPI_CHAR, new_proc, conflict_done_flag, MPI_COMM_WORLD, &request1); 
-                    MPI_Wait(&request1, MPI_STATUS_IGNORE);
-                }
-            }
+            for (int new_proc=0; new_proc<num_procs; new_proc++) { if (new_proc != rank) MPI_Isend( NULL, 0, MPI_CHAR, new_proc, conflict_done_flag, MPI_COMM_WORLD, &request1); }
 
             remove_old_par_moves(move_ticks);
 
-            //MPI_Barrier(MPI_COMM_WORLD);
+            MPI_Barrier(MPI_COMM_WORLD);
             
             while ((!stop_conflict)) {
                 stop_conflict = 1;
@@ -2437,20 +2653,31 @@ class Lattice {
                 }
             }
 
-            //MPI_Barrier(MPI_COMM_WORLD);
+            MPI_Barrier(MPI_COMM_WORLD);
             
             if ((need_reverse) || (interboundary_conflict)) {
                 deincrement_time();
                 reverse_moves_wrapper();
             }
         }
+        
 
-        /*
-        wrapper method containing initialization of all varables in system, timer, and 
-        calls to update state of system and list of moves
+        /**
+        * @brief Wrapper method containing initialization of all variables in the system, timer, and calls to update state.
+        *
+        * This method initializes necessary variables and manages the iterative
+        * process of the Kinetic Monte Carlo simulation. It handles time limits,
+        * updates the state of the system, and logs results to specified data 
+        * structures and folders.
+        *
+        * @param time_lim The time limit for the simulation.
+        * @param start The starting time point of the simulation.
+        * @param folder The folder path for output data.
+        * @param iteration The current iteration of the simulation.
+        * @param rates_i An index representing the rate being processed.
+        * @return A structure containing the results of the KMC simulation.
         */
-        // TO CONNECT TO VNC: ssh -f -N -L xxxx:localhost:yyyy cfrech@ls6.tacc.utexas.edu
-        lattice_return_struct new_kmc_iterator(double time_lim, std::chrono::system_clock::time_point start, std::string folder, int iteration, int rates_i, int last_tick, double last_time) {
+        lattice_return_struct new_kmc_iterator(double time_lim, std::chrono::system_clock::time_point start, std::string folder, int iteration, int rates_i) {
             fprintf(stdout, "%s", "beginning kmc iterations \n\n"); 
 
             // INITIALIZING VARIABLES PRIOR TO BEGINNING FIRST KMC STEP //
@@ -2466,12 +2693,8 @@ class Lattice {
             std::vector<double> all_times; // vector containing trajectory of time elapsed by each type of move
             parallel_get_actions(); // updating list of moves in system
             
-            t = last_time;             
-            int move_ticks = last_tick;
+            int move_ticks = 0;
             double old_time;
-            bool restart = false;
-
-            if (move_ticks != 0) { restart = true; }
 
             int prev_num_vacs = 0;
             int curr_num_vacs = 0;
@@ -2479,13 +2702,94 @@ class Lattice {
 
             std::ostringstream ss;
 
+
             Matrix<int> only_vacancies = vacancies.nonzero();
             std::cout << "rank: " << rank << " move_ticks: " << move_ticks << " vacancies_pos.rows(): " << vacancies_pos.rows() << " vs only_vacancies.rows(): " << only_vacancies.rows() << "\n";
             
             MPI_Barrier(MPI_COMM_WORLD);
+            if (rank == 0) {
+                std::cout << "rank: " << rank << " vacancies_pos.rows: " << vacancies_pos.rows() << "\n";
+                std::cout << "rank: " << rank << " num of vacancies: " << only_vacancies.rows() << "\n";
+                std::cout << "rank: " << rank << " vacancies_pos: \n";
+                vacancies_pos.print();
+                std::cout << "rank: " << rank << " only_vacancies: \n";
+                only_vacancies.print();
+                if (vacancies_pos.rows() != only_vacancies.rows()) { 
+                    std::cout << "rank: " << rank << " vacancies_pos: \n";
+                    vacancies_pos.print();
+                    std::cout << "rank: " << rank << " only_vacancies: \n";
+                    only_vacancies.print();
+                    //std::cout << "rank: " << rank << " vacancies_pos.rows(): " << vacancies_pos.rows() << " mismatch with only_vacancies.rows(): " << only_vacancies.rows() << "\n"; 
+                    //Matrix<int> difference = comparison(vacancies_pos, only_vacancies);
+                    //std::cout << "rank: " << rank << " difference: \n";
+                    //difference.print();
+                    exit(0);
+                }
+            }
+            MPI_Barrier(MPI_COMM_WORLD);
+            if (rank == 1) {
+                std::cout << "rank: " << rank << " vacancies_pos.rows: " << vacancies_pos.rows() << "\n";
+                std::cout << "rank: " << rank << " num of vacancies: " << only_vacancies.rows() << "\n";
+                std::cout << "rank: " << rank << " vacancies_pos: \n";
+                vacancies_pos.print();
+                std::cout << "rank: " << rank << " only_vacancies: \n";
+                only_vacancies.print();
+                if (vacancies_pos.rows() != only_vacancies.rows()) { 
+                    std::cout << "rank: " << rank << " vacancies_pos: \n";
+                    vacancies_pos.print();
+                    std::cout << "rank: " << rank << " only_vacancies: \n";
+                    only_vacancies.print();
+                    //std::cout << "rank: " << rank << " vacancies_pos.rows(): " << vacancies_pos.rows() << " mismatch with only_vacancies.rows(): " << only_vacancies.rows() << "\n"; 
+                    //Matrix<int> difference = comparison(vacancies_pos, only_vacancies);
+                    //std::cout << "rank: " << rank << " difference: \n";
+                    //difference.print();
+                    //exit(0);
+                }
+            }
+            MPI_Barrier(MPI_COMM_WORLD);
+            if (rank == 2) {
+                std::cout << "rank: " << rank << " vacancies_pos.rows: " << vacancies_pos.rows() << "\n";
+                std::cout << "rank: " << rank << " num of vacancies: " << only_vacancies.rows() << "\n";
+                std::cout << "rank: " << rank << " vacancies_pos: \n";
+                vacancies_pos.print();
+                std::cout << "rank: " << rank << " only_vacancies: \n";
+                only_vacancies.print();
+                if (vacancies_pos.rows() != only_vacancies.rows()) { 
+                    std::cout << "rank: " << rank << " vacancies_pos: \n";
+                    vacancies_pos.print();
+                    std::cout << "rank: " << rank << " only_vacancies: \n";
+                    only_vacancies.print();
+                    //std::cout << "rank: " << rank << " vacancies_pos.rows(): " << vacancies_pos.rows() << " mismatch with only_vacancies.rows(): " << only_vacancies.rows() << "\n"; 
+                    //Matrix<int> difference = comparison(vacancies_pos, only_vacancies);
+                    //std::cout << "rank: " << rank << " difference: \n";
+                    //difference.print();
+                    //exit(0);
+                }
+            }
+            MPI_Barrier(MPI_COMM_WORLD);
+            if (rank == 3) {
+                std::cout << "rank: " << rank << " vacancies_pos.rows: " << vacancies_pos.rows() << "\n";
+                std::cout << "rank: " << rank << " num of vacancies: " << only_vacancies.rows() << "\n";
+                std::cout << "rank: " << rank << " vacancies_pos: \n";
+                vacancies_pos.print();
+                std::cout << "rank: " << rank << " only_vacancies: \n";
+                only_vacancies.print();
+                if (vacancies_pos.rows() != only_vacancies.rows()) { 
+                    std::cout << "rank: " << rank << " vacancies_pos: \n";
+                    vacancies_pos.print();
+                    std::cout << "rank: " << rank << " only_vacancies: \n";
+                    only_vacancies.print();
+                    //std::cout << "rank: " << rank << " vacancies_pos.rows(): " << vacancies_pos.rows() << " mismatch with only_vacancies.rows(): " << only_vacancies.rows() << "\n"; 
+                    //Matrix<int> difference = comparison(vacancies_pos, only_vacancies);
+                    //std::cout << "rank: " << rank << " difference: \n";
+                    //difference.print();
+                    //exit(0);
+                }
+            }
             
-            while (elapsed_seconds.count()  < 900) {
-                //if (rank == 0) std::cout << "rank: " << rank << " move_ticks: " << move_ticks << "\n";
+            MPI_Barrier(MPI_COMM_WORLD);
+
+            while (t < 250000) {
 
                 if ( rate_cumsum.size() != 0) {
                     end = std::chrono::system_clock::now(); 
@@ -2497,97 +2801,68 @@ class Lattice {
                     store_time_incr(timestep);
                 }
                 else {
+                    std::cout << "rank: " << rank << " only communicate_rates " << "\n";
                     communicate_rates(); 
                     timestep = new_random_times();
                 }
+                MPI_Barrier(MPI_COMM_WORLD);
 
                 MPI_Request request1;
-                for (int new_proc=0; new_proc<num_procs; new_proc++) { 
-                    if (new_proc != rank) {
-                        MPI_Isend( NULL, 0, MPI_CHAR, new_proc, par_done_tag, MPI_COMM_WORLD, &request1); 
-                        MPI_Wait(&request1, MPI_STATUS_IGNORE);
-                    }
-                }
+                for (int new_proc=0; new_proc<num_procs; new_proc++) { if (new_proc != rank) MPI_Isend( NULL, 0, MPI_CHAR, new_proc, par_done_tag, MPI_COMM_WORLD, &request1); }
                 MPI_Request request2;
-                for (int new_proc=0; new_proc<num_procs; new_proc++) { 
-                    if (new_proc != rank) {
-                        MPI_Isend( NULL, 0, MPI_CHAR, new_proc, ghost_done_tag, MPI_COMM_WORLD, &request2); 
-                        MPI_Wait(&request2, MPI_STATUS_IGNORE);
-                    }
-                }
+                for (int new_proc=0; new_proc<num_procs; new_proc++) { if (new_proc != rank) MPI_Isend( NULL, 0, MPI_CHAR, new_proc, ghost_done_tag, MPI_COMM_WORLD, &request2); }
+                
                 
                 MPI_Barrier(MPI_COMM_WORLD);
                 
+                
                 receive_parallel_comm_helper(move_ticks);
 
-                if (move_ticks % 100000 == 0) {
-                    if (rank == 0) { std::cout << "rank: " << rank << " move_ticks: " << move_ticks << "\n"; }
-                    //Matrix<int> only_vacancies = vacancies.nonzero(); // configuration of vacancies at current timestep
-                    //write_output_parallel(only_vacancies, total_dims, folder, proc_dims[0], proc_dims[1], num_procs, rank, iteration, rates_i, move_ticks, t, get_rank);
+                MPI_Barrier(MPI_COMM_WORLD);
+               
+                if (move_ticks % 10000 == 0) {
+                    Matrix<int> only_vacancies = vacancies.nonzero(); // configuration of vacancies at current timestep
+                    write_output_parallel(only_vacancies, total_dims, folder, proc_dims[0], proc_dims[1], num_procs, rank, iteration, rates_i, move_ticks, t, get_rank);
                 }
+            
+                
+                MPI_Barrier(MPI_COMM_WORLD);
+
 
                 parallel_get_actions();
                 
                 fflush(stdout);
+                MPI_Barrier(MPI_COMM_WORLD);
                 
-                
-                /*
                 curr_num_vacs = sum_vacs_allprocs(only_vacancies, total_dims, proc_dims[0], proc_dims[1], num_procs, rank);
+
+                only_vacancies = vacancies.nonzero();
                 
                 if ((rank == 0) && (prev_num_vacs != curr_num_vacs)) {
                     std::cout << "Change in total number of vacancies in simulation -- curr_num_vacs: " << curr_num_vacs << " prev_num_vacs: " << prev_num_vacs << " move_ticks: " << move_ticks << "\n"; 
-                    if ((move_ticks != 0) && (restart == false)) {
+                    if (move_ticks != 0) {
                         MPI_Barrier(MPI_COMM_WORLD);
                         exit(0);
                     }
                 }
-                */
-                
-                if (restart == true) restart = false;
                 prev_num_vacs = curr_num_vacs;
 
                 elapsed_seconds = end-start;
                 t += timestep;
-                
                 move_ticks ++;
-                old_time = t;
+                old_time = t;                
 
                 // terminating simulation after real-time limit reached
-                if (elapsed_seconds.count() >= 900) { 
-                    std::cout << "rank: " << rank << " elapsed_seconds: " << elapsed_seconds.count() << "\n";
-                    std::cout << "rank: " << rank << " exiting\n";
-                    std::cout << "rank: " << rank << " t: " << t << "\n";
-                    std::cout << "rank: " << rank << " move_ticks: " << move_ticks << "\n";
-                    exit(0);
-
-                    if (rank == 0) {
-                        std::cout << "rank: " << rank << " move_counts_sum: \n";
-                        print_1Dvector(move_counts);
-                        std::cout << "rank: " << rank << " time_count: \n";
-                        print_1Dvector(time_count);
-                    }
+                if (elapsed_seconds.count() >= 172800) {
+                    t = time_lim + 1;
+                    std::cout << "elapsed_seconds: " << elapsed_seconds.count() << "\n";
+                    std::cout << "exiting\n";
+                    std::cout << "t: " << t << "\n";
+                    print_1Dvector(move_counts);
+                    break;
                 }
+                MPI_Barrier(MPI_COMM_WORLD);
             }
-            exit(0);
-            //if (move_ticks >= 75000) {
-            std::cout << "rank: " << rank << " elapsed_seconds: " << elapsed_seconds.count() << "\n";
-            std::cout << "rank: " << rank << " exiting\n";
-            std::cout << "rank: " << rank << " t: " << t << "\n";
-            std::cout << "rank: " << rank << " move_ticks: " << move_ticks << "\n";
-
-            std::vector<int> move_counts_sum = sum_vectors_allprocs(move_counts, num_procs, rank);
-            std::vector<double> time_count_sum = sum_vectors_allprocs(time_count, num_procs, rank);
-
-            if (rank == 0) {
-                std::cout << "rank: " << rank << " move_counts_sum: \n";
-                print_1Dvector(move_counts_sum);
-                std::cout << "rank: " << rank << " time_count: \n";
-                print_1Dvector(time_count);
-            }
-
-            t = time_lim + 1;
-            fflush(stdout);
-            MPI_Barrier(MPI_COMM_WORLD);
 
             all_times.push_back(t);
 
@@ -2595,18 +2870,26 @@ class Lattice {
 
             lattice_return_struct output_vals(move_counts, time_count, all_times);
 
-            //MPI_Barrier(MPI_COMM_WORLD);
+            MPI_Barrier(MPI_COMM_WORLD);
             return output_vals;
         }
 };
 
 /*---------------------------------------------------------------------------*/
 
-/*
-creating entries in rate catalog corresponding to all configurations of
-"len" sites with m vacancies, corresponding to all the binary strings of 
-length "len" with m zeros.
-*/
+/**
+ * @brief Creates entries in the rate catalog corresponding to all configurations of "len" sites with m vacancies.
+ *
+ * This function generates all binary strings of length "len" with exactly m zeros,
+ * representing configurations of vacancies. The output is a vector of integers,
+ * where each integer corresponds to a specific binary configuration.
+ *
+ * @param len The total length of the binary string (total number of sites).
+ * @param m The number of zeros in the binary string (vacancies).
+ * @param size The number of atomic types in the system used to create the encoding.
+ * @return A vector of integers representing all binary configurations of length "len"
+ *         with m zeros.
+ */
 std::vector<int> bin_m_zeros(int len, int m, int size) {
     int num_combos = NCR(len, m);
     std::vector<int> vec(num_combos);
@@ -2632,9 +2915,17 @@ std::vector<int> bin_m_zeros(int len, int m, int size) {
     return vec;
 }
 
-/*
-converting a inputted base-m string (type std::string) to an int
-*/
+/**
+ * @brief Converts an inputted base-m string to an integer.
+ *
+ * This function takes a string representation of a number in base-m,
+ * where each character represents a digit in that base, and converts it
+ * to an integer.
+ *
+ * @param config The input string representing the number in base-m.
+ * @param size The base (m) of the input string.
+ * @return The integer value of the input string in base-m.
+ */
 int base_m_to_int(std::string config, int size) {
 
     int result = 0;
@@ -2650,10 +2941,24 @@ int base_m_to_int(std::string config, int size) {
     return result;
 }
 
-/*
-reads in file corresponding to rate catalog and creates catalog of rates with allowed
-moves, catalog of allowed configurations for moves, and rate catalog corresponding to regions
-*/
+/**
+ * @brief Reads a rate catalog file and creates a catalog of rates, allowed moves,
+ *        and configurations corresponding to regions.
+ *
+ * This function processes an input file that contains information about atom types,
+ * configurations, and associated energies. It constructs and returns a structured
+ * catalog (Matrix) that includes allowed configurations and rates for migration.
+ *
+ * @param catalogfile The path to the rate catalog file.
+ * @param atype_list A list of atom types specified for use in the catalog.
+ * @return A `ratecatalog_struct` containing the catalogs of configurations,
+ *         energies, and regions.
+ *
+ * @note The expected file format is:
+ *       - Line 0: Atom types (e.g., "1:sodium, 2:LLZO, etc...")
+ *       - Line 1: Number of configurations (e.g., "#configs: m")
+ *       - Subsequent lines: Configurations and corresponding energies.
+ */
 ratecatalog_struct updated_create_ratecatalog(std::string catalogfile, std::vector<int> atype_list) {
     
     //FILE FORMAT:
@@ -2909,9 +3214,29 @@ ratecatalog_struct updated_create_ratecatalog(std::string catalogfile, std::vect
     return output_vals;
 }
 
-/*
-creating region objects corresponding to input file
-*/
+/**
+ * @brief Creates a region object based on the provided information.
+ *
+ * This function constructs a `Region` object using data extracted from the 
+ * input vector. It supports different region types, such as "GB" (Grain Boundary)
+ * and "BLOCK" (rectangular prism). The parameters for the region are parsed
+ * from the input string vector.
+ *
+ * @param info A vector of strings containing region information:
+ *             - info[0]: Region ID (e.g., "id:1")
+ *             - info[1]: Region type (e.g., "GB" or "BLOCK")
+ *             - Subsequent elements contain parameters relevant to the region type.
+ * @return A pointer to the newly created `Region` object.
+ *
+ * @note The expected format for `info` elements varies based on the region type.
+ *       - For "GB": 
+ *         - info[2-4]: First set of parameters (3 integers)
+ *         - info[5-7]: Second set of parameters (3 integers)
+ *       - For "BLOCK":
+ *         - info[2-3]: X-dimension parameters (2 integers)
+ *         - info[4-5]: Y-dimension parameters (2 integers)
+ *         - info[6-7]: Z-dimension parameters (2 integers)
+ */
 Region* add_region(std::vector<std::string> info) {
     std::cout << "adding region \n";
     int id = std::stoi(tokenizer(info[0], ":")[0]); // region id number
@@ -2945,9 +3270,28 @@ Region* add_region(std::vector<std::string> info) {
     return new_region;
 }
 
-/*
-populating region_sites FourDArr with values corresponding to custom regions input file
-*/
+/**
+ * @brief Populates the FourDArr `sites` with values corresponding to custom regions from an input file.
+ *
+ * This function reads a custom regions file and updates the specified 
+ * `FourDArr` structure with region IDs based on the coordinates specified 
+ * in the file. The function processes each region from the provided list 
+ * and updates the 4D array for each lattice position.
+ *
+ * @param sites A pointer to a `FourDArr` object to be populated with region IDs.
+ * @param regions A vector of pointers to `Region` objects representing the regions to be processed.
+ * @param custom_reg_idx An index specifying which custom region to draw; currently not used in implementation.
+ * @param dim A vector of integers representing the dimensions of the simulation cell.
+ * @param infile_name The name of the input file containing custom region data.
+ *
+ * @note The input file is expected to contain lines of data formatted for 
+ *       each region's lattice position, which are parsed and used to update 
+ *       the `FourDArr`. The function checks that the coordinates do not exceed 
+ *       the bounds of the simulation cell dimensions.
+ *
+ * @warning If coordinates exceed simulation cell bounds, no updates will be made 
+ *          for those entries, but no exception will be thrown. 
+ */
 void custom_draw_regions(FourDArr* sites, std::vector<Region*> regions, int custom_reg_idx, std::vector<int> dim, std::string infile_name) {
     std::cout << "drawing regions \n";
     Region* region;
@@ -2991,9 +3335,28 @@ void custom_draw_regions(FourDArr* sites, std::vector<Region*> regions, int cust
     }
 }
 
-/*
-populating region_sites FourDArr with values corresponding to input file
-*/
+/**
+ * @brief Populates the `FourDArr` `sites` with values corresponding to defined regions.
+ *
+ * This function iterates through a vector of `Region` objects and updates the 
+ * provided `FourDArr` structure based on the type of each region (either 
+ * "GB" for grain boundaries or "BLOCK" for rectangular prisms). For grain 
+ * boundary regions, the function computes the appropriate coordinates based on 
+ * the slopes and shifts defined for the region. For block regions, the function 
+ * uses the lower and upper bounds to determine the coordinates to populate.
+ *
+ * @param sites A pointer to a `FourDArr` object where the region IDs will be assigned.
+ * @param regions A vector of pointers to `Region` objects that define the regions to draw.
+ * @param dim A vector of integers representing the dimensions of the simulation cell.
+ *
+ * @note The function computes the coordinates for each region based on its 
+ *       specific parameters and assigns the corresponding region ID to those 
+ *       coordinates in the `FourDArr`.
+ *
+ * @warning Ensure that the coordinates calculated do not exceed the bounds of 
+ *          the `FourDArr`, as this implementation does not include checks for 
+ *          out-of-bounds access.
+ */
 void draw_regions(FourDArr* sites, std::vector<Region*> regions, std::vector<int> dim ) {
     Region* region;
     std::vector<int> lo(3);
@@ -3057,10 +3420,33 @@ void draw_regions(FourDArr* sites, std::vector<Region*> regions, std::vector<int
     }
 }
 
-/*
-wrapper function reading input file and creating region objects / populating region_sites FourDArr
-*/
-// std::tuple< int, std::vector<Region*>, FourDArr* > 
+/**
+ * @brief Initializes region objects from input data and populates a `FourDArr` with region sites.
+ *
+ * This wrapper function reads lines from an input file to create `Region` objects
+ * and initializes their corresponding site entries in a `FourDArr`. The function
+ * handles both regions of predefined shape (block or slab) and custom regions based 
+ * defined by a coordinate file depending whether a custom region input file is provided. 
+ * If the number of regions specified exceeds those defined in the input, additional 
+ * block regions are generated with default parameters.
+ *
+ * @param lines A vector of strings containing lines read from the input file.
+ * @param dims A vector of integers representing the dimensions of the simulation space.
+ * @param num_regions The total number of regions to initialize.
+ * @param region_infile The name of the input file for custom regions (can be empty).
+ *
+ * @return add_reg_struct A structure containing information about the initialized regions,
+ *         the number of lines read, and the populated `FourDArr` which indicates which sites
+ *         are specially-defined regions.
+ *
+ * @note The function initializes Region objects based upon the 'num_regions' information 
+ *       and region definitions in the regions input file. It also ensures that any 
+ *       additional regions needed are created with default block parameters.
+ *
+ * @warning Make sure that the input lines are properly formatted to avoid runtime errors 
+ *          when parsing region definitions. Also, ensure that the `region_infile` 
+ *          provided is valid if custom regions are required.
+ */
 add_reg_struct
 init_regions(std::vector<std::string> lines, 
 std::vector<int> dims, int num_regions, std::string region_infile) {
@@ -3133,14 +3519,28 @@ std::vector<int> dims, int num_regions, std::string region_infile) {
     return returnval;
 }
 
-/*
-wrapper function for:
-- reading input files
-- populating vacancy, bc_site, vertex_site, and region_site FourDArr data structures
-corresponding to initial configuration of simulation
--  creating region objects 
-- creating rate catalog for bulk & pre-defined regions
-*/
+/**
+ * @brief Wrapper function for populating the lattice.
+ *
+ * This function reads input files to populate FourDArr and FourDBoolArr 
+ * data structures for vacancy and atomic positions corresponding to the 
+ * initial configuration of the simulation. It also defines regions within
+ * the simulation cell and builds a rate catalog for bulk and pre-defined 
+ * regions.
+ *
+ * @param infile_name The name of the input file to read.
+ * @param catalogfile_name The name of the rate catalog file.
+ * @param region_infile The name of the region input file.
+ * @param vertex_rate The rate for vertex interactions.
+ * @param edge_rate The rate for edge interactions.
+ * @param reg_rates A vector of vectors containing region rates.
+ * @param total_dims A vector containing the total dimensions of the lattice.
+ * @param chunk_bounds A vector of vectors defining the bounds of the chunks.
+ * @param rank The rank of the current process in a parallel computation.
+ * @param procs A vector containing the number of processes in each dimension.
+ *
+ * @return A pointer to the populated Lattice object.
+ */
 Lattice* populate_lattice(std::string infile_name, std::string catalogfile_name, std::string region_infile, double vertex_rate, double edge_rate, 
 std::vector<std::vector<double>> reg_rates, std::vector<int> total_dims, std::vector<std::vector<int>> chunk_bounds, int rank, std::vector<int> procs) {
 
@@ -3235,7 +3635,7 @@ std::vector<std::vector<double>> reg_rates, std::vector<int> total_dims, std::ve
 
     if (lines[read_idx].find(substring) != std::string::npos) {
         read_idx ++;
-        add_reg_struct regions_tuple = init_regions(slice_1Dvec_str(lines, read_idx, (int)lines.size()), dims_int, num_regions, region_infile); 
+        add_reg_struct regions_tuple = init_regions(slice_1Dvec_str_inp(lines, read_idx, (int)lines.size()), dims_int, num_regions, region_infile); 
         std::cout << "region!\n";
         incriment = regions_tuple.get_idx(); temp_regions = regions_tuple.get_regions(); temp_region_sites = regions_tuple.get_region_sites();
     }
@@ -3661,25 +4061,25 @@ std::vector<std::vector<double>> reg_rates, std::vector<int> total_dims, std::ve
     }
 
     new_configs.push_back(lateral_configs);
-    new_energies.push_back(lateral_E); //invalid read
+    new_energies.push_back(lateral_E);
 
     std::vector<int> unsorted_idxs_8bit = arange(0,(int)new_configs[0].size(),1);
     auto comparator_8bit = [new_configs](int idx1, int idx2) {
                 return new_configs[0][idx1] < new_configs[0][idx2];
         };
 
-    std::sort(unsorted_idxs_8bit.begin(), unsorted_idxs_8bit.end(), comparator_8bit); 
+    std::sort(unsorted_idxs_8bit.begin(), unsorted_idxs_8bit.end(), comparator_8bit);
     new_configs[0] = reorder_inp(new_configs[0], unsorted_idxs_8bit);
     new_energies[0] = reorder_inp(new_energies[0], unsorted_idxs_8bit);
 
-    std::vector<int> unsorted_idxs_14bit = arange(0,(int)new_configs[1].size(),1); //invalid read
-    auto comparator_14bit = [new_configs](int idx1, int idx2) { //invalid read
+    std::vector<int> unsorted_idxs_14bit = arange(0,(int)new_configs[1].size(),1);
+    auto comparator_14bit = [new_configs](int idx1, int idx2) {
                 return new_configs[1][idx1] < new_configs[1][idx2];
         };
 
     std::sort(unsorted_idxs_14bit.begin(), unsorted_idxs_14bit.end(), comparator_14bit);
     new_configs[1] = reorder_inp(new_configs[1], unsorted_idxs_14bit);
-    new_energies[1] = reorder_inp(new_energies[1], unsorted_idxs_14bit); //invalid read
+    new_energies[1] = reorder_inp(new_energies[1], unsorted_idxs_14bit);
 
     for (int i=0; i<(int)new_configs[0].size(); i++) {
        configs_111[0][i] = new_configs[0][i]; 
@@ -3704,25 +4104,7 @@ std::vector<std::vector<double>> reg_rates, std::vector<int> total_dims, std::ve
     int num_x_neigh = dims_int[0] + 1;
     int num_y_neigh = dims_int[1] + 1;
 
-    Lattice* new_lattice = new Lattice(dims_int[0], dims_int[1], dims_int[2], vacancies_count, num_regions, nprocs, num_x_neigh, num_y_neigh, total_dims[0], total_dims[1], total_dims[2], rank);
-
-    std::vector < std::vector <int> > temp_diag = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0}, {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}};
-    std::vector < std::vector <int> > temp_edge = {{0,0,1}, {-1,0,0}, {0,-1,0}, {0,1,0}, {1,0,0}, {0,0,-1}};
-    
-    
-    for (int i=0; i < temp_diag.size(); i++) {
-        for (int j=0; j < temp_diag[0].size(); j++) {
-            new_lattice->diag_directions(i,j) =temp_diag[i][j];
-        }
-    } 
-
-    for (int i=0; i < temp_edge.size(); i++) {
-        for (int j=0; j < temp_edge[0].size(); j++) {
-            new_lattice->edge_directions(i,j) = temp_edge[i][j];
-        }
-    } 
-
-
+    Lattice* new_lattice = new Lattice(dims_int[0], dims_int[1], dims_int[2], vacancies_count, num_regions, nprocs, num_x_neigh, num_y_neigh, total_dims[0], total_dims[1], total_dims[2]);
 
     for (size_t i=0; i<2; i++) {
         for (size_t j=0; j<(size_t)dims_int[0]; j++) {
@@ -3735,7 +4117,7 @@ std::vector<std::vector<double>> reg_rates, std::vector<int> total_dims, std::ve
                         new_lattice->bc_sites(0,j,k,l) = (*temp_bc_sites)(0,j,k,l);
                     }
 
-                    new_lattice->region_sites(i,j,k,l) = (*temp_region_sites)(i,j,k,l); // invalid read
+                    new_lattice->region_sites(i,j,k,l) = (*temp_region_sites)(i,j,k,l);
                     new_lattice->vacancies(i,j,k,l) = (*temp_vacancies)(i,j,k,l);
                 }
             }
@@ -3769,6 +4151,7 @@ std::vector<std::vector<double>> reg_rates, std::vector<int> total_dims, std::ve
         }
     }
     
+    new_lattice->rank = rank;
     new_lattice->configs_111 = configs_111;
     new_lattice->configs_100 = configs_100;
 
@@ -3829,7 +4212,7 @@ std::vector<std::vector<double>> reg_rates, std::vector<int> total_dims, std::ve
     nonzero_vacs.print();
     
     new_lattice->proc_neighbors = temp_proc_neighbors;
-    /*
+   
     std::cout << "proc_neighbors: \n";
     new_lattice->proc_neighbors.print();
 
@@ -3892,9 +4275,9 @@ std::vector<std::vector<double>> reg_rates, std::vector<int> total_dims, std::ve
         std::cout << "rank: " << rank << " new_lattice->proc_pos_y_neighbors: \n";
         new_lattice->proc_pos_y_neighbors.print();
     }
-    */
-    new_lattice->void_threshold = 3;
-    new_lattice->void_barrier = 5.06e7;
+    
+    new_lattice->void_threshold = 2;
+    new_lattice->void_barrier = 2761667.608;
 
 
     delete temp_111_catalog;
@@ -3916,3 +4299,4 @@ std::vector<std::vector<double>> reg_rates, std::vector<int> total_dims, std::ve
     return new_lattice;
 }
 /*---------------------------------------------------------------------------*/
+

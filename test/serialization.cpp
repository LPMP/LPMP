#include "test.h"
#include "serialization.hxx"

using namespace LPMP;

template<typename T>
const INDEX archive_size(T&& e)
{
   allocate_archive ar;
   ar(e);
   return ar.size();
}

int main()
{

   // INDEX valued types
   INDEX i=10;
   std::array<INDEX,2> a{20,30};
   std::vector<INDEX> v{30,40,50};
   INDEX p [4] = {60,70,80,90};

   // REAL valued types
   REAL i_r=10;
   std::array<REAL,2> a_r{20,30};
   std::vector<REAL> v_r{30,40,50};
   REAL p_r [4] = {60,70,80,90};
   LPMP::vector<REAL> r_r(3); v_r[0] = 10.0; v_r[1] = 20.0; v_r[2] = 30.0;
   matrix<REAL> m_r(5,3);
   m_r(0,0) = 1.0; m_r(1,0) = 2.0; m_r(2,0) = 3.0; m_r(3,0) = 4.0; m_r(4,0) = 5.0;
   m_r(0,1) = 11.0; m_r(1,1) = 12.0; m_r(2,1) = 13.0; m_r(3,1) = 14.0; m_r(4,1) = 15.0;
   m_r(0,2) = 111.0; m_r(1,2) = 112.0; m_r(2,2) = 113.0; m_r(3,2) = 114.0; m_r(4,2) = 115.0;

   allocate_archive a_ar;
   a_ar(i);
   a_ar(a);
   a_ar(v);
   a_ar( binary_data<INDEX>(p,4) );
   a_ar(r_r);
   a_ar(m_r);

    { // allocate archive
      test(archive_size(i) == sizeof(INDEX));
      test(archive_size(a) == 2*sizeof(INDEX));
      test(archive_size(v) == 3*sizeof(INDEX));
      test(archive_size(binary_data<INDEX>(p, 4)) == 4*sizeof(INDEX));

      test(archive_size(i_r) == sizeof(REAL));
      test(archive_size(a_r) == 2*sizeof(REAL));
      test(archive_size(v_r) == 3*sizeof(REAL));
      test(archive_size(binary_data<REAL>(p_r, 4)) == 4*sizeof(REAL));
      test(archive_size(r_r) == 3*sizeof(REAL));
      test(archive_size(m_r) == 15*sizeof(REAL));
   }


   { // individual saving
      serialization_archive ar(a_ar);

      save_archive s_ar(ar);
      s_ar(i);
      s_ar(a);
      s_ar(v);
      s_ar( binary_data<INDEX>(p,4) );
      s_ar(r_r);
      s_ar(m_r);

      load_archive l_ar(ar);

      decltype(i) i_test;
      l_ar(i_test);
      test(i_test == i);

      decltype(a) a_test;
      l_ar(a_test);
      test(a_test == a);

      decltype(v) v_test(v.size());
      l_ar(v_test);
      test(v_test == v);

      decltype(p) p_test = {0,0,0,0};
      l_ar( binary_data<INDEX>(p_test, 4) );
      test(p_test[0] == p[0]);
      test(p_test[1] == p[1]);
      test(p_test[2] == p[2]);
      test(p_test[3] == p[3]);

      decltype(r_r) r_test(3);
      l_ar( r_test );
      test(r_test == r_r);

      decltype(m_r) m_test(m_r.dim1(), m_r.dim2());
      l_ar( m_test );
      test(m_test == m_r);
   }

   { // collective saving
      serialization_archive ar(a_ar);

      save_archive s_ar(ar);
      s_ar(i, a, v, binary_data<INDEX>(p,4), r_r, m_r);

      load_archive l_ar(ar);

      decltype(i) i_test;
      decltype(a) a_test;
      decltype(v) v_test(v.size());
      decltype(p) p_test = {0,0,0,0};
      decltype(r_r) r_test(r_r.size());
      decltype(m_r) m_test(m_r.dim1(), m_r.dim2());

      l_ar(i_test, a_test, v_test, binary_data<INDEX>(p_test,4), r_test, m_test);

      test(i_test == i);
      test(a_test == a);
      test(v_test == v);
      test(p_test[0] == p[0]);
      test(p_test[1] == p[1]);
      test(p_test[2] == p[2]);
      test(p_test[3] == p[3]);
      test(r_test == r_r);
      test(m_test == m_r);
   } 
}


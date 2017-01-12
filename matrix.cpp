
#include "matrix.h"
#include "eig.h"

int lat_pack(const Lattice_Vector &x) {
  int i, site = 0;

  for (i = 0;i<D;i++) {
    site = site+Lattice_Map[i]*x.get(i);
  }

  return(site);
}

void lat_unpack(const int n, Lattice_Vector &x) {
  int current, i;

  current = n;
  for (i = D-1;i>= 0;i--) {
    x.set(i, current/Lattice_Map[i]);
    current = current-Lattice_Map[i]*x.get(i);
  }

  return;
}

int pack(const int flavor, const Lattice_Vector &x, const int color) {
  int i;
  i =(SITES*NUMGEN)*flavor+NUMGEN*lat_pack(x)+color;

  return(i);
}

void unpack(const int i, int &flavor, Lattice_Vector &x, int &color) {
  int tmp;
  flavor = i/(SITES*NUMGEN);
  tmp = i-flavor*(SITES*NUMGEN);
  lat_unpack(tmp/NUMGEN, x);
  color = tmp-(tmp/NUMGEN)*NUMGEN;
  return;
}


int flavor(const int mu, const int nu) {
  int dum;

  if (mu == 0) {
    dum = nu;}
  if (mu == 1) {
    dum = 5+nu-2;}
  if (mu == 2) {
    dum = 8+nu-3;}
  if (mu == 3) {
    dum = 10;}
  return(dum);
}


int find_index(const int i, const int j, int col[LEN], int row[LEN+1]) {
  int k;

  //cout << "i and j are " << i << "\t" << j << "\n" << flush;

  for (k = row[i];k<row[i]+num_in_row[i];k++) {
    if (col[k]== j) {
      //cout << "found col " << j << "already with index " << k << "\n" << flush;
      return(k);}
  }
  // new col
  // check that lies in allowed range
  k = row[i]+num_in_row[i];
  if (k<row[i+1]) {
    //cout << "new nonzero at " << k << "\n" << flush;
    num_in_row[i]++;
    //cout << "number in row " << i << " is " << num_in_row[i] << "\n" << flush;
    return(k);
  }
  else {
    cout << "error building matrix: num_in_row exceeds nonzeroes \n";
    exit(1);
  }
}


void build_vector(const Twist_Fermion &t, Complex v[LEN]) {
  int sites, mu, nu, index, a, num_chis;
  Lattice_Vector x;

  num_chis =(NUMLINK*(NUMLINK-1))/2;

  sites = 0;
  while(loop_over_lattice(x, sites)) {
    for (a = 0;a<NUMGEN;a++) {

      index = pack(0, x, a);
      v[index]= t.getS().get(x).get(a);

      for (mu = 0;mu<NUMLINK;mu++) {
        for (nu = mu+1;nu<NUMLINK;nu++) {
          index = pack(flavor(mu, nu), x, a);
          v[index]= t.getC().get(x, mu, nu).get(a);}}

      for (mu = 0;mu<NUMLINK;mu++) {
        index = pack((num_chis+1+mu), x, a);
        v[index]= t.getL().get(x, mu).get(a);}
    }
  }
  return;
}

Twist_Fermion extract_vector(Complex v[LEN]) {
  int sites, mu, nu, index, a, num_chis;
  Lattice_Vector x;
  Site_Field s = Site_Field();
  Link_Field l = Link_Field();
  Plaq_Field p = Plaq_Field();
  Afield atmp;
  Twist_Fermion t;

  num_chis =(NUMLINK*(NUMLINK-1))/2;
  sites = 0;
  while(loop_over_lattice(x, sites)) {

    for (a = 0;a<NUMGEN;a++) {
      index = pack(0, x, a);
      atmp.set(a, v[index]);}

    s.set(x, atmp);

    for (mu = 0;mu<NUMLINK;mu++) {
      for (nu = mu+1;nu<NUMLINK;nu++) {

        for (a = 0;a<NUMGEN;a++) {
          index = pack(flavor(mu, nu), x, a);
          atmp.set(a, v[index]);}

        p.set(x, mu, nu, atmp);
        p.set(x, nu, mu,-1.0*atmp);
      }}

    for (mu = 0;mu<NUMLINK;mu++) {
      for (a = 0;a<NUMGEN;a++) {
        index = pack(mu+num_chis+1, x, a);
        atmp.set(a, v[index]);}

      l.set(x, mu, atmp);}

  }

  t.setS(s);
  t.setL(l);
  t.setC(p);

  return(t);
}



void build_sparse_matrix(const Adjoint_Links &V, const Gauge_Field &U,
    Complex m[], int col[], int row[]) {
  int sites, mu, nu, a, b, c, ix, i, j1, j2, d, e, l, k;
  Lattice_Vector x, e_mu, e_nu, e_a, e_b, e_c;
  static int first_time = 1;
  static ofstream f_op;
  Complex W, Z;
  UPlaq_Field P;

  //if (first_time) {
  //  f_op.open("sparse_fermion_op");
  //  if (f_op.bad()) {
  //  cout << "failed to open fermion op file\n" << flush ;}
  //}

  int num_chis =(NUMLINK*(NUMLINK-1))/2;

  for (i = 0;i<LEN;i++) {
    row[i]= NONZEROES*i;
    num_in_row[i]= 0;
  }
  row[LEN]= NONZEROES*LEN;

  // set off initially
  for (i = 0;i<(LEN*NONZEROES);i++) {
    col[i]=-1;
    m[i]= Complex();
  }

  // eta Dminus psi_mu
  sites = 0;
  while(loop_over_lattice(x, sites)) {
    for (a = 0;a<NUMGEN;a++) {
      i = pack(0, x, a);

      for (mu = 0;mu<NUMLINK;mu++) {
        e_mu = Lattice_Vector(mu);
        for (b = 0;b<NUMGEN;b++) {
          j1 = pack(mu+num_chis+1, x, b);
          j2 = pack(mu+num_chis+1, x-e_mu, b);
          ix = find_index(i, j1, col, row);
          m[ix]= m[ix]+
            0.5*conjug(V.get(x, mu).get(a, b));
          col[ix]= j1;
          ix = find_index(i, j2, col, row);
          m[ix]= m[ix]-
            0.5*conjug(V.get(x-e_mu, mu).get(b, a))*BC(x,-e_mu);
          col[ix]= j2;

        }
      }
    }
  }

  // SIMON: susy det term modification
  if (GO) {
    P = Plaq(U);
    sites = 0;
    while(loop_over_lattice(x, sites)) {
      for (mu = 0;mu<NUMLINK;mu++) {
        for (nu = 0;nu<NUMLINK;nu++) {
          if (mu == nu) continue;
          e_mu = Lattice_Vector(mu);

          W = det(P.get(x, mu, nu))-Complex(1.0, 0.0);

          Z = sqrt(1.0*NCOLOR)*Complex(0.0, 1.0)*(Complex(1.0, 0.0)+W);

          i = pack(0, x, NUMGEN-1);

          for (b = 0;b<NUMGEN;b++) {
            j1 = pack(num_chis+1+mu, x, b);
            ix = find_index(i, j1, col, row);
            m[ix]= m[ix]+G*0.5*Z*Tr(inverse(U.get(x, mu))*Lambda[b]);
            col[ix]= j1;

            j1 = pack(num_chis+1+nu, x+e_mu, b);
            ix = find_index(i, j1, col, row);
            m[ix]= m[ix]+G*0.5*Z*Tr(inverse(U.get(x+e_mu, nu))*Lambda[b])*BC(x, e_mu);
            col[ix]= j1;
          }

          j1 = pack(0, x, NUMGEN-1);

          for (b = 0;b<NUMGEN;b++) {
            i = pack(num_chis+1+mu, x, b);
            ix = find_index(i, j1, col, row);
            m[ix]= m[ix]-G*0.5*Z*Tr(inverse(U.get(x, mu))*Lambda[b]);
            col[ix]= j1;

            i = pack(num_chis+1+nu, x+e_mu, b);
            ix = find_index(i, j1, col, row);
            m[ix]= m[ix]-G*0.5*Z*Tr(inverse(U.get(x+e_mu, nu))*Lambda[b])*BC(x, e_mu);
            col[ix]= j1;
          }

        }}}}

  // psi_mu Dplus eta

  sites = 0;
  while(loop_over_lattice(x, sites)) {
    for (mu = 0;mu<NUMLINK;mu++) {
      e_mu = Lattice_Vector(mu);

      for (a = 0;a<NUMGEN;a++) {
        i = pack(mu+num_chis+1, x, a);

        for (b = 0;b<NUMGEN;b++) {
          j1 = pack(0, x+e_mu, b);
          j2 = pack(0, x, b);
          ix = find_index(i, j1, col, row);
          m[ix]= m[ix]+
            0.5*conjug(V.get(x, mu).get(a, b))*BC(x, e_mu);
          col[ix]= j1;
          ix = find_index(i, j2, col, row);
          m[ix]= m[ix]-
            0.5*conjug(V.get(x, mu).get(b, a));
          col[ix]= j2;
        }

      }
    }
  }


  // psi_mu Dminus chi
  sites = 0;
  while(loop_over_lattice(x, sites)) {
    for (nu = 0;nu<NUMLINK;nu++) {
      e_nu = Lattice_Vector(nu);

      for (a = 0;a<NUMGEN;a++) {

        for (mu = nu+1;mu<NUMLINK;mu++) {

          e_mu = Lattice_Vector(mu);

          for (b = 0;b<NUMGEN;b++) {
            i = pack(nu+num_chis+1, x, a);
            j1 = pack(flavor(nu, mu), x, b);
            j2 = pack(flavor(nu, mu), x-e_mu, b);
            ix = find_index(i, j1, col, row);
            m[ix]= m[ix]-
              V.get(x+e_nu, mu).get(a, b);
            col[ix]= j1;
            ix = find_index(i, j2, col, row);
            m[ix]= m[ix]+
              V.get(x-e_mu, mu).get(b, a)*BC(x,-e_mu);
            col[ix]= j2;

            i = pack(mu+num_chis+1, x, a);
            j1 = pack(flavor(nu, mu), x, b);
            j2 = pack(flavor(nu, mu), x-e_nu, b);
            ix = find_index(i, j1, col, row);
            m[ix]= m[ix]+
              V.get(x+e_mu, nu).get(a, b);
            col[ix]= j1;
            ix = find_index(i, j2, col, row);
            m[ix]= m[ix]-
              V.get(x-e_nu, nu).get(b, a)*BC(x,-e_nu);
            col[ix]= j2;
          }
        }
      }

    }
  }

  // chi Dplus psi_mu
  sites = 0;
  while(loop_over_lattice(x, sites)) {
    for (mu = 0;mu<NUMLINK;mu++) {
      for (nu = mu+1;nu<NUMLINK;nu++) {

        e_mu = Lattice_Vector(mu);
        e_nu = Lattice_Vector(nu);

        for (a = 0;a<NUMGEN;a++) {
          i = pack(flavor(mu, nu), x, a);

          for (b = 0;b<NUMGEN;b++) {
            j1 = pack(nu+num_chis+1, x+e_mu, b);
            j2 = pack(nu+num_chis+1, x, b);
            ix = find_index(i, j1, col, row);
            m[ix]= m[ix]+
              V.get(x, mu).get(a, b)*BC(x, e_mu);
            col[ix]= j1;
            ix = find_index(i, j2, col, row);
            m[ix]= m[ix]-
              V.get(x+e_nu, mu).get(b, a);
            col[ix]= j2;

            j1 = pack(mu+num_chis+1, x+e_nu, b);
            j2 = pack(mu+num_chis+1, x, b);
            ix = find_index(i, j1, col, row);
            m[ix]= m[ix]-
              V.get(x, nu).get(a, b)*BC(x, e_nu);
            col[ix]= j1;
            ix = find_index(i, j2, col, row);
            m[ix]= m[ix]+
              V.get(x+e_mu, nu).get(b, a);
            col[ix]= j2;

          }
        }

      }
    }
  }

  // Q-closed pieces
  if (NUMLINK == 5) {
    sites = 0;
    while(loop_over_lattice(x, sites)) {
      for (d = 0;d<NUMLINK;d++) {
        for (e = d+1;e<NUMLINK;e++) {

          for (l = 0;l<NUMGEN;l++) {
            i = pack(flavor(d, e), x, l);

            for (c = 0;c<NUMLINK;c++) {
              if ((c == d)||(c == e)) continue;

              e_c = Lattice_Vector(c);

              for (a = 0;a<NUMLINK;a++) {
                e_a = Lattice_Vector(a);
                if ((a == c)||(a == d)||(a == e)) continue;
                for (b = a+1;b<NUMLINK;b++) {
                  e_b = Lattice_Vector(b);
                  if ((b == c)||(b == d)||(b == e)) continue;

                  for (k = 0;k<NUMGEN;k++) {

                    j1 = pack(flavor(a, b), x-e_a-e_b, k);
                    j2 = pack(flavor(a, b), x-e_a-e_b-e_c, k);

                    ix = find_index(i, j1, col, row);
                    m[ix]= m[ix]+0.5*perm[d][e][c][a][b]*
                      conjug(V.get(x-e_a-e_b-e_c, c).get(l, k))*BC(x,-e_a,-e_b);
                    col[ix]= j1;
                    ix = find_index(i, j2, col, row);
                    m[ix]= m[ix]-0.5*perm[d][e][c][a][b]*
                      conjug(V.get(x-e_c, c).get(k, l))*BC(x,-e_a,-e_b,-e_c);
                    col[ix]= j2;
                  }}
              }

            }}}}}


    sites = 0;
    while(loop_over_lattice(x, sites)) {
      for (a = 0;a<NUMLINK;a++) {
        e_a = Lattice_Vector(a);
        for (b = a+1;b<NUMLINK;b++) {
          e_b = Lattice_Vector(b);

          for (l = 0;l<NUMGEN;l++) {
            i = pack(flavor(a, b), x, l);

            for (c = 0;c<NUMLINK;c++) {
              if ((c == a)||(c == b)) continue;

              e_c = Lattice_Vector(c);

              for (d = 0;d<NUMLINK;d++) {

                if ((d == c)||(d == a)||(d == b)) continue;
                for (e = d+1;e<NUMLINK;e++) {

                  if ((e == c)||(e == a)||(e == b)) continue;
                  for (k = 0;k<NUMGEN;k++) {

                    j1 = pack(flavor(d, e), x+e_a+e_b+e_c, k);
                    j2 = pack(flavor(d, e), x+e_a+e_b, k);

                    ix = find_index(i, j1, col, row);
                    m[ix]= m[ix]+0.5*perm[a][b][c][d][e]*
                      conjug(V.get(x+e_a+e_b, c).get(l, k))*BC(x, e_a, e_b, e_c);
                    col[ix]= j1;
                    ix = find_index(i, j2, col, row);
                    m[ix]= m[ix]-0.5*perm[a][b][c][d][e]*
                      conjug(V.get(x-e_c, c).get(k, l))*BC(x, e_a, e_b);
                    col[ix]= j2;

                  }}
              }

            }}}}}
  }

  // remove blanks ...
  int n2[LEN], c2[LEN*NONZEROES];
  Complex m2[LEN*NONZEROES];

  k = 0;
  l = 0;
  for (i = 0;i<LEN;i++) {
    l+= num_in_row[i];
    n2[i]= 0;
    for (int j = row[i];j<row[i+1];j++) {
      if (col[j]!=(-1)) {
        m2[k]= m[j];
        c2[k]= col[j];
        k++;
        n2[i]++;}
    }
  }

  //cout << "Total number of nonzeroes " << k << "\t" << l << endl;
  TOTALNONZEROES = k;
  row[0]= 0;
  for (i = 0;i<LEN;i++) {
    if (n2[i]!= num_in_row[i]) {cout << "uh oh ..." << endl;}
    row[i+1]= row[i]+n2[i];
  }

  for (i = 0;i<LEN;i++) {
    for (int j = row[i];j<row[i+1];j++) {
      if (c2[j]==(-1)) {cout << "oops i = " << i << "\t" << "j is " << j << endl;
        cout << "row[i] is " << row[i] << "\t" << "row[i+1] is" << row[i+1] << endl;}
      col[j]= c2[j];
      m[j]= m2[j];
    }}



  // order by columns within a row block
  col_order(m, col, row);

  //for (i = 0;i<= LEN;i++) {
  //f_op << row[i] << endl;
  //}
  //for (i = 0;i<(NONZEROES*LEN);i++) {
  //f_op << col[i] << endl;
  //}
  //for (i = 0;i<(NONZEROES*LEN);i++) {
  //f_op << m[i] << endl;
  //}
  //f_op.close();

  if (first_time) {
    cout << "Maximum number nonzeroes in row " << num_in_row[3] << "\n" << flush;
    cout << "Total number of nonzeroes " << TOTALNONZEROES << endl;
    first_time = 0;}

  return;
}

// simple bubble sort for each row
void mysort(Complex d[], int key[], int n) {
  int i, newn;
  Complex c;
  int p;

  do{
    newn = 0;
    for (i = 0;i<n-1;i++) {
      if (key[i]>key[i+1]) {p = key[i];c = d[i];key[i]= key[i+1];d[i]= d[i+1];
        key[i+1]= p;d[i+1]= c;}
      newn = i+1;}

    n = newn;
  }while(n>1);

  return;
}


void col_order(Complex m[], int col[], int row[]) {
  Complex d[LEN];
  int key[LEN];
  int i, j;

  for (i = 0;i<LEN;i++) {

    for (j = 0;j<row[i+1]-row[i];j++) {
      d[j]= m[row[i]+j];
      key[j]= col[row[i]+j];}

    mysort(d, key, row[i+1]-row[i]);

    for (j = 0;j<row[i+1]-row[i];j++) {
      m[row[i]+j]= d[j];
      col[row[i]+j]= key[j];}

  }

  // test row
  //for (i = 0;i<LEN;i++) {
  //for (j = row[i];j<row[i]-1;j++) {
  //if (col[j]>col[j+1]) {cout << "error in col order " << endl;}
  //}}

  return;}


  Complex dot(Complex a[], Complex b[]) {
    Complex dum = Complex(0.0, 0.0);
    for (int i = 0;i<LEN;i++) {
      dum = dum+conjug(a[i])*b[i];
    }
    return(dum);
  }


void sparse_mult(Complex m[], int col[], int row[], int sign, Complex v[], Complex p[]) {
  int i, j;

  if (sign ==(-1)) {
    for (i = 0;i<LEN;i++) {
      v[i]=-1.0*conjug(v[i]);}
  }

  for (i = 0;i<LEN;i++) {
    p[i]= Complex();
    for (j = row[i];j<row[i+1];j++) {
      if (col[j]!=(-1)) {
        p[i]= p[i]+m[j]*v[col[j]];
      }
    }
  }

  if (sign ==(-1)) {
    for (i = 0;i<LEN;i++) {
      p[i]= conjug(p[i]);}
  }

  return;
}

#ifdef FULLMATRIX
void full_fermion_op(const Gauge_Field &U, Complex M[LEN][LEN]) {
  int sites, a, b, j1, j2, i, mu, nu, l, m, d, e, c, index, num_chis;
  Lattice_Vector x, e_mu, e_nu, e_a, e_b, e_c;
  Gauge_Field Udag, Udtr;
  UPlaq_Field P;
  Complex W, Z;

  Udag = Adj(U);


  for (j1 = 0;j1<LEN;j1++) {
    for (j2 = 0;j2<LEN;j2++) {
      M[j1][j2]= Complex();
    }}


  // eta dmubar psi_mu

  num_chis =(NUMLINK*(NUMLINK-1))/2;
  //cout<< "number of chis " << num_chis << "\n" << flush;

  sites = 0;
  while(loop_over_lattice(x, sites)) {
    for (a = 0;a<NUMGEN;a++) {

      i = pack(0, x, a);

      for (mu = 0;mu<NUMLINK;mu++) {
        e_mu = Lattice_Vector(mu);
        for (b = 0;b<NUMGEN;b++) {

          j1 = pack(num_chis+1+mu, x, b);
          j2 = pack(num_chis+1+mu, x-e_mu, b);
          M[i][j1]= M[i][j1]+Tr(Lambda[b]*Udag.get(x, mu)*Lambda[a]);
          M[i][j2]= M[i][j2]-Tr(Lambda[a]*Udag.get(x-e_mu, mu)*Lambda[b])*BC(x,-e_mu);
        }
      }

    }
  }

  // psi_mu dbarT eta
  sites = 0;
  while(loop_over_lattice(x, sites)) {
    for (a = 0;a<NUMGEN;a++) {
      for (mu = 0;mu<NUMLINK;mu++) {
        e_mu = Lattice_Vector(mu);
        i = pack(mu+num_chis+1, x, a);

        for (b = 0;b<NUMGEN;b++) {
          j1 = pack(0, x+e_mu, b);
          j2 = pack(0, x, b);
          M[i][j1]= M[i][j1]+Tr(Lambda[b]*Udag.get(x, mu)*Lambda[a])*BC(x, e_mu);
          M[i][j2]= M[i][j2]-Tr(Lambda[a]*Udag.get(x, mu)*Lambda[b]);
        }}}}

  // chi_{munu} dmu psi_nu
  sites = 0;
  while(loop_over_lattice(x, sites)) {

    for (mu = 0;mu<NUMLINK;mu++) {
      e_mu = Lattice_Vector(mu);

      for (nu = mu+1;nu<NUMLINK;nu++) {
        e_nu = Lattice_Vector(nu);

        index = flavor(mu, nu);
        //cout << "index is " << index << "\n" << flush;
        for (a = 0;a<NUMGEN;a++) {
          i = pack(index, x, a);

          for (b = 0;b<NUMGEN;b++) {

            j1 = pack(num_chis+1+nu, x+e_mu, b);
            j2 = pack(num_chis+1+nu, x, b);
            M[i][j1]= M[i][j1]+Tr(Lambda[a]*U.get(x, mu)*Lambda[b])*BC(x, e_mu);
            M[i][j2]= M[i][j2]-Tr(Lambda[b]*U.get(x+e_nu, mu)*Lambda[a]);
            j1 = pack(num_chis+1+mu, x+e_nu, b);
            j2 = pack(num_chis+1+mu, x, b);
            M[i][j1]= M[i][j1]-Tr(Lambda[a]*U.get(x, nu)*Lambda[b])*BC(x, e_nu);
            M[i][j2]= M[i][j2]+Tr(Lambda[b]*U.get(x+e_mu, nu)*Lambda[a]);

          }}}}}

  // psi_nu dmu chi_{munu}
  sites = 0;
  while(loop_over_lattice(x, sites)) {

    for (nu = 0;nu<NUMLINK;nu++) {
      e_nu = Lattice_Vector(nu);

      for (a = 0;a<NUMGEN;a++) {

        for (mu = nu+1;mu<NUMLINK;mu++) {
          e_mu = Lattice_Vector(mu);

          index = flavor(nu, mu);
          //cout << "index is " << index << "\n" << flush;
          for (b = 0;b<NUMGEN;b++) {

            i = pack(num_chis+1+nu, x, a);
            j1 = pack(index, x, b);
            j2 = pack(index, x-e_mu, b);
            M[i][j1]= M[i][j1]-Tr(Lambda[a]*U.get(x+e_nu, mu)*Lambda[b]);
            M[i][j2]= M[i][j2]+Tr(Lambda[b]*U.get(x-e_mu, mu)*Lambda[a])*BC(x,-e_mu);
            i = pack(num_chis+1+mu, x, a);
            j1 = pack(index, x, b);
            j2 = pack(index, x-e_nu, b);
            M[i][j1]= M[i][j1]+Tr(Lambda[a]*U.get(x+e_mu, nu)*Lambda[b]);
            M[i][j2]= M[i][j2]-Tr(Lambda[b]*U.get(x-e_nu, nu)*Lambda[a])*BC(x,-e_nu);
          }
        }
      }
    }}

  // Q-closed term
  if (NUMLINK == 5) {
    sites = 0;
    while(loop_over_lattice(x, sites)) {
      for (d = 0;d<NUMLINK;d++) {
        for (e = d+1;e<NUMLINK;e++) {

          for (l = 0;l<NUMGEN;l++) {
            i = pack(flavor(d, e), x, l);

            for (c = 0;c<NUMLINK;c++) {
              if ((c == d)||(c == e)) continue;

              e_c = Lattice_Vector(c);

              for (a = 0;a<NUMLINK;a++) {
                e_a = Lattice_Vector(a);
                if ((a == c)||(a == d)||(a == e)) continue;
                for (b = a+1;b<NUMLINK;b++) {
                  e_b = Lattice_Vector(b);
                  if ((b == c)||(b == d)||(b == e)) continue;

                  for (m = 0;m<NUMGEN;m++) {

                    j1 = pack(flavor(a, b), x-e_a-e_b, m);
                    j2 = pack(flavor(a, b), x-e_a-e_b-e_c, m);

                    M[i][j1]= M[i][j1]+0.5*perm[d][e][c][a][b]*
                      Tr(Lambda[m]*Udag.get(x-e_a-e_b-e_c, c)*Lambda[l])*BC(x,-e_a,-e_b);
                    M[i][j2]= M[i][j2]-0.5*perm[d][e][c][a][b]*
                      Tr(Lambda[l]*Udag.get(x-e_c, c)*Lambda[m])*BC(x,-e_a,-e_b,-e_c);

                  }}
              }

            }}}}}


    sites = 0;
    while(loop_over_lattice(x, sites)) {
      for (a = 0;a<NUMLINK;a++) {
        e_a = Lattice_Vector(a);
        for (b = a+1;b<NUMLINK;b++) {
          e_b = Lattice_Vector(b);

          for (l = 0;l<NUMGEN;l++) {
            i = pack(flavor(a, b), x, l);

            for (c = 0;c<NUMLINK;c++) {
              if ((c == a)||(c == b)) continue;

              e_c = Lattice_Vector(c);

              for (d = 0;d<NUMLINK;d++) {

                if ((d == c)||(d == a)||(d == b)) continue;
                for (e = d+1;e<NUMLINK;e++) {

                  if ((e == c)||(e == a)||(e == b)) continue;
                  for (m = 0;m<NUMGEN;m++) {

                    j1 = pack(flavor(d, e), x+e_a+e_b+e_c, m);
                    j2 = pack(flavor(d, e), x+e_a+e_b, m);

                    M[i][j1]= M[i][j1]+0.5*perm[a][b][c][d][e]*
                      Tr(Lambda[m]*Udag.get(x+e_a+e_b, c)*Lambda[l])*BC(x, e_a, e_b, e_c);
                    M[i][j2]= M[i][j2]-0.5*perm[a][b][c][d][e]*
                      Tr(Lambda[l]*Udag.get(x-e_c, c)*Lambda[m])*BC(x, e_a, e_b);

                  }}
              }

            }}}}}
  }


  // SIMON: susy det term modification
  if (NUMGEN ==(NCOLOR*NCOLOR)) {
    P = Plaq(U);
    sites = 0;
    while(loop_over_lattice(x, sites)) {
      for (mu = 0;mu<NUMLINK;mu++) {
        for (nu = 0;nu<NUMLINK;nu++) {
          if (mu == nu) continue;
          e_mu = Lattice_Vector(mu);

          W = det(P.get(x, mu, nu))-Complex(1.0, 0.0);

          Z = sqrt(1.0*NCOLOR)*Complex(0.0, 1.0)*(Complex(1.0, 0.0)+W);
          i = pack(0, x, NUMGEN-1);

          for (b = 0;b<NUMGEN;b++) {
            j1 = pack(num_chis+1+mu, x, b);
            M[i][j1]= M[i][j1]+G*0.5*Z*Tr(inverse(U.get(x, mu))*Lambda[b]);

            j1 = pack(num_chis+1+nu, x+e_mu, b);
            M[i][j1]= M[i][j1]+G*0.5*Z*Tr(inverse(U.get(x+e_mu, nu))*Lambda[b])*BC(x, e_mu);
          }

          j1 = pack(0, x, NUMGEN-1);

          for (b = 0;b<NUMGEN;b++) {
            i = pack(num_chis+1+mu, x, b);
            M[i][j1]= M[i][j1]-G*0.5*Z*Tr(inverse(U.get(x, mu))*Lambda[b]);

            i = pack(num_chis+1+nu, x+e_mu, b);
            M[i][j1]= M[i][j1]-G*0.5*Z*Tr(inverse(U.get(x+e_mu, nu))*Lambda[b])*BC(x, e_mu);
          }

        }}}
  }

  //test antisymmetric nature
  for (j1 = 0;j1<LEN;j1++) {
    for (j2 = j1+1;j2<LEN;j2++) {
      //cout << "M[j1][j2] is " << M[j1][j2] << "\n" <<flush;
      Complex dum = M[j1][j2]+M[j2][j1];
      if (dum.norm()>0.00000001) {
        cout << "problem " << j1 << "\t" << j2 << "\t" <<
          M[j1][j2] << "\t" << M[j2][j1] << "\n" << flush;}
    }
  }

  return;
}


double cnorm(const Complex &c) {
  double dum;
  dum = c.real()*c.real()+c.imag()*c.imag();
  return(sqrt(dum));
}

/*
void eigenvalues(Complex M[LEN][LEN]) {
  static int firsttime = 1;
  static ofstream f_feigen, f_det;
  int i, j, ok, c1, c2, c3;
  char c4;
  Complex temp[2*LEN+1];
  double b[2*LEN], dummy[2], work[4*LEN];
  double at[2*LEN*LEN];

  if (firsttime) {
   f_feigen.open("feigen", ios::app);
   if (f_feigen.bad()) {
   cout << "failed to open feigen file\n" << flush ;}
   f_det.open("detarg", ios::app);
   if (f_det.bad()) {
     cout << "failed to open argdet file\n" << flush;}
     firsttime = 0;
   }

   if (1) {
   for (i = 0;i<LEN;i++) {
   for (j = 0;j<LEN;j++) {
   at[2*(j+LEN*i)]= M[j][i].real();
   at[2*(j+LEN*i)+1]= M[j][i].imag();
   }
   }

   c1 = LEN;
   c2 = 2*LEN;
   c3 = 1;
   c4 ='N';

   zgeev_(&c4,&c4,&c1, at,&c1, b, dummy,&c3, dummy,&c3, work,&c2, work,&ok);

   for (j = 0;j<2*LEN;j = j+2) {
   temp[j/2+1]= Complex(b[j], b[j+1]);}

   ceigsrt(temp, LEN);
   int SUB = 0;

   double arg[LEN+1];
   Complex phase[LEN+1];

   for (i = 1;i<=(LEN-SUB);i++) {
   arg[i]= 0.0;
   for (j = 1;j<= i;j++) {
   arg[i]= arg[i]+atan2(temp[j].imag(), temp[j].real());
   }
   phase[i]= Complex(cos(arg[i]), sin(arg[i]));}


   Complex tmp;
   int count = 0;
   for (i = 1;i<=(LEN);i++) {
   if (temp[i].norm()<0.000001) {count++;}
//f_feigen << i << "\t" << phase[i] << endl;
f_feigen << temp[i] << endl;}

f_feigen << "\n";

cout << "number of zeromodes  is " << count << endl;

f_det << Complex(cos(arg[LEN]), sin(arg[LEN])) << endl;
  }
  return;
}
*/

Complex Pfaffian(Complex M[LEN][LEN]) {
  int i, j, k, jpiv, interchange = 1;
  static int firsttime = 1;
  double pivot, totalangle, angle, cosine, sine, mag;
  Complex dum, scale, f;
  //        cout << "in Pfaffian\n" << flush;
  static ofstream f_pfaff;

  if (firsttime) {
    f_pfaff.open("pfaffian", ios::app);
    if (f_pfaff.bad()) {
      cout << "failed to open pfaffian file\n" << flush ;}
    firsttime = 0;
  }

  // Loop over all rows in steps of 2
  for (i = 0; i < LEN - 2; i += 2) {
    // Find column whose ith component is biggest to use as pivot
    pivot = cnorm(M[i][i + 1]);
    jpiv = i + 1;
    for (j = i + 2; j < LEN; j++) {
      if (cnorm(M[i][j]) > pivot) {
        pivot = cnorm(M[i][j]);
        jpiv = j;
      }
    }

    // If necessary, interchange M[j][i + 1] <--> M[j][jpiv]
    //                       and M[i + 1][j] <--> M[jpiv][j]
    if (jpiv != i + 1) {
      interchange *= -1;
      // Interchange column(i + 1) with column(jpiv)
      for (j = 0; j < LEN; j++) {
        dum = M[j][i + 1];
        M[j][i + 1] = M[j][jpiv];
        M[j][jpiv] = dum;
      }
      // Interchange row(i + 1) with row(jpiv)
      for (j = 0; j < LEN; j++) {
        dum = M[i + 1][j];
        M[i + 1][j] = M[jpiv][j];
        M[jpiv][j] = dum;
      }
    }

    // Progressively zero out elements j = i + 2, ..., LEN - 1 of row M[i]
    // Also zero out elements along the corresponding column M[:][i]
    for (j = i + 2; j < LEN; j++) {
      scale = M[i][j] / M[i][i + 1];
      for (k = 0; k < LEN; k++) {
        M[k][j] = M[k][j] - scale * M[k][i + 1];
        M[j][k] = M[j][k] - scale * M[i + 1][k];
      }
    }

    // Progressively zero out elements j = i + 2, ..., LEN - 1 of row M[i + 1]
    // Also zero out elements along the corresponding column M[:][i + 1]
    for (j = i + 2; j < LEN; j++) {
      scale = M[i + 1][j] / M[i + 1][i];
      for (k = 0; k < LEN; k++) {
        M[k][j] = M[k][j] - scale * M[k][i];
        M[j][k] = M[j][k] - scale * M[i][k];
      }
    }
  }

  f = Complex(1.0, 0.0);
  mag = 0.0;
  totalangle = 0.0;
  double TWOPI = 4.0 * acos(0.0);
  for (i = 0; i < LEN; i += 2) {
    cosine = M[i][i + 1].real() / M[i][i + 1].norm();
    sine = M[i][i + 1].imag() / M[i][i + 1].norm();
    if (cosine >= 0.0) {
      if (sine >= 0.0)
        angle = acos(cosine);
      else
        angle = TWOPI - acos(cosine);
    }
    else {    // cosine < 0.0
      if (sine >= 0.0)
        angle = acos(cosine);
      else
        angle = TWOPI - acos(cosine);
    }
    mag += log(M[i][i + 1].norm());
    totalangle += angle;
    f = f * M[i][i + 1];
  }

  cout << "Pfaffian is " << Complex(cos(totalangle)*interchange, sin(totalangle)*interchange) << endl;

  f_pfaff << mag << "\t" << cos(totalangle)*interchange << "\t"
    << sin(totalangle)*interchange << "\n" << flush;

  return (f * interchange);
}
#endif

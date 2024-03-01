#include "dm_ead.h"
#include "../d_basis_functions.h"
#include <vector>
#include <iomanip>


/* Index ii iterates rhoi array
 * rhoi array is split 50/50 into density and derivative
 * of the embedding function
 */
DM_EAD::DM_EAD(Config &config):
    verbose(config.get<int>("VERBOSE")),
    rc(config.get<double>("RCUTMB"))
{

    if (!config.get<bool>("INITMB")) return;

    get_grid(config,"CGRIDMB",cgrid);
    get_grid(config,"SGRIDMB",sgrid);

    if (cgrid.size()!=sgrid.size()) {
        throw std::runtime_error("SGRID2B and CGRID2B arrays differ in size.\n");
    }

    Lmax = config.get<int>("AGRIDMB");
    s=sgrid.size()*(Lmax+1);
    rhoisize = gen_atomic_orbitals(Lmax);
    rhoisize *= sgrid.size();

    if (verbose) {
        std::cout << std::endl;
        std::cout << "SGRID (SGRIDMB): " << sgrid.size() << std::endl;
        for (auto e:sgrid) std::cout << e << "  ";
        std::cout << std::endl;

        std::cout << "CGRID (CGRIDMB): " << cgrid.size() << std::endl;
        for (auto L:cgrid) std::cout << L << "  ";
        std::cout << std::endl;
        std::cout << "rhoisize: " << rhoisize << std::endl;
    }
}

void DM_EAD::calc_aed(
        rho_type& rhoi,
        aed_type2 &aed)
{
    size_t Lshift = 0;
    size_t ii=0;
    for (int L=0; L<=Lmax; ++L) {
        for (size_t o=0; o<orbitals[L].size(); ++o) {
            const double f = fac[L][o];
            for (size_t c=0; c<cgrid.size(); c++) {
                aed(c+Lshift+fidx) += f*rhoi(ii)*rhoi(ii);
                rhoi(ii+rhoisize) = 2.0*f*rhoi(ii);
                ii++;
            }
        }
        Lshift += sgrid.size();
    }
}
int DM_EAD::calc_dXijdri_dXjidri(
        const double rij,
        const double,
        const Vec3d &vec_ij,
        const double fc_ij,
        const double fcp_ij,
        rho_type& rhoi,
        rho_type& rhoj,
        fd_type &fd_ij,
        const double wi,
        const double wj)
{
    const double x = vec_ij[0];
    const double y = vec_ij[1];
    const double z = vec_ij[2];

    size_t Lshift = 0;

    size_t ii=rhoisize;
    for (int L=0; L<=Lmax; ++L) {
        for (size_t o=0; o<orbitals[L].size(); ++o) {
            const int lx = orbitals[L][o][0];
            const int ly = orbitals[L][o][1];
            const int lz = orbitals[L][o][2];

            const double powxlxij = my_pow(x,lx);
            const double powylyij = my_pow(y,ly);
            const double powzlzij = my_pow(z,lz);

            const double powxlxji = my_pow(-x,lx);
            const double powylyji = my_pow(-y,ly);
            const double powzlzji = my_pow(-z,lz);

            double txij=0.0,tyij=0.0,tzij=0.0;
            double txji=0.0,tyji=0.0,tzji=0.0;
            if (lx!=0) {
                txij = lx*my_pow(x,lx-1)*powylyij*powzlzij;
                txji = lx*my_pow(-x,lx-1)*powylyji*powzlzji;
            }
            if (ly!=0) {
                tyij = ly*my_pow(y,ly-1)*powxlxij*powzlzij;
                tyji = ly*my_pow(-y,ly-1)*powxlxji*powzlzji;
            }
            if (lz!=0) {
                tzij = lz*my_pow(z,lz-1)*powylyij*powxlxij;
                tzji = lz*my_pow(-z,lz-1)*powylyji*powxlxji;
            }
            const double pqrij = powxlxij*powylyij*powzlzij;
            const double pqrji = powxlxji*powylyji*powzlzji;

            for (size_t c=0; c<cgrid.size(); c++) {
                double gauss = G(rij,sgrid[c],cgrid[c],fc_ij);
                double dgauss = dG(rij,sgrid[c],cgrid[c],fc_ij,fcp_ij);

                const double term2ijx=pqrij*dgauss*x/rij;
                const double term2ijy=pqrij*dgauss*y/rij;
                const double term2ijz=pqrij*dgauss*z/rij;

                const double term2jix=-pqrji*dgauss*x/rij;
                const double term2jiy=-pqrji*dgauss*y/rij;
                const double term2jiz=-pqrji*dgauss*z/rij;

                const double term1ijx = gauss*txij;
                const double term1ijy = gauss*tyij;
                const double term1ijz = gauss*tzij;

                const double term1jix = gauss*txji;
                const double term1jiy = gauss*tyji;
                const double term1jiz = gauss*tzji;

                fd_ij(c+Lshift+fidx,0) +=
                    rhoi(ii)*(term1ijx + term2ijx)*wj
                    -rhoj(ii)*(term1jix + term2jix)*wi
                    ;

                fd_ij(c+Lshift+fidx,1) +=
                    rhoi(ii)*(term1ijy + term2ijy)*wj
                    -rhoj(ii)*(term1jiy + term2jiy)*wi
                    ;

                fd_ij(c+Lshift+fidx,2) += 
                    rhoi(ii)*(term1ijz + term2ijz)*wj
                    -rhoj(ii)*(term1jiz + term2jiz)*wi
                    ;
                ii++;
            }
        }
        Lshift+=sgrid.size();
    }

    return 1;

}
int DM_EAD::calc_dXijdri(
        const double rij,
        const double,
        const Vec3d &vec_ij,
        const double fc_ij,
        const double fcp_ij,
        rho_type& rhoi,
        fd_type &fd_ij)
{
    const double x = vec_ij[0];
    const double y = vec_ij[1];
    const double z = vec_ij[2];

    size_t Lshift = 0;

    size_t ii=rhoisize;
    for (int L=0; L<=Lmax; ++L) {
        for (size_t o=0; o<orbitals[L].size(); ++o) {
            const int lx = orbitals[L][o][0];
            const int ly = orbitals[L][o][1];
            const int lz = orbitals[L][o][2];

            const double powxlxij = my_pow(x,lx);
            const double powylyij = my_pow(y,ly);
            const double powzlzij = my_pow(z,lz);

            double txij=0.0,tyij=0.0,tzij=0.0;
            if (lx!=0) {
                txij = lx*my_pow(x,lx-1)*powylyij*powzlzij;
            }
            if (ly!=0) {
                tyij = ly*my_pow(y,ly-1)*powxlxij*powzlzij;
            }
            if (lz!=0) {
                tzij = lz*my_pow(z,lz-1)*powylyij*powxlxij;
            }
            const double pqrij = powxlxij*powylyij*powzlzij;

            for (size_t c=0; c<cgrid.size(); c++) {
                double gauss = G(rij,sgrid[c],cgrid[c],fc_ij);
                double dgauss = dG(rij,sgrid[c],cgrid[c],fc_ij,fcp_ij);

                const double term2ijx=pqrij*dgauss*x/rij;
                const double term2ijy=pqrij*dgauss*y/rij;
                const double term2ijz=pqrij*dgauss*z/rij;

                const double term1ijx = gauss*txij;
                const double term1ijy = gauss*tyij;
                const double term1ijz = gauss*tzij;

                fd_ij(c+Lshift+fidx,0) +=
                    rhoi(ii)*(term1ijx + term2ijx)
                    ;

                fd_ij(c+Lshift+fidx,1) +=
                    rhoi(ii)*(term1ijy + term2ijy)
                    ;

                fd_ij(c+Lshift+fidx,2) += 
                    rhoi(ii)*(term1ijz + term2ijz)
                    ;
                ii++;
            }
        }
        Lshift+=sgrid.size();
    }

    return 1;

}
size_t DM_EAD::size() {
    return s;
}
std::string DM_EAD::label() {
    return lab;
}

void DM_EAD::init_rhoi(rho_type& rhoi)
{
    rhoi.resize(2*rhoisize);
    rhoi.set_zero();
}
void DM_EAD::calc_rho(
        const double rij,
        const double,
        const double fc_ij,
        const Vec3d &vec_ij,
        rho_type& rhoi)
{
    size_t ii=0;
    for (int L=0; L<=Lmax; ++L) {
        for (size_t o=0; o<orbitals[L].size(); ++o) {
            const size_t lx = orbitals[L][o][0];
            const size_t ly = orbitals[L][o][1];
            const size_t lz = orbitals[L][o][2];
            double t = (my_pow(vec_ij[0],lx)*my_pow(vec_ij[1],ly)*my_pow(vec_ij[2],lz));
            for (size_t c=0; c<cgrid.size(); c++) {
                rhoi(ii++) += G(rij,sgrid[c],cgrid[c],fc_ij)*t;
            }
        }
    }
}

size_t DM_EAD::gen_atomic_orbitals(int Lmax) {
    orbitals.resize(Lmax+1);
    fac.resize(Lmax+1);

    size_t count = 0;
    for (int L = 0; L <= Lmax; ++L) {
        for (int lx = 0; lx <= L; ++lx) {
            for (int ly = 0; ly <= L; ++ly) {
                for (int lz = 0; lz <= L; ++lz) {
                    if (lx+ly+lz == L) {
                        std::vector<int> o = {lx, ly, lz};
                        orbitals[L].push_back(o);
                        const double f = fact(L)/(fact(lx)*fact(ly)*fact(lz))/my_pow(4*rc,L);
                        fac[L].push_back(f);
                        count++;
                    }
                }
            }
        }
    }
    return count;
}
double DM_EAD::fact(int n)
{
    double f=1.0;
    while (n) {
        f *= n;
        --n;
    }
    return f;
}
double DM_EAD::my_pow(double x, int n){
    double r = 1.0;

    while(n > 0){
        r *= x;
        --n;
    }
    return r;
}
size_t DM_EAD::rhoi_size() {
    return rhoisize;
}
size_t DM_EAD::rhoip_size() {
    return rhoisize;
}

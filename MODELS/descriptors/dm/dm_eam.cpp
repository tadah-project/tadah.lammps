#include "dm_eam.h"


DM_EAM::DM_EAM(Config &c):
    verbose(c.get<int>("VERBOSE")) 
{
    if (!c.get<bool>("INITMB")) {
        s=0;
        return;
    }
    ef.file_path = c("SETFL")[0];
    read_setfl();

    if (std::abs(ef.rcut - c.get<double>("RCUTMB")) > std::numeric_limits<double>::min()) {
        if (verbose) {
            std::cout << "Config file cutoff and setfl file cutoff differ: "
                << c.get<double>("RCUTMB") << " " << ef.rcut << std::endl;
            std::cout << "Enforcing SETFL file cutoff: " << ef.rcut << std::endl;
        }
        c.remove("RCUTMB");
        c.add("RCUTMB", ef.rcut);
        c.postprocess();
    }


    frho_spline.resize(ef.nrho+1, std::vector<double>(7));
    rhor_spline.resize(ef.nr+1, std::vector<double>(7));
    z2r_spline.resize(ef.nr+1, std::vector<double>(7));


    gen_splines(ef.nrho, ef.drho, ef.frho, frho_spline);
    gen_splines(ef.nr, ef.dr, ef.rhor, rhor_spline);
    gen_splines(ef.nr, ef.dr, ef.z2r, z2r_spline);

}

void DM_EAM::calc_aed(
		rho_type& rho,
		aed_type2 &aed)
{
    //std::cout << "rho[i]:   " << rho(0) << std::endl;
    double p = rho(0)*ef.rdrho + 1.0;
    int m = static_cast<int> (p);
    m = std::max(1,std::min(m,ef.nrho-1));
    p -= m;
    p = std::min(p,1.0);
    //std::vector<double> coeff = frho_spline[m];
    double phi = ((frho_spline[m][3]*p + frho_spline[m][4])*p + frho_spline[m][5])*p + frho_spline[m][6];

    rho(1) = (frho_spline[m][0]*p + frho_spline[m][1])*p + frho_spline[m][2];    // lammps fp[i]
    //std::cout << "fp[i]:   " << rho(1) << std::endl;

    if (rho(0) > ef.rhomax) phi += rho(1) * (rho(0)-ef.rhomax);
    //phi *= scale[type[i]][type[i]];

    aed(fidx) += phi;


}
int DM_EAM::calc_dXijdri_dXjidri(
		const double rij,
		const double,
		const Vec3d &,
		const double,
		const double,
		rho_type& rhoi,
		rho_type& rhoj,
		fd_type &fd_ij,
        const double,
        const double)
{
    double r = rij;
    //double recip = 1.0/r;
    double p = r*ef.rdr + 1.0;
    int m = static_cast<int> (p);
    m = std::min(m,ef.nr-1);
    p -= m;
    p = std::min(p,1.0);
    //std::vector<double> coeff = rhor_spline[m];
    double rhoip = (rhor_spline[m][0]*p + rhor_spline[m][1])*p + rhor_spline[m][2];
    double rhojp = (rhor_spline[m][0]*p + rhor_spline[m][1])*p + rhor_spline[m][2];

    // here rho_i[1] is a derivative of rho_i at local atom, similarly for rho_j
    //double psip = rho_i[1]*rhojp + rho_j[1]*rhoip;
    double psip = rhoi(1)*rhojp + rhoj(1)*rhoip;

    fd_ij(fidx,0) = psip; // requires *delij*recip in the main loop
    return 0;

}
int DM_EAM::calc_dXijdri(
		const double rij,
		const double,
		const Vec3d &,
		const double,
		const double,
		rho_type& rhoi,
		fd_type &fd_ij)
{
    double r = rij;
    //double recip = 1.0/r;
    double p = r*ef.rdr + 1.0;
    int m = static_cast<int> (p);
    m = std::min(m,ef.nr-1);
    p -= m;
    p = std::min(p,1.0);
    double rhojp = (rhor_spline[m][0]*p + rhor_spline[m][1])*p + rhor_spline[m][2];

    // here rho_i[1] is a derivative of rho_i at local atom
    double psip = rhoi(1)*rhojp;

    fd_ij(fidx,0) = psip; // requires *delij*recip in the main loop
    return 0;

}
size_t DM_EAM::size() {
	return s;
}
std::string DM_EAM::label() {
	return lab;
}

void DM_EAM::read_setfl()
{
    std::string line;
    std::ifstream in_file(ef.file_path);
    if (!in_file.good())
        throw std::runtime_error("SETFL file does not exists.\n");

    if (in_file.is_open()) {
        if (verbose)
            std::cout << "<DM_EAM> Reading setfl: " << ef.file_path << std::endl;
        // skip ficgridt 3 comment lines
        getline (in_file,line);
        getline (in_file,line);
        getline (in_file,line);
        // skip number of types and types
        getline (in_file,line);
        // read 5th line
        in_file >> ef.nrho >> ef.drho >> ef.nr >> ef.dr >> ef.rcut;
        in_file >> ef.atomic_number >> ef.atomic_mass >> ef.lattice_param >> ef.lattice;
        ef.rdr = 1.0/ef.dr;
        ef.rdrho = 1.0/ef.drho;
        // prepare arrays
        ef.frho.resize(ef.nrho);
        ef.rhor.resize(ef.nr);
        ef.z2r.resize(ef.nr);
        // read all data
        for (int i=0; i<ef.nrho; ++i) in_file >> ef.frho[i];
        for (int i=0; i<ef.nr; ++i) in_file >> ef.rhor[i];
        for (int i=0; i<ef.nr; ++i) in_file >> ef.z2r[i];
        in_file.close();
    }
    else {
        std::cout << "<DM_EAM> Unable to open file: " << ef.file_path << std::endl;;
    }
    // get max value of rho
    //ef.rhomax = *std::max_element(ef.rhor.begin(), ef.rhor.end());
    ef.rhomax = (ef.nrho-1) * ef.drho;

}
void DM_EAM::gen_splines(int &n, double &delta, std::vector<double> &f,
		std::vector<std::vector<double>> &spline)
{
	// in lammps f is n+1, here is size n
	for (int m=1; m<=n; m++) spline[m][6] = f[m-1];

	spline[1][5] = spline[2][6] - spline[1][6];
	spline[2][5] = 0.5 * (spline[3][6]-spline[1][6]);
	spline[n-1][5] = 0.5 * (spline[n][6]-spline[n-2][6]);
	spline[n][5] = spline[n][6] - spline[n-1][6];

	for (int m = 3; m <= n-2; m++)
		spline[m][5] = ((spline[m-2][6]-spline[m+2][6]) +
				8.0*(spline[m+1][6]-spline[m-1][6])) / 12.0;

	for (int m = 1; m <= n-1; m++) {
		spline[m][4] = 3.0*(spline[m+1][6]-spline[m][6]) - 2.0*spline[m][5] - spline[m+1][5];
		spline[m][3] = spline[m][5] + spline[m+1][5] - 2.0*(spline[m+1][6]-spline[m][6]);
	}

	spline[n][4] = 0.0;
	spline[n][3] = 0.0;

	for (int m = 1; m <= n; m++) {
		spline[m][2] = spline[m][5]/delta;
		spline[m][1] = 2.0*spline[m][4]/delta;
		spline[m][0] = 3.0*spline[m][3]/delta;
	}
}
void DM_EAM::init_rhoi(rho_type& rhoi)
{
    rhoi.resize(2);
    rhoi.set_zero();
}
void DM_EAM::calc_rho(
		const double rij,
        const double,
        const double,
        const Vec3d &,
        rho_type& rho)
{
    double r = rij;
    double p = r*ef.rdr + 1.0;
    int m = static_cast<int> (p);
    m = std::min(m,ef.nr-1);
    p -= m;
    p = std::min(p,1.0);
    std::vector<double> coeff = rhor_spline[m];
    rho(0) += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

}
size_t DM_EAM::rhoi_size() {
    return 1;
}
size_t DM_EAM::rhoip_size() {
    return 1;
}

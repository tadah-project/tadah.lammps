#include "d2_eam.h"


D2_EAM::D2_EAM(Config &c):
    verbose(c.get<int>("VERBOSE"))
{
    if (!c.get<bool>("INIT2B")) {
        s=0;
        return;
    }
	ef.file_path = c("SETFL")[0];
	read_setfl();

	if (std::abs(ef.rcut - c.get<double>("RCUT2B")) > std::numeric_limits<double>::min()) {
        if (verbose) {
            std::cout << std::endl;
            std::cout << "Config file cutoff and setfl file cutoff differ: "
                << c.get<double>("RCUT2B") << " " << ef.rcut << std::endl;
            std::cout << "Enforcing SETFL file cutoff: " << ef.rcut << std::endl;
        }
		c.remove("RCUT2B");
		c.add("RCUT2B", ef.rcut);
        c.postprocess();
	}


	frho_spline.resize(ef.nrho+1, std::vector<double>(7));
	rhor_spline.resize(ef.nr+1, std::vector<double>(7));
	z2r_spline.resize(ef.nr+1, std::vector<double>(7));


	gen_splines(ef.nrho, ef.drho, ef.frho, frho_spline);
	gen_splines(ef.nr, ef.dr, ef.rhor, rhor_spline);
	gen_splines(ef.nr, ef.dr, ef.z2r, z2r_spline);

}

void D2_EAM::calc_aed(
		const double rij,
		const double,
		const double,
		aed_type2 &aed)
{

	double r = rij;
	const double recip = 1.0/r;
	double p = r*ef.rdr + 1.0;
	int m = static_cast<int> (p);
	m = std::min(m,ef.nr-1);
	p -= m;
	p = std::min(p,1.0);
	double z2 = ((z2r_spline[m][3]*p + z2r_spline[m][4])*p + z2r_spline[m][5])*p + z2r_spline[m][6];
	aed(fidx) += z2*recip;

}
void D2_EAM::calc_dXijdri(
		const double rij,
		const double,
		const double,
		const double,
		fd_type &fd_ij)
{

	const double r = rij;
	const double recip = 1.0/r;
	double p = r*ef.rdr + 1.0;
	int m = static_cast<int> (p);
	m = std::min(m,ef.nr-1);
	p -= m;
	p = std::min(p,1.0);
	double z2p = (z2r_spline[m][0]*p + z2r_spline[m][1])*p + z2r_spline[m][2];
	double z2 = ((z2r_spline[m][3]*p + z2r_spline[m][4])*p + z2r_spline[m][5])*p + z2r_spline[m][6];
	//double phi = z2*recip;
	//double phip = z2p*recip - z2*recip*recip;

	fd_ij(fidx,0) = (z2p*recip - z2*recip*recip);
	//force_fp_ij[first_idx] = 0.5*(z2p*recip - z2*recip*recip);

}
void D2_EAM::calc_all(
		const double rij,
		const double,
		const double,
		const double,
		aed_type2 &aed,
		fd_type &fd_ij)
{
	double r = rij;
	const double recip = 1.0/r;
	double p = r*ef.rdr + 1.0;
	int m = static_cast<int> (p);
	m = std::min(m,ef.nr-1);
	p -= m;
	p = std::min(p,1.0);
	double z2 = ((z2r_spline[m][3]*p + z2r_spline[m][4])*p + z2r_spline[m][5])*p + z2r_spline[m][6];
	aed(fidx) += z2*recip;

	double z2p = (z2r_spline[m][0]*p + z2r_spline[m][1])*p + z2r_spline[m][2];
	// 0.5 b/c full neighbour list is used TODO check
	fd_ij(fidx,0) = (z2p*recip - z2*recip*recip);


}
size_t D2_EAM::size() {
	return s;
}
std::string D2_EAM::label() {
	return lab;
}

void D2_EAM::read_setfl()
{
    std::string line;
    std::ifstream in_file(ef.file_path);
    if (!in_file.good())
        throw std::runtime_error("SETFL file does not exists.\n");

    if (in_file.is_open()) {
        if (verbose)
            std::cout << std::endl << "<D2_EAM> Reading setfl: " << ef.file_path << std::endl;
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
        if (verbose) std::cout << "<D2_EAM> Unable to open file: " << ef.file_path << std::endl;;
    }
}
void D2_EAM::gen_splines(int &n, double &delta, std::vector<double> &f,
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


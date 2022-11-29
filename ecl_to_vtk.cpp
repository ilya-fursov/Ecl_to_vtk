#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cassert>

//------------------------------------------------------------------------------------------
class Exception : public std::exception
{
protected:
	std::string msg;
public:
	Exception(std::string s) : msg(s){};
	~Exception() noexcept {};
	const char *what() const noexcept {return msg.c_str();};
};
//------------------------------------------------------------------------------------------
long StoL(std::string s, bool &complete);		// complete = true if whole string is a number
double StoD(std::string s, bool &complete);		// complete = true if whole string is a number
long StoL(std::string s);		// throws exception if whole string is not a number
double StoD(std::string s);		// throws exception if whole string is not a number
//------------------------------------------------------------------------------------------
void ReadGrids(const char *file, size_t len, std::vector<std::vector<double>> &data, std::vector<std::string> S1, std::string S2);	// reads a number of grids from file
															// allocates and fills "data" of size S1.size(), with data[i].size() = len
															// S1[i], S2 - start and end of grid[i] which is loaded to data[i]
void ReadGrid(const char *file, size_t len, std::vector<double> &data, std::string S1, std::string S2);		// same as above; reads one grid
bool ReadTokenComm(FILE *F, char **str, bool &new_line, char *str0, const int str0_len);
															// reads a token from the file (delimited by ' ', '\t', '\r', '\n'), dropping "--..." comments
															// returns true on success, false on failure/EOF
															// the token is saved to str
															// set new_line = true in the first call, then the function will manage it
															// str0 is a working array, it should have been allocated
int StrIndex(const std::string &S, const std::vector<std::string> &VECS);	// index of S in VECS[], -1 if not found
inline bool scan_two(const char *str, size_t &cnt, double &d, bool &expect_scan_two);	// parses "cnt*d", returns 'true' on success, updates 'expect_scan_two'
inline bool scan_one(const char *str, double &d, bool &expect_scan_two);				// parses "d", returns 'true' on success, updates 'expect_scan_two'
bool is_active(size_t cell_ind, int poro_ind, int act_ind, const std::vector<std::vector<double>> &props);
//------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	try
	{
		if (argc < 9 || argc%2 == 0)
		{
			std::cerr << "This is a program for converting eclipse-style corner point grid to the legacy VTK ascii format\n"
					 	 "	Expected command line arguments:\n"
						 "	outfile.vtk NX NY NZ COORD_file ZCORN_file PROP1_file PROP1 [PROPi_file PROPi]\n"
						 "	Where outfile.vtk - output file name\n"
						 "	NX, NY, NZ - grid dimensions\n"
						 "	XXX_file - file names for ZCORN, COORD (defining the mesh) and different properties (cell attributes)\n"
						 "	PROP1,..PROPi,.. - corresponding property names\n";
			return 1;
		}

		const size_t Nx = StoL(argv[2]);
		const size_t Ny = StoL(argv[3]);
		const size_t Nz = StoL(argv[4]);

		const int Nprop = (argc-7)/2;
		assert(Nprop >= 1);

		// main input data
		std::vector<double> coord;
		std::vector<double> zcorn;
		std::vector<std::vector<double>> Props(Nprop);	// different grids with properties

		int poro_ind = -1;	// indices of PORO and ACTNUM within "Props", -1 if not present
		int act_ind = -1;

		const size_t coord_size = 6*(Nx+1)*(Ny+1);
		const size_t zcorn_size = 8*Nx*Ny*Nz;
		const size_t prop_size = Nx*Ny*Nz;

		// read the input files
		ReadGrid(argv[5], coord_size, coord, "COORD", "/");
		std::cout << "Loaded COORD\n";

		ReadGrid(argv[6], zcorn_size, zcorn, "ZCORN", "/");
		std::cout << "Loaded ZCORN\n";

		for (int i = 0; i < Nprop; i++)
		{
			ReadGrid(argv[7+2*i], prop_size, Props[i], argv[8+2*i], "/");
			std::cout << "Loaded " << argv[8+2*i] << "\n";

			if ((std::string)argv[8+2*i] == "PORO")
				poro_ind = i;
			if ((std::string)argv[8+2*i] == "ACTNUM")
				act_ind = i;
		}
														// TODO NOTE: inactive cells will not be output
		std::vector<bool> active(prop_size, false);		// indicates if cell is active
		size_t count_act = 0;
		for (size_t c = 0; c < prop_size; c++)
		{
			active[c] = is_active(c, poro_ind, act_ind, Props);
			if (active[c])
				count_act++;
		}

		// calculate vertex coordinates for each cell
		std::vector<double> cell_coor(zcorn_size*3);	// ORDER: (x,y,z) for 8 vertices of the 1st cell, (x,y,z) for 8 vertices of the second cell,...
														// Vertex order in a cell: as in VTK_VOXEL
														// CELLS: i - fastest, k - slowest
		for (size_t k = 0; k < Nz; k++)
			for (size_t j = 0; j < Ny; j++)
				for (size_t i = 0; i < Nx; i++)			// consider cell (i, j, k)
				{
					size_t p[4];
					p[0] = j*(Nx+1) + i;				// global indices of the four pillars
					p[1] = j*(Nx+1) + i+1;
					p[2] = (j+1)*(Nx+1) + i;
					p[3] = (j+1)*(Nx+1) + i+1;
					assert(p[3]*6 + 5 < coord_size);

					size_t v[8];
					v[0] = 2*i + 4*Nx*j + 8*Nx*Ny*k;	// global indices in "zcorn" of the vertices
					v[1] = 2*i+1 + 4*Nx*j + 8*Nx*Ny*k;
					v[2] = 2*(i+Nx) + 4*Nx*j + 8*Nx*Ny*k;
					v[3] = 2*(i+Nx)+1 + 4*Nx*j + 8*Nx*Ny*k;

					v[4] = 2*i + 4*Nx*(j+Ny) + 8*Nx*Ny*k;
					v[5] = 2*i+1 + 4*Nx*(j+Ny) + 8*Nx*Ny*k;
					v[6] = 2*(i+Nx) + 4*Nx*(j+Ny) + 8*Nx*Ny*k;
					v[7] = 2*(i+Nx)+1 + 4*Nx*(j+Ny) + 8*Nx*Ny*k;
					assert(v[7] < zcorn_size);

					double xpill[8];				// order: 4 pillars x_up, 4 pillars x_down
					double ypill[8];
					double zpill[8];

					bool use_z[8];					// if pillar is not vertical, z values from COORD will be used
					for (int n = 0; n < 4; n++)		// n - pillar number
					{
						use_z[n] = use_z[n+4] = true;

						xpill[n] = coord[p[n]*6];
						ypill[n] = coord[p[n]*6+1];
						zpill[n] = coord[p[n]*6+2];

						xpill[n+4] = coord[p[n]*6+3];
						ypill[n+4] = coord[p[n]*6+4];
						zpill[n+4] = coord[p[n]*6+5];

						if (xpill[n] == xpill[n+4] && ypill[n] == ypill[n+4])
							use_z[n] = use_z[n+4] = false;
					}

					// fill the final cell vertices
					for (int n = 0; n < 4; n++)		// n - vertex in the upper plane, n+4 - in the lower plane
					{
						size_t ind = Nx*Ny*k + Nx*j + i;
						const double x0 = xpill[n];
						const double y0 = ypill[n];
						const double z0 = zpill[n];
						const double x1 = xpill[n+4];
						const double y1 = ypill[n+4];
						const double z1 = zpill[n+4];

						if (use_z[n])
						{
							cell_coor[24*ind + 3*n] = x0 + (x1-x0)/(z1-z0)*(zcorn[v[n]] - z0);
							cell_coor[24*ind + 3*n+1] = y0 + (y1-y0)/(z1-z0)*(zcorn[v[n]] - z0);

							cell_coor[24*ind + 3*n + 12] = x0 + (x1-x0)/(z1-z0)*(zcorn[v[n+4]] - z0);
							cell_coor[24*ind + 3*n + 13] = y0 + (y1-y0)/(z1-z0)*(zcorn[v[n+4]] - z0);
						}
						else
						{
							cell_coor[24*ind + 3*n] = x0;
							cell_coor[24*ind + 3*n+1] = y0;

							cell_coor[24*ind + 3*n + 12] = x0;
							cell_coor[24*ind + 3*n + 13] = y0;
						}
						cell_coor[24*ind + 3*n+2] = zcorn[v[n]];
						cell_coor[24*ind + 3*n+14] = zcorn[v[n+4]];
					}
				}

		// output to file
		FILE *fout = fopen(argv[1], "w");
		if (fout == NULL)
			throw Exception((std::string)"Cannot open file for writing: " + argv[1]);

		fprintf(fout, "# vtk DataFile Version 2.0\n");
		fprintf(fout, "Unstructured Grid from COORD & ZCORN\n");
		fprintf(fout, "ASCII\n");
		fprintf(fout, "DATASET UNSTRUCTURED_GRID\n\n");
		fprintf(fout, "POINTS %zu float\n", count_act*8);

		for (size_t i = 0; i < zcorn_size; i++)			// i <-> vertex
			if (active[i/8])
				fprintf(fout, "%-14.9g\t%-14.9g\t%-14.9g\n", cell_coor[i*3], cell_coor[i*3+1], cell_coor[i*3+2]);

		fprintf(fout, "\nCELLS %zu %zu\n", count_act, count_act*9);
		size_t count_1 = 0;
		for (size_t i = 0; i < prop_size; i++)
			if (active[i])
			{
				fprintf(fout, "8");
				for (size_t j = 0; j < 8; j++)
				{
					fprintf(fout, "\t%zu", count_1);
					count_1++;
				}
				fprintf(fout, "\n");
			}

		fprintf(fout, "\nCELL_TYPES %zu\n", count_act);
		for (size_t i = 0; i < count_act; i++)
			fprintf(fout, "11\n");

		fprintf(fout, "\nCELL_DATA %zu\n", count_act);
		for (int p = 0; p < Nprop; p++)
		{
			fprintf(fout, "\nSCALARS %s float 1\n", argv[8+2*p]);
			fprintf(fout, "LOOKUP_TABLE default\n");
			for (size_t i = 0; i < prop_size; i++)
				if (active[i])
					fprintf(fout, "%g\n", Props[p][i]);
		}

		fclose(fout);

		printf("Active cells fraction: %f, ACTNUM: %sfound, PORO: %sfound\n",
				double(count_act)/prop_size, (act_ind == -1 ? "not " : ""), (poro_ind == -1 ? "not " : ""));
	}
	catch (const Exception &e)
	{
		std::cerr << "ERROR: " << e.what() << "\n";
		return 1;
	}
	return 0;
}
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------
long StoL(std::string s, bool &complete)
{
	const char *input = s.c_str();	// 6.10.2013, using C++98 conversion
	char *end;
	double res = strtol(input, &end, 10);
	if (end == input || *end != '\0')
		complete = false;
	else
		complete = true;

	return res;
}
//------------------------------------------------------------------------------------------
double StoD(std::string s, bool &complete)
{
	const char *input = s.c_str();	// 6.10.2013, using C++98 conversion
	char *end;
	double res = strtod(input, &end);
	if (end == input || *end != '\0')
		complete = false;
	else
		complete = true;

	return res;
}
//------------------------------------------------------------------------------------------
long StoL(std::string s)
{
	bool complete;
	double res = StoL(s, complete);
	if (!complete)
		throw Exception("Cannot convert string '" + s + "' to long int");

	return res;
}
//------------------------------------------------------------------------------------------
double StoD(std::string s)
{
	bool complete;
	double res = StoD(s, complete);
	if (!complete)
		throw Exception("Cannot convert string '" + s + "' to double");

	return res;
}
//------------------------------------------------------------------------------------------
void ReadGrids(const char *file, size_t len, std::vector<std::vector<double>> &data, std::vector<std::string> S1, std::string S2)	// reads a number of grids from file
{															// allocates and fills "data" of size S1.size(), with data[i].size() = len
	const int Buff = 4096;									// S1[i], S2 - start and end of grid[i] which is loaded to data[i]
	char strmsg[Buff], str0[Buff];
	char *str = 0;
	FILE *File = fopen(file, "r");

	data = std::vector<std::vector<double>>(S1.size(), std::vector<double>(len));

	size_t GridCount = 0;	// counts total grids already read
	size_t ValCount = 0;	// counts values read for the current grid
	size_t c = 0;			// index within currently read grid
	bool seek_beg = true;
	bool new_line = true;
	int ind = -1;			// index within the grid names array

	if (File == 0)
		throw Exception((std::string)"Cannot open " + file + "\n");

	bool expect_scan_two = true;
	while (ReadTokenComm(File, &str, new_line, str0, Buff))	// reads a token to str, ignoring comments
	{
		std::string S = str;
		if (seek_beg)
		{
			ind = StrIndex(S, S1);
			if (ind != -1)					// found the starting string for grid-i
			{
				seek_beg = false;
				ValCount = 0;
				c = 0;
			}
			continue;
		}
		if (!seek_beg && S == S2)			// found the ending string
		{
			seek_beg = true;
			if (ValCount < len)
			{
				fclose(File);
				File = 0;
				sprintf(strmsg, " grid contains less values (%zu) than NI*NJ*NK (%zu)\n", ValCount, len);
				throw Exception(std::string(file) + ": " + S1[ind] + std::string(strmsg));
			}

			GridCount++;
			if (GridCount >= S1.size())		// all grids have been read
				break;

			continue;
		}
		if (!seek_beg)						// reading the main data
		{
			size_t cnt;
			double d;
			bool err = false;

			if (expect_scan_two)			// scan_two and scan_one are invoked based on the expectations (based on the previous successful scan)
			{
				if (!scan_two(str, cnt, d, expect_scan_two) && !scan_one(str, d, expect_scan_two))
					err = true;
			}
			else
			{
				if (!scan_one(str, d, expect_scan_two) && !scan_two(str, cnt, d, expect_scan_two))
					err = true;
			}

			if (!err)
			{
				if (!expect_scan_two)
					cnt = 1;

				ValCount += cnt;
				if (ValCount > len)			// too many values encountered
				{
					fclose(File);
					File = 0;
					sprintf(strmsg, " grid contains more values than NI*NJ*NK = %zu\n", len);
					throw Exception(std::string(file) + ": " + S1[ind] + std::string(strmsg));
				}

				if (expect_scan_two)
					for (size_t i = 0; i < cnt; i++)
					{
						data[ind][c] = d;
						c++;
					}
				else
				{
					data[ind][c] = d;
					c++;
				}
			}
			else							// error reading values
			{
				fclose(File);
				File = 0;
				sprintf(strmsg, " grid contains non-numeric symbol %s\n", str);
				throw Exception(std::string(file) + ": " + S1[ind] + std::string(strmsg));
			}
		}
	}
	fclose(File);
	File = 0;

	if (!seek_beg)
	{
		if (ValCount < len)
		{
			assert(ind != -1);
			sprintf(strmsg, " grid contains less values (%zu) than NI*NJ*NK (%zu)\n", ValCount, len);
			throw Exception(std::string(file) + ": " + S1[ind] + std::string(strmsg));
		}
		GridCount++;
	}

	if (GridCount < S1.size())
	{
		sprintf(strmsg, "Only %zu grid(s) found out of %zu\n", GridCount, S1.size());
		throw Exception(std::string(file) + ": " + std::string(strmsg));
	}
}
//------------------------------------------------------------------------------------------
void ReadGrid(const char *file, size_t len, std::vector<double> &data, std::string S1, std::string S2)		// same as above; reads one grid
{
	const int Buff = 4096;
	char strmsg[Buff], str0[Buff];
	char *str = 0;
	FILE *File = fopen(file, "r");

	data = std::vector<double>(len);

	size_t GridCount = 0;	// counts total grids already read
	size_t ValCount = 0;	// counts values read for the current grid
	size_t c = 0;			// index within currently read grid
	bool seek_beg = true;
	bool new_line = true;

	if (File == 0)
		throw Exception((std::string)"Cannot open " + file + "\n");

	bool expect_scan_two = true;
	while (ReadTokenComm(File, &str, new_line, str0, Buff))	// reads a token to str, ignoring comments
	{
		std::string S = str;
		if (seek_beg)
		{
			if (S == S1)					// found the starting string
			{
				seek_beg = false;
				ValCount = 0;
				c = 0;
			}
			continue;
		}
		if (!seek_beg && S == S2)			// found the ending string
		{
			seek_beg = true;
			if (ValCount < len)
			{
				fclose(File);
				File = 0;
				sprintf(strmsg, " grid contains less values (%zu) than expected (%zu)\n", ValCount, len);
				throw Exception(std::string(file) + ": " + S1 + std::string(strmsg));
			}

			GridCount++;					// the grid has been read
			break;
		}
		if (!seek_beg)						// reading the main data
		{
			size_t cnt;
			double d;
			bool err = false;

			if (expect_scan_two)			// scan_two and scan_one are invoked based on the expectations (based on the previous successful scan)
			{
				if (!scan_two(str, cnt, d, expect_scan_two) && !scan_one(str, d, expect_scan_two))
					err = true;
			}
			else
			{
				if (!scan_one(str, d, expect_scan_two) && !scan_two(str, cnt, d, expect_scan_two))
					err = true;
			}

			if (!err)
			{
				if (!expect_scan_two)
					cnt = 1;

				ValCount += cnt;
				if (ValCount > len)			// too many values encountered
				{
					fclose(File);
					File = 0;
					sprintf(strmsg, " grid contains more values than expected (%zu)\n", len);
					throw Exception(std::string(file) + ": " + S1 + std::string(strmsg));
				}

				if (expect_scan_two)
					for (size_t i = 0; i < cnt; i++)
					{
						data[c] = d;
						c++;
					}
				else
				{
					data[c] = d;
					c++;
				}
			}
			else							// error reading values
			{
				fclose(File);
				File = 0;
				sprintf(strmsg, " grid contains non-numeric symbol %s\n", str);
				throw Exception(std::string(file) + ": " + S1 + std::string(strmsg));
			}
		}
	}
	fclose(File);
	File = 0;

	if (!seek_beg)
	{
		if (ValCount < len)
		{
			sprintf(strmsg, " grid contains less values (%zu) than expected (%zu)\n", ValCount, len);
			throw Exception(std::string(file) + ": " + S1 + std::string(strmsg));
		}
		assert(ValCount == len);
		GridCount++;
	}

	if (GridCount < 1)
	{
		sprintf(strmsg, "Grid %.500s not found\n", S1.c_str());
		throw Exception(std::string(file) + ": " + std::string(strmsg));
	}
}
//------------------------------------------------------------------------------------------
bool ReadTokenComm(FILE *F, char **str, bool &new_line, char *str0, const int str0_len)
{																	// reads a token from the file (delimited by ' ', '\t', '\r', '\n'), dropping "--..." comments
	*str = 0;														// returns true on success, false on failure/EOF
																	// the token is saved to "str"
	static const char COMM[] = "--";		// comment beginning	// set new_line = true in the first call, then the function will manage it
	static const char DELIM[] = " \t\r\n";	// delimiters			// str0 is a working array (stores a line), it should have been allocated

	while (*str == 0)
	{
		if (new_line)
		{
			if (fgets(str0, str0_len, F) != 0)	// read the line
			{
				// remove the comment
				char *comm_ind = strstr(str0, COMM);
				if (comm_ind != 0)		// comment found
					comm_ind[0] = 0;	// set end-of-line at the comment start

				new_line = false;

				// get the first token
				*str = strtok(str0, DELIM);
			}
			else
				return false;
		}
		else
			*str = strtok(0, DELIM);

		if (*str == 0)
			new_line = true;
	}

	return true;
}
//------------------------------------------------------------------------------------------
int StrIndex(const std::string &S, const std::vector<std::string> &VECS)	// index of S in VECS[], -1 if not found
{
	for (size_t i = 0; i < VECS.size(); i++)
		//if (S.find(VECS[i]) != std::string::npos)		-- older version: VECS[i] can be only a substring in S
		if (S == VECS[i])
			return i;

	return -1;
}
//------------------------------------------------------------------------------------------
inline bool scan_two(const char *str, size_t &cnt, double &d, bool &expect_scan_two)	// parses "cnt*d", returns 'true' on success, updates 'expect_scan_two'
{
	char swork[8];
	swork[0] = '\0';

	int read = sscanf(str, "%zu*%lg%5s", &cnt, &d, swork);
	if (read == 2 && swork[0] == '\0')				// RPT*VAL successfully read
	{
		expect_scan_two = true;
		return true;
	}
	else
		return false;
}
//------------------------------------------------------------------------------------------
inline bool scan_one(const char *str, double &d, bool &expect_scan_two)					// parses "d", returns 'true' on success, updates 'expect_scan_two'
{
	char swork[8];
	swork[0] = '\0';

	int read = sscanf(str, "%lg%5s", &d, swork);
	if (read == 1 && swork[0] == '\0') 				// VAL successfully read
	{
		expect_scan_two = false;
		return true;
	}
	else
		return false;
}
//------------------------------------------------------------------------------------------
// indicates if cell cell_ind = "i, j, k" is active based on PORO > 0 and ACTIND > 0
bool is_active(size_t cell_ind, int poro_ind, int act_ind, const std::vector<std::vector<double>> &props)
{
	bool res = true;
	if (poro_ind != -1 && props[poro_ind][cell_ind] == 0)
		res = false;
	else if (act_ind != -1 && props[act_ind][cell_ind] == 0)
		res = false;

	return res;
}
//------------------------------------------------------------------------------------------




/*
C++ executable which can be called from MATLAB to speed up the conversion of DLD files to (t,x,y) format, 
and the computation of correlation functions from this data, using settings defined by a MATLAB front panel.

MATLAB: The relevant Matlab front panel can be found in
C:\Users\Delay.AMPLPC21\Documents\MATLAB\DLD_front_panel_c
To run it, open the m-file and click run (F5) to start the front panel

ARMADILLO: To use this linear algebra library, download from http://arma.sourceforge.net/
To install, copy include folder from download into the folder for VC++ 2012 (similarly for VC++ 2010 etc)
C:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\include

NOTE: ppl.h is compatible with the Armadillo, but amp.h isn't.
amp is for GPU stuff anyway which isnt appropriate here.
*/

#include "convert_and_correlate_header.h"
#include <stdexcept>
using namespace std;

struct config_vals_struct
{
	char filepath[500];
	int firstfile;
	int lastfile;
	int num_CPU;
	int num_GPU;
	int convert_DLD_txy;
	double t_min;
	double t_max;
	double x_min;
	double x_max;
	double y_min;
	double y_max;
	int axis;
	int order;
	double t_bin;
	double x_bin;
	double y_bin;
	int num_bins;
	int tof_plot;
	int spatial_plot;
	int nvf;
	int fit_remove_t0;
	int fit_temp;
	int files_per_chunk;
	double temp_window_pc;
	double num_window_pc;
	double t_bin_plot;
	int bins_per_pixel;
	int dead_time; 
	int exclusion_active;
	double t_min_exc;
	double t_max_exc;
	double x_min_exc;
	double x_max_exc;
	double y_min_exc;
	double y_max_exc;
	char filepath_browse[500];
	char filepath_reload[500];
	int TF_fit;
	double f_x;
	double f_y;
	double f_z;
	double temp_guess_TF;
	double TF_radius_guess;
	double cond_percent_guess;
	int force_recalc_txy;
	int ready_for_gn;
	int num_files_in_denom;
	int group_denom_temp;
};

struct channel_triggers_struct
{
	double ch0;
	double ch1;
	double ch2;
	double ch3;
	double rubbish;
};

int read_config(struct config_vals_struct *p_config_vals)
{
	ifstream config_file("config123.txt");
	if(!config_file)
	{
		ifstream config_file("C:\\Users\\Andrew\\Documents\\Visual Studio 2010\\Projects\\DLD_convert_and_correlations_for_MATLAB\\config123.txt");
		if(!config_file)
		{
			printf("ERROR: cannot find config file in current folder \n");
			return 1;
		}
	}

	for (int line = 0; line < 49; line++)
	{
		char buffer[500];
		config_file.getline(buffer,500,';');
		config_file.ignore(500,'\n');

		switch(line)
		{
		case 0:
			strncpy(p_config_vals->filepath, buffer, 500);
			cout << buffer << endl;
			break;
		case 1:
			p_config_vals->firstfile = atoi(buffer);
			break;
		case 2:
			p_config_vals->lastfile = atoi(buffer);
			break;
		case 3:
			p_config_vals->num_CPU = atoi(buffer);
			break;
		case 4:
			p_config_vals->num_GPU = atoi(buffer);
			break;
		case 5:
			p_config_vals->convert_DLD_txy = atoi(buffer);
			break;
		case 6:
			p_config_vals->t_min = atof(buffer);
			break;
		case 7:
			p_config_vals->t_max = atof(buffer);
			break;
		case 8:
			p_config_vals->x_min = atof(buffer)/1000.0;
			break;
		case 9:
			p_config_vals->x_max = atof(buffer)/1000.0;
			break;
		case 10:
			p_config_vals->y_min = atof(buffer)/1000.0;
			break;
		case 11:
			p_config_vals->y_max = atof(buffer)/1000.0;
			break;
		case 12:
			p_config_vals->axis = atoi(buffer);
			break;
		case 13:
			p_config_vals->order = atoi(buffer);
			break;
		case 14:
			p_config_vals->t_bin = atof(buffer)/1000.0;
			break;
		case 15:
			p_config_vals->x_bin = atof(buffer)/1000.0;
			break;
		case 16:
			p_config_vals->y_bin = atof(buffer)/1000.0;
			break;
		case 17:
			p_config_vals->num_bins = atoi(buffer);
			break;
		case 18:
			p_config_vals->tof_plot = atoi(buffer);
			break;
		case 19:
			p_config_vals->spatial_plot = atoi(buffer);
			break;
		case 20:
			p_config_vals->nvf = atoi(buffer);
			break;
		case 21:
			p_config_vals->fit_remove_t0 = atoi(buffer);
			break;
		case 22:
			p_config_vals->fit_temp = atoi(buffer);
			break;
		case 23:
			p_config_vals->files_per_chunk = atoi(buffer);
			break;
		case 24:
			p_config_vals->temp_window_pc = atof(buffer);
			break;
		case 25:
			p_config_vals->num_window_pc = atof(buffer);
			break;
		case 26:
			p_config_vals->t_bin_plot = atof(buffer)/1000.0;
			break;
		case 27:
			p_config_vals->bins_per_pixel = atoi(buffer);
			break;
		case 28:
			p_config_vals->dead_time = atoi(buffer);
			break;
		case 29:
			p_config_vals->exclusion_active = atoi(buffer);
			break;
		case 30:
			p_config_vals->t_min_exc = atof(buffer);
			break;
		case 31:
			p_config_vals->t_max_exc = atof(buffer);
			break;
		case 32:
			p_config_vals->x_min_exc = atof(buffer);
			break;
		case 33:
			p_config_vals->x_max_exc = atof(buffer);
			break;
		case 34:
			p_config_vals->y_min_exc = atof(buffer);
			break;
		case 35:
			p_config_vals->y_max_exc = atof(buffer);
			break;
		case 36:
			strncpy(p_config_vals->filepath_browse, buffer, 500);
			break;
		case 37:
			strncpy(p_config_vals->filepath_reload, buffer, 500);
			break;
		case 38:
			p_config_vals->TF_fit = atoi(buffer);
			break;
		case 39:
			p_config_vals->f_x = atof(buffer);
			break;
		case 40:
			p_config_vals->f_y = atof(buffer);
			break;
		case 41:
			p_config_vals->f_z = atof(buffer);
			break;
		case 42:
			p_config_vals->temp_guess_TF = atof(buffer);
			break;
		case 43:
			p_config_vals->TF_radius_guess = atof(buffer);
			break;
		case 44:
			p_config_vals->cond_percent_guess = atof(buffer);
			break;
		case 45:
			p_config_vals->force_recalc_txy = atoi(buffer);
			break;
		case 46:
			p_config_vals->ready_for_gn = atoi(buffer);
			break;
		case 47:
			p_config_vals->num_files_in_denom = atoi(buffer);
			break;
		case 48:
			p_config_vals->group_denom_temp = atoi(buffer);
			break;
		}
	}

	return 0;
}

int DLD_txy(struct config_vals_struct *p_config_vals, struct channel_triggers_struct *p_channel_triggers)
{
	/*
	Meant to be reproduction of Sean's dld_read_5channels_reconst_multi_a.m 
	Took out some superfluous stuff, always normalises time, only looks for 4 channel triggers
	Added in the 0.61 radian rotation.
	*/

	int size_nvf_vec = p_config_vals->lastfile - p_config_vals->firstfile + 1;	//REMOVE
	long* nvf = new long[size_nvf_vec];
	long* p_nvf = &nvf[0];

	char filepath_temp [200];
	char buffer [14];	// reading 12 digit numbers

	long size_TOF_vec = (long) ((p_config_vals->t_max - p_config_vals->t_min)*100000);	//limit to 10 microsec accuracy - saturates at 1MHz.

	long* TOF = new long[size_TOF_vec];
	long* p_TOF = &TOF[0];

	long* spatial_image = new long[368449];	// resolution limit is 304 bins = 8cm/132micron, 304^2 bins
	long* p_spatial_image = &spatial_image[0];


	if (p_config_vals->tof_plot == 1)
	{
		for (long i = 0; i<size_TOF_vec; i++)
		{
			TOF[i] = 0;
		}
	}

	if (p_config_vals->spatial_plot == 1)
	{
		for(int i = 0; i<368449; i++)
		{
			spatial_image[i] = 0;
		}
	}

	for (int file = p_config_vals->firstfile; file <= p_config_vals->lastfile; file++)
	{
		strncpy(filepath_temp,p_config_vals->filepath, 200);
		strcat(filepath_temp,itoa(file,buffer,10));
		strcat(filepath_temp,".txt");
		//cout << filepath_temp << '\n';

		ifstream file_to_read(filepath_temp);
		if (!file_to_read)
		{
			printf("ERROR: cannot find DLD raw file in current folder \n");
			cout << filepath_temp << '\n';
			file_to_read.close();
			//return 1;
			continue;
		}

		//determine number of rows in file to allocate arrays
		long begin,end,numrows;
		begin = file_to_read.tellg();
		file_to_read.seekg (0, ios::end);
		end = file_to_read.tellg();
		numrows = (end-begin)/16;
		file_to_read.seekg (0, ios::beg);	// send back to beginning for read

		//cout << "Number of rows to read is " << numrows << '\n';

		char temp_read[20];	// restrict to 9 digits, even though full data has 10, as long is 2^32 = 4 billion

		// if we avoid using Armadillo
		int* channel = new int[numrows];
		int* p_channel = &channel[0];
		long long* timestamp = new long long[numrows];
		long long* p_timestamp = &timestamp[0];

		//umat DLD_raw_data = umat(numrows, 2);	//spare Arma version

		// using Armadillo - DOG SLOW!!! DO NOT USE!!!

		/*mat DLD_raw_data = mat(numrows, 2);
		seconds1 = time(NULL);
		mat DLD_raw_data; DLD_raw_data.load(file_to_read, csv_ascii);
		seconds2 = time(NULL);
		cout << "Armadillo load seconds: " << seconds2 - seconds1 << endl;
		cout << DLD_raw_data(0,1) << endl;
		file_to_read.seekg (0, ios::beg);
		*/

		const int digits_to_read = 12;					// can set to 9 if you want 32 bit numbers

		long row_index = 0;
		//char temp_timestamp_long[digits_to_read] = "";
		long long timestamp_long = 0;
		long channel_number = 0;
		long long timestamp_0 = 0;
		long write_row = 0;
		int rubbish = 0;							// many files have a few dozen garbage triggers directly after initial pulse
		long ch0_long = 0;							// temporarily hold number of counts per trigger. convert to double later
		long ch1_long = 0;							// if left as integer, gets really big, but double wont ++ properly
		long ch2_long = 0;
		long ch3_long = 0;

		// not used!!!
		double timestamp_0_d = 0;					// using double precision to get all 12 digits (long is 9)
		double timestamp_d = 0;


		for(long row_index = 0; row_index < numrows*10; row_index++)
		{
			if (!file_to_read)
			{
				if (file_to_read.eof())
				{
					//cout << "end of file" << '\n';
					break;
				}
				else
				{
					cout << "Fail: " << file_to_read.fail() << '\n';
					cout << "Bad: " << file_to_read.bad() << '\n';
					printf("ERROR: DLD raw file failed bit \n");
					file_to_read.close();
					//return 1;
					continue;
				}
			}

			file_to_read.getline(temp_read,4,',');
			channel_number = atol(temp_read);
			file_to_read.getline(temp_read,14,'\n');
			timestamp_d = atof(temp_read);
			timestamp_long = (long long) timestamp_d;

			if (row_index == 0)
			{
				if (channel_number == 5)
				{
					timestamp_0 = timestamp_long;
					continue;
				}
				else
				{
					printf("ERROR: DLD raw trigger timestamp is wrong \n");
					file_to_read.close();
					return 1;				
				}
			}

			timestamp_long = timestamp_long - timestamp_0;

			if (write_row > 0 && timestamp_long < 0)
			{
				rubbish++;
				continue;
			}

			if (timestamp_long >= 0 && channel_number >= 0  && channel_number < 4 )
			{
				switch(channel_number)
				{
				case 0:
					ch0_long++;
					break;
				case 1:
					ch1_long++;
					break;
				case 2:
					ch2_long++;
					break;
				case 3:
					ch3_long++;
					break;
				}

				channel[write_row] = channel_number;
				timestamp[write_row] = timestamp_long;
				//cout << "write line" << endl;

				//DLD_raw_data(write_row,0) = channel_number;
				//DLD_raw_data(write_row,1) = timestamp_long;

				write_row++;
			}
			else if (isnan_agm(timestamp_long) != 0 || isnan_agm(channel_number) != 0)
			{
				//cout << "got here" << endl;
				break;
			}
			/*else
			{
			cout << "skipped a line" << endl;
			}*/
			if (write_row >= numrows)
			{
				//cout << "exit loop" << endl;
				break;
			}
		}
		//cout << "Channel: " << channel[write_row-1] << ", Timestamp: "<< temp_read << ", Timestamp adjust: "<< timestamp[write_row-1] + timestamp_0 << '\n';



		p_channel_triggers->ch0 = p_channel_triggers->ch0 + (double)ch0_long;
		p_channel_triggers->ch1 = p_channel_triggers->ch1 + (double)ch1_long;
		p_channel_triggers->ch2 = p_channel_triggers->ch2 + (double)ch2_long;
		p_channel_triggers->ch3 = p_channel_triggers->ch3 + (double)ch3_long;
		p_channel_triggers->rubbish = p_channel_triggers->rubbish + (double)rubbish;

		// can sort, but data already appears sorted as is.
		/*umat DLD_raw_data_sorted = sort(DLD_raw_data,0,1);
		cout << DLD_raw_data_sorted.rows(write_row-5,write_row-1) << endl;
		cout << "cf." << endl;
		cout << DLD_raw_data.rows(write_row-5,write_row-1) << endl;*/

		long number_detections = write_row;										// Will equal [5*n + 1 2] for n detections ideally, the +1 is a master trigger to throw out
		int number_successes = 0;                                               // Tally successful hits
		long which_row = 0;                                                     // Index of matrix row to write to
		long tolerance_keep = 200;                              // Tolerance in what data to keep - should be 200 bins of 25ps
		int search_no = 36;                                                     // Seach over 9 (=36/4) complete hits
		int T_spread = 0;                                                       // Spreads in times reconstructed
		int T_sum_spread = 0;
		int T_sum_x = 0;
		int T_sum_y = 0;

		int T_sum_tol_keep = T_sum + tolerance_keep;

		//umat txy_output = umat(long(write_row/4), 5);

		// should use these eventually, get rid of other mantrices
		/*long* three_channel_out = new long[long(3*write_row/4)];
		long* p_three_channel_out = &three_channel_out[0];

		umat three_channel_out_arma = umat(long(write_row/4), 3);

		for (int i = 0; i < long(3*write_row/4); i++)
		{
		three_channel_out[i] = -1;
		}*/


		long dummy = 0;
		long dummy2 = 0;
		long dummy3 = 0;

		long long x1 = 0;
		long long x2 = 0;
		long long y1 = 0;
		long long y2 = 0;
		long how_many_x = 0;
		long how_many_y = 0;
		double tx = 0;
		double t = 0;
		long x2_index = 0;

		double x_pre_rot;
		double y_pre_rot;
		double x;
		double y;

		FILE *out_txy;
		char outfile_txy_handle[200] = "";
		strncpy(outfile_txy_handle,p_config_vals->filepath, 150);
		strcat(outfile_txy_handle,"_txy_forc_AGM_");
		strcat(outfile_txy_handle,itoa(file,buffer,10));
		strcat(outfile_txy_handle,".txt");
		out_txy = fopen(outfile_txy_handle,"w");

		for (long count = 1; count < write_row; count++)
		{
			if (channel[count] == 0)
			{
				dummy = 0;
				how_many_x = 0;

				//cout << "x1 " << timestamp[count] + timestamp_0 << endl;
				//system("pause");

				for (long count2 = 1; count2 <= search_no; count2++)
				{
					dummy = dummy+count2*powl(-1,count2);

					if ((count + dummy) <= 0 || (count + dummy) > write_row)
					{
						continue;
					}
					if (abs(timestamp[count + dummy] - timestamp[count]) > T_sum)	//FIX
					{
						continue;
					}

					if (channel[count + dummy] == 1)
					{
						if (timestamp[count] == timestamp[count + dummy])
						{
							continue;
						}

						if (abs(timestamp[count] - timestamp[count + dummy]) < T_sum_tol_keep)
						{
							how_many_x++;

							if (how_many_x > 1)
							{
								break;
							}
							else
							{
								x1 = timestamp[count];
								x2 = timestamp[count + dummy];
								which_row++;
								tx = 0.5*(x1+x2-T_sum);

								//cout << count << ", " << 1.0*tx/4294967296 << endl;

								x2_index = count + dummy;
							}
						}
					}
				}

				if (how_many_x != 1)
				{
					continue;
				}

				dummy2 = 0;
				how_many_y = 0;

				for (int count3 = 1; count3 <= search_no; count3++)
				{
					dummy2 = dummy2+count3*powl(-1,count3);

					if ((count + dummy2) <= 0 || (count + dummy2) > write_row)
					{
						continue;
					}
					if (abs(timestamp[count + dummy2] - timestamp[count]) > T_sum)
					{
						continue;
					}

					if (channel[count + dummy2] == 2)
					{
						dummy3 = 0;

						for (int count4 = 1; count4 <= search_no; count4++)
						{
							dummy3 = dummy3+count4*powl(-1,count4);

							if ((count + dummy3) <= 0 || (count + dummy3) > write_row)
							{
								continue;
							}
							if (abs(timestamp[count + dummy3] - timestamp[count]) > T_sum)
							{
								continue;
							}

							if (channel[count + dummy3] == 3)
							{
								if (timestamp[count + dummy2] == timestamp[count + dummy3])
								{
									continue;
								}

								if (abs(timestamp[count + dummy2] - timestamp[count + dummy3]) < T_sum_tol_keep && abs(x1+x2-timestamp[count + dummy2]-timestamp[count + dummy3]) < 2*tolerance_keep)
								{

									how_many_y++;
									if (how_many_y > 1)
									{
										break;
									}
									else
									{
										y1 = timestamp[count + dummy2];
										y2 = timestamp[count + dummy3];

										t = (tx + (0.5*(y1+y2-T_sum)))/2;

										channel[count] = -1;   //now that we've read the data, NaN the entries in dld_output_sorted we've read for later removal
										channel[x2_index] = -1;
										channel[count+dummy2] = -1;
										channel[count+dummy3] = -1;

										number_successes++;

										//cout << x1 << ", " << x2 << ", " << y1 << ", " << y2 << endl;
										//system("pause");

										x_pre_rot = (x1 - x2)*v_perp_x*bin_time;
										y_pre_rot = (y1 - y2)*v_perp_y*bin_time;
										x = x_pre_rot*cos(0.61) - y_pre_rot*sin(0.61);
										y = x_pre_rot*sin(0.61) + y_pre_rot*cos(0.61);
										t=abs(t*bin_time);

										// put hits in TOF
										if (t >= p_config_vals->t_min && t < p_config_vals->t_max && x >= p_config_vals->x_min && x <= p_config_vals->x_max && y >= p_config_vals->y_min && y <= p_config_vals->y_max)
										{
											if (p_config_vals->tof_plot == 1)
											{
												TOF[long((t - p_config_vals->t_min)*100000)]++;
											}
											if (p_config_vals->spatial_plot == 1)
											{
												if ((long((x - p_config_vals->x_min)/0.00013158*607) + long((y - p_config_vals->y_min)/0.00013158)) < 368449)
												{
													//cout << 368449 - (long((x - p_config_vals->x_min)/0.00013158*607) + long((y - p_config_vals->y_min)/0.00013158)) << endl;
													spatial_image[long((x - p_config_vals->x_min)/0.000132)*607 + long((y - p_config_vals->y_min)/0.000132)]++;
												}
											}
										}

										fprintf(out_txy,"%.9f,%.9f,%.9f\n",t,x,y);

										// put hit in TOF and spatial image

									}
								}
							}
						}
					}
				}
				if (how_many_y == 0)
				{
					which_row--;
				}
			}
		}
		cout << "File number " << file << " with this many hits "<< number_successes << endl;
		nvf[file-1] = number_successes;

		fclose(out_txy);
		file_to_read.close();
		delete [] channel;
		delete [] timestamp;
	}

	FILE *out_TOF;
	char outfile_TOF_handle[200] = "";
	strncpy(outfile_TOF_handle,p_config_vals->filepath, 150);
	strcat(outfile_TOF_handle,"_TOF.txt");
	out_TOF = fopen(outfile_TOF_handle,"wt");
	for (long i = 0; i < size_TOF_vec; i++)
	{
		fprintf(out_TOF,"%lu; \n",TOF[i]);
	}
	fclose(out_TOF);
	delete [] TOF;

	//cout << "TOF done" << endl;

	FILE *out_spatial;
	char outfile_spatial_handle[200] = "";
	strncpy(outfile_spatial_handle,p_config_vals->filepath, 150);
	strcat(outfile_spatial_handle,"_spatial.txt");
	out_spatial = fopen(outfile_spatial_handle,"wt");
	//cout << "about to do space" << endl;
	for (long i = 0; i < 368449; i++)
	{
		fprintf(out_spatial,"%lu,",spatial_image[i]);
	}
	fclose(out_spatial);
	delete [] spatial_image;

	FILE *out_nvf;
	char outfile_nvf_handle[200] = "";
	strncpy(outfile_nvf_handle,p_config_vals->filepath, 150);
	strcat(outfile_nvf_handle,"_nvf.txt");
	out_nvf = fopen(outfile_nvf_handle,"wt");
	for (int i = 0; i < size_nvf_vec; i++)
	{
		fprintf(out_nvf,"%lu; \n",nvf[i]);
	}
	fclose(out_nvf);
	delete [] nvf;

	return 0;
}

int crunch_gn_single_core(struct config_vals_struct *p_config_vals, int hits_this_file, 
						  int hist_axis, int bin_axis_1, int bin_axis_2,
						  double hist_binsize, double bin_1_binsize, double bin_2_binsize,
						  double* p_txy_this_file,
						  uword* p_g2_0, uword* p_g2_1, uword* p_g2_2, uword* p_g2_3, uword* p_g2_m1, uword* p_g2_m2, uword* p_g2_m3,
						  uword* p_g3_0, uword* p_g3_1, uword* p_g3_2, uword* p_g3_3, uword* p_g3_m1, uword* p_g3_m2, uword* p_g3_m3)
{
	double dh_12, d1_12, d2_12, dh_13, d1_13, d2_13, dh_23, d1_23, d2_23;
	int hist_bin_12, hist_bin_13, hist_bin_23;
	long number_pairs = 0, number_triples = 0;	// count these at the original bin size

	for (int hit1 = 0; hit1 < hits_this_file - 2; hit1++)
	{
		for (int hit2 = hit1+1; hit2 < hits_this_file - 1; hit2++)
		{
			dh_12 = *(p_txy_this_file+hit2+hist_axis*hits_this_file) - *(p_txy_this_file+hit1+hist_axis*hits_this_file);
			d1_12 = *(p_txy_this_file+hit2+bin_axis_1*hits_this_file) - *(p_txy_this_file+hit1+bin_axis_1*hits_this_file);
			d2_12 = *(p_txy_this_file+hit2+bin_axis_2*hits_this_file) - *(p_txy_this_file+hit1+bin_axis_2*hits_this_file);

			if (hist_axis == 0)	// time
			{
				if (dh_12 > p_config_vals->num_bins * hist_binsize)	// hit off end of histogram
				{
					break;
				}
				if (dh_12 < 0) // prevent double counting
				{
					continue;
				}
				if (abs(d1_12) > bin_1_binsize * 1.6 || abs(d2_12) > bin_2_binsize * 1.6) // too far apart for biggest bin
				{
					continue;
				}
			}

			if (hist_axis == 1 || hist_axis == 2)	// x or y
			{
				if (d1_12 > bin_1_binsize * 1.6)	// will be time
				{
					break;
				}
				if (d1_12 < 0)		// prevent double counting
				{
					continue;
				}
				if (abs(d2_12) > bin_2_binsize * 1.6)
				{
					continue;
				}
				if (abs(dh_12) > p_config_vals->num_bins * hist_binsize)
				{
					continue;
				}
			}

			hist_bin_12 = ceil(abs(dh_12) / hist_binsize) - 1;

			if (hist_bin_12 < p_config_vals->num_bins && hist_bin_12 >= 0)
			{
				if (d1_12 < bin_1_binsize * 1.6 && d2_12 < bin_2_binsize * 1.6)
				{
					*(p_g2_3 + hist_bin_12) += 1;
					if (d1_12 < bin_1_binsize * 1.4 && d2_12 < bin_2_binsize * 1.4)
					{
						*(p_g2_2 + hist_bin_12) += 1;
						if (d1_12 < bin_1_binsize * 1.2 && d2_12 < bin_2_binsize * 1.2)
						{
							*(p_g2_1 + hist_bin_12) += 1;
							if (d1_12 < bin_1_binsize && d2_12 < bin_2_binsize)
							{		
								*(p_g2_0 + hist_bin_12) += 1;
								number_pairs++;
								if (d1_12 < bin_1_binsize * 0.8 && d2_12 < bin_2_binsize * 0.8)
								{
									*(p_g2_m1 + hist_bin_12) += 1;
									if (d1_12 < bin_1_binsize * 0.6 && d2_12 < bin_2_binsize * 0.6)
									{
										*(p_g2_m2 + hist_bin_12) += 1;
										if (d1_12 < bin_1_binsize * 0.4 && d2_12 < bin_2_binsize * 0.4)
										{
											*(p_g2_m3 + hist_bin_12) += 1;
										}
									}
								}
							}
						}
					}
				}

				if (p_config_vals->order > 2)
				{
					for (int hit3 = hit2+1; hit3 < hits_this_file; hit3++)
					{
						dh_23 = *(p_txy_this_file+hit3+hist_axis*hits_this_file) - *(p_txy_this_file+hit2+hist_axis*hits_this_file);
						d1_23 = *(p_txy_this_file+hit3+bin_axis_1*hits_this_file) - *(p_txy_this_file+hit2+bin_axis_1*hits_this_file);
						d2_23 = *(p_txy_this_file+hit3+bin_axis_2*hits_this_file) - *(p_txy_this_file+hit2+bin_axis_2*hits_this_file);

						if (hist_axis == 0)	// time
						{
							if ((dh_23 > p_config_vals->num_bins * hist_binsize))	// hit off end of histogram
							{
								break;
							}
							if (dh_23 < 0) // prevent double counting
							{
								continue;
							}
							if (abs(d1_23) > bin_1_binsize * 1.6 || abs(d2_23) > bin_2_binsize * 1.6) // too far apart for biggest bin
							{
								continue;
							}
						}

						if (hist_axis == 1 || hist_axis == 2)	// x or y
						{
							if (d1_23 > bin_1_binsize * 1.6)	// will be time
							{
								break;
							}
							if (d1_23 < 0)		// prevent double counting
							{
								continue;
							}
							if (abs(d2_23) > bin_2_binsize * 1.6)
							{
								continue;
							}
							if (abs(dh_23) > p_config_vals->num_bins * hist_binsize)
							{
								continue;
							}
						}

						hist_bin_23 = ceil(abs(dh_23) / hist_binsize) - 1;

						if (hist_bin_23 < p_config_vals->num_bins && hist_bin_23 >= 0)
						{
							if (d1_23 < bin_1_binsize * 1.6 && d2_23 < bin_2_binsize * 1.6)
							{
								*(p_g3_3 + hist_bin_12 + p_config_vals->num_bins*hist_bin_23) += 1;
								if (d1_23 < bin_1_binsize * 1.4 && d2_23 < bin_2_binsize * 1.4)
								{
									*(p_g3_2 + hist_bin_12 + p_config_vals->num_bins*hist_bin_23) += 1;
									if (d1_23 < bin_1_binsize * 1.2 && d2_23 < bin_2_binsize * 1.2)
									{
										*(p_g3_1 + hist_bin_12 + p_config_vals->num_bins*hist_bin_23) += 1;
										if (d1_23 < bin_1_binsize && d2_23 < bin_2_binsize)
										{		
											*(p_g3_0 + hist_bin_12 + p_config_vals->num_bins*hist_bin_23) += 1;
											number_triples++;
											if (d1_23 < bin_1_binsize * 0.8 && d2_23 < bin_2_binsize * 0.8)
											{
												*(p_g3_m1 + hist_bin_12 + p_config_vals->num_bins*hist_bin_23) += 1;
												if (d1_23 < bin_1_binsize * 0.6 && d2_23 < bin_2_binsize * 0.6)
												{
													*(p_g3_m2 + hist_bin_12 + p_config_vals->num_bins*hist_bin_23) += 1;
													if (d1_23 < bin_1_binsize * 0.4 && d2_23 < bin_2_binsize * 0.4)
													{
														*(p_g3_m3 + hist_bin_12 + p_config_vals->num_bins*hist_bin_23) += 1;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	//cout << "Number pairs single core: " << number_pairs << endl;
	//cout << "Number triplets single core: " << number_triples << endl;

	return 0;
}

int crunch_gn_multi_core(struct config_vals_struct *p_config_vals, int hits_this_file, 
						 int hist_axis, int bin_axis_1, int bin_axis_2,
						 double hist_binsize, double bin_1_binsize, double bin_2_binsize,
						 double* p_txy_this_file,
						 uword* p_g2_0, uword* p_g2_1, uword* p_g2_2, uword* p_g2_3, uword* p_g2_m1, uword* p_g2_m2, uword* p_g2_m3,
						 uword* p_g3_0, uword* p_g3_1, uword* p_g3_2, uword* p_g3_3, uword* p_g3_m1, uword* p_g3_m2, uword* p_g3_m3)
{
	double dh_12, d1_12, d2_12, dh_13, d1_13, d2_13, dh_23, d1_23, d2_23;
	int hist_bin_12, hist_bin_13, hist_bin_23;
	long number_pairs = 0, number_triples = 0;	// count these at the original bin size

	parallel_for (0, hits_this_file - 2, [hist_axis, bin_axis_1, bin_axis_2,
		hist_binsize, bin_1_binsize, bin_2_binsize, 
		p_g2_3, p_g2_2, p_g2_1, p_g2_0, p_g2_m1, p_g2_m2, p_g2_m3,
		p_g3_3, p_g3_2, p_g3_1, p_g3_0, p_g3_m1, p_g3_m2, p_g3_m3,
		hits_this_file, p_config_vals, p_txy_this_file] (int hit1)
	//for (int hit1 = 0; hit1 < hits_this_file - 2; hit1++)
	{
		double dh_12, d1_12, d2_12;
		double dh_23, d1_23, d2_23;
		int hist_bin_12, hist_bin_23;
		int number_pairs, number_triples; 

		for (int hit2 = hit1+1; hit2 < hits_this_file - 1; hit2++)
		{
			dh_12 = *(p_txy_this_file+hit2+hist_axis*hits_this_file) - *(p_txy_this_file+hit1+hist_axis*hits_this_file);
			d1_12 = *(p_txy_this_file+hit2+bin_axis_1*hits_this_file) - *(p_txy_this_file+hit1+bin_axis_1*hits_this_file);
			d2_12 = *(p_txy_this_file+hit2+bin_axis_2*hits_this_file) - *(p_txy_this_file+hit1+bin_axis_2*hits_this_file);

			if (hist_axis == 0)	// time
			{
				if (dh_12 > p_config_vals->num_bins * hist_binsize)	// hit off end of histogram
				{
					break;
				}
				if (dh_12 < 0) // prevent double counting
				{
					continue;
				}
				if (abs(d1_12) > bin_1_binsize * 1.6 || abs(d2_12) > bin_2_binsize * 1.6) // too far apart for biggest bin
				{
					continue;
				}
			}

			if (hist_axis == 1 || hist_axis == 2)	// x or y
			{
				if (d1_12 > bin_1_binsize * 1.6)	// will be time
				{
					break;
				}
				if (d1_12 < 0)		// prevent double counting
				{
					continue;
				}
				if (abs(d2_12) > bin_2_binsize * 1.6)
				{
					continue;
				}
				if (abs(dh_12) > p_config_vals->num_bins * hist_binsize)
				{
					continue;
				}
			}

			int hist_bin_12 = ceil(abs(dh_12) / hist_binsize) - 1;

			if (hist_bin_12 < p_config_vals->num_bins && hist_bin_12 >= 0)
			{
				if (d1_12 < bin_1_binsize * 1.6 && d2_12 < bin_2_binsize * 1.6)
				{
					*(p_g2_3 + hist_bin_12) += 1;
					if (d1_12 < bin_1_binsize * 1.4 && d2_12 < bin_2_binsize * 1.4)
					{
						*(p_g2_2 + hist_bin_12) += 1;
						if (d1_12 < bin_1_binsize * 1.2 && d2_12 < bin_2_binsize * 1.2)
						{
							*(p_g2_1 + hist_bin_12) += 1;
							if (d1_12 < bin_1_binsize && d2_12 < bin_2_binsize)
							{		
								*(p_g2_0 + hist_bin_12) += 1;
								number_pairs++;
								if (d1_12 < bin_1_binsize * 0.8 && d2_12 < bin_2_binsize * 0.8)
								{
									*(p_g2_m1 + hist_bin_12) += 1;
									if (d1_12 < bin_1_binsize * 0.6 && d2_12 < bin_2_binsize * 0.6)
									{
										*(p_g2_m2 + hist_bin_12) += 1;
										if (d1_12 < bin_1_binsize * 0.4 && d2_12 < bin_2_binsize * 0.4)
										{
											*(p_g2_m3 + hist_bin_12) += 1;
										}
									}
								}
							}
						}
					}
				}

				if (p_config_vals->order > 2)
				{
					for (int hit3 = hit2+1; hit3 < hits_this_file; hit3++)
					{
						dh_23 = *(p_txy_this_file+hit3+hist_axis*hits_this_file) - *(p_txy_this_file+hit2+hist_axis*hits_this_file);
						d1_23 = *(p_txy_this_file+hit3+bin_axis_1*hits_this_file) - *(p_txy_this_file+hit2+bin_axis_1*hits_this_file);
						d2_23 = *(p_txy_this_file+hit3+bin_axis_2*hits_this_file) - *(p_txy_this_file+hit2+bin_axis_2*hits_this_file);

						if (hist_axis == 0)	// time
						{
							if ((dh_23 > p_config_vals->num_bins * hist_binsize))	// hit off end of histogram
							{
								break;
							}
							if (dh_23 < 0) // prevent double counting
							{
								continue;
							}
							if (abs(d1_23) > bin_1_binsize * 1.6 || abs(d2_23) > bin_2_binsize * 1.6) // too far apart for biggest bin
							{
								continue;
							}
						}

						if (hist_axis == 1 || hist_axis == 2)	// x or y
						{
							if (d1_23 > bin_1_binsize * 1.6)	// will be time
							{
								break;
							}
							if (d1_23 < 0)		// prevent double counting
							{
								continue;
							}
							if (abs(d2_23) > bin_2_binsize * 1.6)
							{
								continue;
							}
							if (abs(dh_23) > p_config_vals->num_bins * hist_binsize)
							{
								continue;
							}
						}

						int hist_bin_23 = ceil(abs(dh_23) / hist_binsize) - 1;

						if (hist_bin_23 < p_config_vals->num_bins && hist_bin_23 >= 0)
						{
							if (d1_23 < bin_1_binsize * 1.6 && d2_23 < bin_2_binsize * 1.6)
							{
								*(p_g3_3 + hist_bin_12 + p_config_vals->num_bins*hist_bin_23) += 1;
								if (d1_23 < bin_1_binsize * 1.4 && d2_23 < bin_2_binsize * 1.4)
								{
									*(p_g3_2 + hist_bin_12 + p_config_vals->num_bins*hist_bin_23) += 1;
									if (d1_23 < bin_1_binsize * 1.2 && d2_23 < bin_2_binsize * 1.2)
									{
										*(p_g3_1 + hist_bin_12 + p_config_vals->num_bins*hist_bin_23) += 1;
										if (d1_23 < bin_1_binsize && d2_23 < bin_2_binsize)
										{		
											*(p_g3_0 + hist_bin_12 + p_config_vals->num_bins*hist_bin_23) += 1;
											number_triples++;
											if (d1_23 < bin_1_binsize * 0.8 && d2_23 < bin_2_binsize * 0.8)
											{
												*(p_g3_m1 + hist_bin_12 + p_config_vals->num_bins*hist_bin_23) += 1;
												if (d1_23 < bin_1_binsize * 0.6 && d2_23 < bin_2_binsize * 0.6)
												{
													*(p_g3_m2 + hist_bin_12 + p_config_vals->num_bins*hist_bin_23) += 1;
													if (d1_23 < bin_1_binsize * 0.4 && d2_23 < bin_2_binsize * 0.4)
													{
														*(p_g3_m3 + hist_bin_12 + p_config_vals->num_bins*hist_bin_23) += 1;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	);
	//cout << "Number pairs single core: " << number_pairs << endl;
	//cout << "Number triplets single core: " << number_triples << endl;

	return 0;
}

int calculate_gn(struct config_vals_struct *p_config_vals)
{
	// t_bin, x_bin, y_bin are in ms or mm
	// axis is 1=t, 2=x, 3=y
	// order is 0 (don't do g(n) - but you wont be here in this case), 2 = g(2), 3 = g(3)
	// num_bins can be in any of the three axes to plot g(n) along

	int hist_axis, bin_axis_1, bin_axis_2;
	double hist_binsize, bin_1_binsize, bin_2_binsize;

	switch (p_config_vals->axis)
	{
	case 1:	// bin in t
		hist_axis = 0;
		bin_axis_1 = 1;
		bin_axis_2 = 2;
		hist_binsize = (p_config_vals->t_bin);
		bin_1_binsize = (p_config_vals->x_bin);
		bin_2_binsize = (p_config_vals->y_bin);
		break;
	case 2:	// bin in x
		hist_axis = 1;
		bin_axis_1 = 0;
		bin_axis_2 = 2;
		hist_binsize = (p_config_vals->x_bin);
		bin_1_binsize = (p_config_vals->t_bin);
		bin_2_binsize = (p_config_vals->y_bin);
		break;
	case 3:	// bin in y
		hist_axis = 2;
		bin_axis_1 = 0;
		bin_axis_2 = 1;
		hist_binsize = (p_config_vals->y_bin);
		bin_1_binsize = (p_config_vals->t_bin);
		bin_2_binsize = (p_config_vals->x_bin);
		break;
	}

	double dh_12, d1_12, d2_12, dh_13, d1_13, d2_13, dh_23, d1_23, d2_23;
	int hist_bin_12, hist_bin_13, hist_bin_23;
	long number_pairs = 0, number_triples = 0;	// count these at the original bin size

	// the underscore is how big the transverse bins are
	//     _0 is the specified bin size
	//     _1 is 20% larger bin
	//     _m2 is 40% smaller bin
	// and so on

	// numerators
	umat g2_0(p_config_vals->num_bins,1);
	g2_0.fill(0);	
	uword* p_g2_0 = g2_0.memptr();
	umat g2_1(p_config_vals->num_bins,1);
	g2_1.fill(0);
	uword* p_g2_1 = g2_1.memptr();
	umat g2_2(p_config_vals->num_bins,1);
	g2_2.fill(0);
	uword* p_g2_2 = g2_2.memptr();
	umat g2_3(p_config_vals->num_bins,1);
	g2_3.fill(0);
	uword* p_g2_3 = g2_3.memptr();
	umat g2_m1(p_config_vals->num_bins,1);
	g2_m1.fill(0);
	uword* p_g2_m1 = g2_m1.memptr();
	umat g2_m2(p_config_vals->num_bins,1);
	g2_m2.fill(0);
	uword* p_g2_m2 = g2_m2.memptr();
	umat g2_m3(p_config_vals->num_bins,1);
	g2_m3.fill(0);
	uword* p_g2_m3 = g2_m3.memptr();

	umat g3_0(p_config_vals->num_bins,p_config_vals->num_bins);
	g3_0.fill(0);	
	uword* p_g3_0 = g3_0.memptr();
	umat g3_1(p_config_vals->num_bins,p_config_vals->num_bins);
	g3_1.fill(0);
	uword* p_g3_1 = g3_1.memptr();
	umat g3_2(p_config_vals->num_bins,p_config_vals->num_bins);
	g3_2.fill(0);
	uword* p_g3_2 = g3_2.memptr();
	umat g3_3(p_config_vals->num_bins,p_config_vals->num_bins);
	g3_3.fill(0);
	uword* p_g3_3 = g3_3.memptr();
	umat g3_m1(p_config_vals->num_bins,p_config_vals->num_bins);
	g3_m1.fill(0);
	uword* p_g3_m1 = g3_m1.memptr();
	umat g3_m2(p_config_vals->num_bins,p_config_vals->num_bins);
	g3_m2.fill(0);
	uword* p_g3_m2 = g3_m2.memptr();
	umat g3_m3(p_config_vals->num_bins,p_config_vals->num_bins);
	g3_m3.fill(0);
	uword* p_g3_m3 = g3_m3.memptr();

	//denominators
	umat g2_den_0(p_config_vals->num_bins,1);
	g2_den_0.fill(0);	
	uword* p_g2_den_0 = g2_den_0.memptr();
	umat g2_den_1(p_config_vals->num_bins,1);
	g2_den_1.fill(0);
	uword* p_g2_den_1 = g2_den_1.memptr();
	umat g2_den_2(p_config_vals->num_bins,1);
	g2_den_2.fill(0);
	uword* p_g2_den_2 = g2_den_2.memptr();
	umat g2_den_3(p_config_vals->num_bins,1);
	g2_den_3.fill(0);
	uword* p_g2_den_3 = g2_den_3.memptr();
	umat g2_den_m1(p_config_vals->num_bins,1);
	g2_den_m1.fill(0);
	uword* p_g2_den_m1 = g2_den_m1.memptr();
	umat g2_den_m2(p_config_vals->num_bins,1);
	g2_den_m2.fill(0);
	uword* p_g2_den_m2 = g2_den_m2.memptr();
	umat g2_den_m3(p_config_vals->num_bins,1);
	g2_den_m3.fill(0);
	uword* p_g2_den_m3 = g2_den_m3.memptr();

	umat g3_den_0(p_config_vals->num_bins,p_config_vals->num_bins);
	g3_den_0.fill(0);	
	uword* p_g3_den_0 = g3_den_0.memptr();
	umat g3_den_1(p_config_vals->num_bins,p_config_vals->num_bins);
	g3_den_1.fill(0);	
	uword* p_g3_den_1 = g3_den_1.memptr();
	umat g3_den_2(p_config_vals->num_bins,p_config_vals->num_bins);
	g3_den_2.fill(0);	
	uword* p_g3_den_2 = g3_den_2.memptr();
	umat g3_den_3(p_config_vals->num_bins,p_config_vals->num_bins);
	g3_den_3.fill(0);	
	uword* p_g3_den_3 = g3_den_3.memptr();
	umat g3_den_m1(p_config_vals->num_bins,p_config_vals->num_bins);
	g3_den_m1.fill(0);	
	uword* p_g3_den_m1 = g3_den_m1.memptr();
	umat g3_den_m2(p_config_vals->num_bins,p_config_vals->num_bins);
	g3_den_m2.fill(0);	
	uword* p_g3_den_m2 = g3_den_m2.memptr();
	umat g3_den_m3(p_config_vals->num_bins,p_config_vals->num_bins);
	g3_den_m3.fill(0);	
	uword* p_g3_den_m3 = g3_den_m3.memptr();

	int size_file_list = p_config_vals->lastfile - p_config_vals->firstfile + 1;
	long* file_list = new long[size_file_list];
	long* nvf_windowed = new long[size_file_list];
	double* t0 = new double[size_file_list];

	umat nvf_windowed_arma(size_file_list,1); 
	nvf_windowed_arma.fill(0);	

	char buffer[100];

	long cumulative_nvf = 0;
	int number_files_on_list = 0;

	char filepath_nvf [200];
	strncpy(filepath_nvf,p_config_vals->filepath, 200);
	strcat(filepath_nvf,"_nvf_windowed.txt");

	ifstream file_nvf(filepath_nvf);
	if (!file_nvf)
	{
		printf("ERROR: cannot find file nvf_windowed in current folder \n");
		cout << filepath_nvf << '\n';
		file_nvf.close();
		return 1;
	}



	for (int line = 0; line < (p_config_vals->lastfile - p_config_vals->firstfile + 1); line++)
	{
		file_nvf.getline(buffer,100,'\n');
		nvf_windowed[line] = atoi(buffer);
		nvf_windowed_arma(line) = atoi(buffer);
		t0[line]=0;
	}	
	file_nvf.close();

	long size_data_arrays = 10000000;	// 10 million counts in denominator txy array

	if (p_config_vals->num_files_in_denom > 0)
	{
		size_data_arrays = p_config_vals->num_files_in_denom * nvf_windowed_arma.max();	
	}

	mat txy_all(size_data_arrays,3);
	txy_all.fill(999.999);					// fill with rubbish values

	if (p_config_vals->temp_window_pc != 0 || p_config_vals->num_window_pc != 0)
	{
		char filepath_temp [200];
		strncpy(filepath_temp,p_config_vals->filepath, 200);
		strcat(filepath_temp,"_file_list.txt");

		ifstream file_to_read(filepath_temp);
		if (!file_to_read)
		{
			printf("ERROR: cannot find file list in current folder \n");
			cout << filepath_temp << '\n';
			file_to_read.close();
			return 1;
		}

		for (int line = 0; line < (p_config_vals->lastfile - p_config_vals->firstfile + 1); line++)
		{
			file_to_read.getline(buffer,100,'\n');
			file_list[line] = atoi(buffer);
			number_files_on_list++;
		}	
		file_to_read.close();
	}
	else
	{
		for (int i = p_config_vals->firstfile; i <= p_config_vals->lastfile; i++)
		{
			file_list[i-1] = i;
		}
		number_files_on_list = p_config_vals->lastfile - p_config_vals->firstfile + 1;
	}

	if (p_config_vals->fit_remove_t0 == 1)
	{
		cout << "feature currently being added..." << endl;

		char filepath_temp [200];
		strncpy(filepath_temp,p_config_vals->filepath, 200);
		strcat(filepath_temp,"_t0.txt");

		ifstream file_to_read(filepath_temp);
		if (!file_to_read)
		{
			printf("ERROR: cannot find file list in current folder \n");
			cout << filepath_temp << '\n';
			file_to_read.close();
			return 1;
		}

		for (int line = 0; line < (p_config_vals->lastfile - p_config_vals->firstfile + 1); line++)
		{
			file_to_read.getline(buffer,100,';');
			t0[line] = atof(buffer);
			//cout << t0[line] << endl;
		}	
		file_to_read.close();
	}

	int line = 0;

	for (int line_dummy = 0; line_dummy < number_files_on_list * 2; line_dummy++)
	{
		if (line >= number_files_on_list )
		{
			break;
		}

		char filepath_temp_data_t [200];
		strncpy(filepath_temp_data_t,p_config_vals->filepath, 200);
		strcat(filepath_temp_data_t,"_txy_forc_AGM_correlations_fwrite_t");
		strcat(filepath_temp_data_t,itoa(file_list[line],buffer,10));
		strcat(filepath_temp_data_t,".bin");

		char filepath_temp_data_x [200];
		strncpy(filepath_temp_data_x,p_config_vals->filepath, 200);
		strcat(filepath_temp_data_x,"_txy_forc_AGM_correlations_fwrite_x");
		strcat(filepath_temp_data_x,itoa(file_list[line],buffer,10));
		strcat(filepath_temp_data_x,".bin");

		char filepath_temp_data_y [200];
		strncpy(filepath_temp_data_y,p_config_vals->filepath, 200);
		strcat(filepath_temp_data_y,"_txy_forc_AGM_correlations_fwrite_y");
		strcat(filepath_temp_data_y,itoa(file_list[line],buffer,10));
		strcat(filepath_temp_data_y,".bin");

		ifstream file_to_read_data(filepath_temp_data_t);
		if (!file_to_read_data)
		{
			printf("ERROR: cannot find data file for correlations in current folder \n");
			cout << filepath_temp_data_t << '\n';
			file_to_read_data.close();
			line++;
			//return 1;
			continue;
		}

		cout << filepath_temp_data_t << " with " << nvf_windowed[file_list[line] - 1] << " hits"<< '\n';

		if (cumulative_nvf + nvf_windowed[file_list[line] - 1] > size_data_arrays)
		{
			cout << "txy_all array for denominator exceeding " << itoa(size_data_arrays,buffer,10) << " elements" << endl;
			cout << "will calculate demoninator and reset txy_all array" << endl;

			// calculate denominators, add, and reset txy_all

			if (p_config_vals->num_CPU < 2 || p_config_vals->num_CPU > 8)	// single core
			{
				uvec indices = sort_index(txy_all.col(0));
				mat txy_sorted(cumulative_nvf,3);
				txy_sorted.fill(999.999);

				for (int m = 0; m < cumulative_nvf; m++)
				{
					txy_sorted.row(m) = txy_all.row(indices(m));
				}

				double* p_txy_sorted = txy_sorted.memptr();

				int denominator_result = crunch_gn_single_core(p_config_vals, cumulative_nvf, 
					hist_axis, bin_axis_1, bin_axis_2,
					hist_binsize, bin_1_binsize, bin_2_binsize,
					p_txy_sorted,
					p_g2_den_0, p_g2_den_1, p_g2_den_2, p_g2_den_3, p_g2_den_m1, p_g2_den_m2, p_g2_den_m3,
					p_g3_den_0, p_g3_den_1, p_g3_den_2, p_g3_den_3, p_g3_den_m1, p_g3_den_m2, p_g3_den_m3);

				cumulative_nvf = 0;
				continue;
			}
			else if (p_config_vals->num_CPU > 1)	// multi core
			{
				uvec indices = sort_index(txy_all.col(0));
				mat txy_sorted(cumulative_nvf,3);
				txy_sorted.fill(999.999);

				for (int m = 0; m < cumulative_nvf; m++)
				{
					txy_sorted.row(m) = txy_all.row(indices(m));
				}

				double* p_txy_sorted = txy_sorted.memptr();

				int denominator_result = crunch_gn_multi_core(p_config_vals, cumulative_nvf, 
					hist_axis, bin_axis_1, bin_axis_2,
					hist_binsize, bin_1_binsize, bin_2_binsize,
					p_txy_sorted,
					p_g2_den_0, p_g2_den_1, p_g2_den_2, p_g2_den_3, p_g2_den_m1, p_g2_den_m2, p_g2_den_m3,
					p_g3_den_0, p_g3_den_1, p_g3_den_2, p_g3_den_3, p_g3_den_m1, p_g3_den_m2, p_g3_den_m3);

				cumulative_nvf = 0;
				continue;
			}
			continue;
		}

		time_t s1, s2, s3, s4;
		s1 = time(NULL);

		long t_index = 0;
		long x_index = 0;
		long y_index = 0;

		mat txy_this_file(nvf_windowed[file_list[line] - 1],3);

		// read t vector
		FILE *ptr_myfile_t;
		double data_to_read_t;
		errno_t err_t;
		err_t = fopen_s(&ptr_myfile_t,filepath_temp_data_t,"rb");
		if (!ptr_myfile_t)
		{
			std::cout << ("Unable to open t file!");
			line++;
			continue;
			//return 1;
		}

			//if (p_config_vals->fit_remove_t0 == 1)
			//{
			//	cout << t0[line] << endl;
			//	system("pause");
			//}

		while (nvf_windowed[file_list[line] - 1] - t_index > 0) 
		{		
			fread(&data_to_read_t,sizeof(data_to_read_t),1,ptr_myfile_t);
			if (p_config_vals->fit_remove_t0 == 1)
			{
				txy_this_file(t_index,0) = data_to_read_t - t0[line];
				txy_all(t_index + cumulative_nvf,0) = data_to_read_t - t0[line];
			}
			if (p_config_vals->fit_remove_t0 == 0)
			{
				txy_this_file(t_index,0) = data_to_read_t;
				txy_all(t_index + cumulative_nvf,0) = data_to_read_t;
			}
			t_index++;
		}
		fclose(ptr_myfile_t);

		// read x vector
		FILE *ptr_myfile_x;
		double data_to_read_x;
		errno_t err_x;
		err_x = fopen_s(&ptr_myfile_x,filepath_temp_data_x,"rb");
		if (!ptr_myfile_x)
		{
			std::cout << ("Unable to open x file!");
			line++;
			continue;
			//return 1;
		}

		while (nvf_windowed[file_list[line] - 1] - x_index > 0) 
		{		
			fread(&data_to_read_x,sizeof(data_to_read_x),1,ptr_myfile_x);
			txy_this_file(x_index,1) = data_to_read_x;
			txy_all(x_index + cumulative_nvf,1) = data_to_read_x;
			x_index++;
		}
		fclose(ptr_myfile_x);

		// read y vector
		FILE *ptr_myfile_y;
		double data_to_read_y;
		errno_t err_y;
		err_y = fopen_s(&ptr_myfile_y,filepath_temp_data_y,"rb");
		if (!ptr_myfile_y)
		{
			std::cout << ("Unable to open y file!");
			line++;
			continue;
			//return 1;
		}

		while (nvf_windowed[file_list[line] - 1] - y_index > 0) 
		{		
			fread(&data_to_read_y,sizeof(data_to_read_y),1,ptr_myfile_y);
			txy_this_file(y_index,2) = data_to_read_y;
			txy_all(y_index + cumulative_nvf,2) = data_to_read_y;
			y_index++;
		}
		fclose(ptr_myfile_y);

		// calculate numerators

		double* p_txy_this_file = txy_this_file.memptr();

		if (p_config_vals->num_CPU < 2 || p_config_vals->num_CPU > 8)	// single core
		{
			int numerator_result = crunch_gn_single_core(p_config_vals, nvf_windowed[file_list[line] - 1], 
				hist_axis, bin_axis_1, bin_axis_2,
				hist_binsize, bin_1_binsize, bin_2_binsize,
				p_txy_this_file,
				p_g2_0, p_g2_1, p_g2_2, p_g2_3, p_g2_m1, p_g2_m2, p_g2_m3,
				p_g3_0, p_g3_1, p_g3_2, p_g3_3, p_g3_m1, p_g3_m2, p_g3_m3);
		}
		else if (p_config_vals->num_CPU > 1)	// multi core
		{
			int numerator_result = crunch_gn_multi_core(p_config_vals, nvf_windowed[file_list[line] - 1], 
				hist_axis, bin_axis_1, bin_axis_2,
				hist_binsize, bin_1_binsize, bin_2_binsize,
				p_txy_this_file,
				p_g2_0, p_g2_1, p_g2_2, p_g2_3, p_g2_m1, p_g2_m2, p_g2_m3,
				p_g3_0, p_g3_1, p_g3_2, p_g3_3, p_g3_m1, p_g3_m2, p_g3_m3);
		}

		cumulative_nvf = cumulative_nvf + nvf_windowed[file_list[line]-1];
		file_to_read_data.close();
		line++;
	}

	//calculate denominators

	uvec indices = sort_index(txy_all.col(0));
	mat txy_sorted(cumulative_nvf,3);
	txy_sorted.fill(999.999);

	for (int m = 0; m < cumulative_nvf; m++)
	{
		txy_sorted.row(m) = txy_all.row(indices(m));
	}

	double* p_txy_sorted = txy_sorted.memptr();

	if (p_config_vals->num_CPU < 2 || p_config_vals->num_CPU > 8)	// single core
	{
		int denominator_result = crunch_gn_single_core(p_config_vals, cumulative_nvf, 
			hist_axis, bin_axis_1, bin_axis_2,
			hist_binsize, bin_1_binsize, bin_2_binsize,
			p_txy_sorted,
			p_g2_den_0, p_g2_den_1, p_g2_den_2, p_g2_den_3, p_g2_den_m1, p_g2_den_m2, p_g2_den_m3,
			p_g3_den_0, p_g3_den_1, p_g3_den_2, p_g3_den_3, p_g3_den_m1, p_g3_den_m2, p_g3_den_m3);
	}
	else if (p_config_vals->num_CPU > 1)	// multi core
	{
		int denominator_result = crunch_gn_multi_core(p_config_vals, cumulative_nvf, 
			hist_axis, bin_axis_1, bin_axis_2,
			hist_binsize, bin_1_binsize, bin_2_binsize,
			p_txy_sorted,
			p_g2_den_0, p_g2_den_1, p_g2_den_2, p_g2_den_3, p_g2_den_m1, p_g2_den_m2, p_g2_den_m3,
			p_g3_den_0, p_g3_den_1, p_g3_den_2, p_g3_den_3, p_g3_den_m1, p_g3_den_m2, p_g3_den_m3);
	}

	// save outputs
	char filepath_save [200];
	strncpy(filepath_save,p_config_vals->filepath, 200);

	g2_0.save(std::string(filepath_save) + "_g2_0.txt", raw_ascii);
	g2_1.save(std::string(filepath_save) + "_g2_1.txt", raw_ascii);
	g2_2.save(std::string(filepath_save) + "_g2_2.txt", raw_ascii);
	g2_3.save(std::string(filepath_save) + "_g2_3.txt", raw_ascii);
	g2_m1.save(std::string(filepath_save) + "_g2_m1.txt", raw_ascii);
	g2_m2.save(std::string(filepath_save) + "_g2_m2.txt", raw_ascii);
	g2_m3.save(std::string(filepath_save) + "_g2_m3.txt", raw_ascii);
	g2_den_0.save(std::string(filepath_save) + "_g2_den_0.txt", raw_ascii);
	g2_den_1.save(std::string(filepath_save) + "_g2_den_1.txt", raw_ascii);
	g2_den_2.save(std::string(filepath_save) + "_g2_den_2.txt", raw_ascii);
	g2_den_3.save(std::string(filepath_save) + "_g2_den_3.txt", raw_ascii);
	g2_den_m1.save(std::string(filepath_save) + "_g2_den_m1.txt", raw_ascii);
	g2_den_m2.save(std::string(filepath_save) + "_g2_den_m2.txt", raw_ascii);
	g2_den_m3.save(std::string(filepath_save) + "_g2_den_m3.txt", raw_ascii);

	if (p_config_vals->order > 2)
	{
		g3_0.save(std::string(filepath_save) + "_g3_0.txt", raw_ascii);
		g3_1.save(std::string(filepath_save) + "_g3_1.txt", raw_ascii);
		g3_2.save(std::string(filepath_save) + "_g3_2.txt", raw_ascii);
		g3_3.save(std::string(filepath_save) + "_g3_3.txt", raw_ascii);
		g3_m1.save(std::string(filepath_save) + "_g3_m1.txt", raw_ascii);
		g3_m2.save(std::string(filepath_save) + "_g3_m2.txt", raw_ascii);
		g3_m3.save(std::string(filepath_save) + "_g3_m3.txt", raw_ascii);
		g3_den_0.save(std::string(filepath_save) + "_g3_den_0.txt", raw_ascii);
		g3_den_1.save(std::string(filepath_save) + "_g3_den_1.txt", raw_ascii);
		g3_den_2.save(std::string(filepath_save) + "_g3_den_2.txt", raw_ascii);
		g3_den_3.save(std::string(filepath_save) + "_g3_den_3.txt", raw_ascii);
		g3_den_m1.save(std::string(filepath_save) + "_g3_den_m1.txt", raw_ascii);
		g3_den_m2.save(std::string(filepath_save) + "_g3_den_m2.txt", raw_ascii);
		g3_den_m3.save(std::string(filepath_save) + "_g3_den_m3.txt", raw_ascii);
	}

	delete [] file_list;
	delete [] nvf_windowed;
	return 0;
}


int _tmain(int argc, _TCHAR* argv[])
{
	seconds_start = time(NULL);

	// read config file written by MATLAB front panel into a struct
	struct config_vals_struct config_vals;
	int read_config_result = read_config(&config_vals);

	if ( config_vals.convert_DLD_txy == 1 && config_vals.ready_for_gn == 0)
	{
		struct channel_triggers_struct channel_triggers;
		channel_triggers.ch0 = 0;
		channel_triggers.ch1 = 0;
		channel_triggers.ch2 = 0;
		channel_triggers.ch3 = 0;
		channel_triggers.rubbish = 0;

		int convert_DLD_txy_result = DLD_txy(&config_vals,&channel_triggers);

		FILE *out_stats;
		char outfile_stats_handle[200] = "";
		strncpy(outfile_stats_handle,config_vals.filepath, 150);
		strcat(outfile_stats_handle,"_channels.txt");
		out_stats = fopen(outfile_stats_handle,"wt");
		fprintf(out_stats,"%.9f; \n",channel_triggers.ch0);
		fprintf(out_stats,"%.9f; \n",channel_triggers.ch1);
		fprintf(out_stats,"%.9f; \n",channel_triggers.ch2);
		fprintf(out_stats,"%.9f; \n",channel_triggers.ch3);
		fprintf(out_stats,"%.9f; \n",channel_triggers.rubbish);
		fclose(out_stats);
	}

	if (config_vals.ready_for_gn == 1)
	{
		int calculate_gn_result = calculate_gn(&config_vals);
	}


	seconds_end = time(NULL);

	cout << "\nTotal calculation took " << seconds_end - seconds_start << " seconds.\n" << endl;

	//system("pause");
	return 0;
}
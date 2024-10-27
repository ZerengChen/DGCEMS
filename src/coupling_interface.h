#pragma once
#include<iostream>
#include <mpi.h>
#include "MeshUnion.h"

extern int* dg_demo_comp_id;
extern int* decomp_id, *grid_h2d_id;
extern int* time_step;
extern int* coupling_freq;
extern MeshUnion *meshunion;
// extern float *HS;
// extern float *DIR;
// extern float *RTP;
// extern float *QB;
// extern float *Hm;
// extern float *Um;
// extern float *Vm;

// #define LINK_WITHOUT_UNDERLINE
/*1*/
#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void initialize_ccpl_mgrs
#else
extern "C" void initialize_ccpl_mgrs_
#endif
();

/*2*/
#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_root_component(MPI_Fint *f_comm, const char *comp_name, 
 const char *local_comp_type, const char *annotation, int *comp_id,
 int *enabled_in_parent_coupling_gen, int *change_dir, const char *executable_name);
#else
extern "C" void register_root_component_(MPI_Fint *f_comm, const char *comp_name, const char *local_comp_type, const char *annotation, int *comp_id,
 int *enabled_in_parent_coupling_gen, int *change_dir, const char *executable_name);
#endif

/*3*/
#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void set_component_time_step
#else
extern "C" void set_component_time_step_
#endif
(int *comp_id, int *time_step_in_second, const char *annotation);

/*4*/
#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_h2d_grid_with_global_data
#else
extern "C" void register_h2d_grid_with_global_data_
#endif
(int *comp_id, int *grid_id, const char *grid_name, const char *edge_type, const char *coord_unit, const char *cyclic_or_acyclic, const char *data_type, int *dim_size1, int *dim_size2, int *size_center_lon, int *size_center_lat, int *size_mask, int *size_area, int *size_vertex_lon, int *size_vertex_lat, char *min_lon, char *max_lon, char *min_lat, char *max_lat, char *center_lon, char *center_lat, int *mask, char *area, char *vertex_lon, char *vertex_lat, const char *annotation);


/*5*/
#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_parallel_decomposition
#else
extern "C" void register_parallel_decomposition_
#endif
(int *decomp_id, int *grid_id, int *num_local_cells, int *array_size, const int *local_cells_global_indx, const char *decomp_name, const char *annotation);


/*6*/
#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_external_field_instance
#else
extern "C" void register_external_field_instance_
#endif
(int *field_instance_id, const char *field_name, long *data_buffer_ptr, int *field_size, int *decomp_id, int *comp_or_grid_id, int *buf_mark, int *usage_tag, const char *unit, const char *data_type, const char *annotation);


/*7*/
#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void define_single_timer
#else
extern "C" void define_single_timer_
#endif
(int *comp_id, int *timer_id, const char *freq_unit, int *freq_count, int *local_lag_count, int *remote_lag_count, const char *annotation);

/*8*/
#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void register_inout_interface
#else
extern "C" void register_inout_interface_
#endif
(const char *interface_name, int *interface_id, int *import_or_export, int *num_fields, int *field_ids, int *timer_id, int *inst_or_aver, const char *annotation, int *array_size1);

/*9*/
#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void execute_inout_interface_with_name
#else
extern "C" void execute_inout_interface_with_name_
#endif
(int *comp_id, const char *interface_name, int *bypass_timer, int *field_update_status, int *size_field_update_status, int *num_dst_fields, const char *annotation);

/*10*/
#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void finalize_ccpl
#else
extern "C" void finalize_ccpl_
#endif
(int *to_finalize_MPI, const char *annotation);

/*11*/
#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void advance_component_time
#else
extern "C" void advance_component_time_
#endif
(int *comp_id, const char *annotation);

/*12*/
#ifdef LINK_WITHOUT_UNDERLINE
extern "C" void ccpl_end_registration
#else
extern "C" void ccpl_end_registration_
#endif
(int *comp_id, const char * annotation);

// extern "C"
// {
//     int ccpl_interface_mod_mp_ccpl_register_component_(int *parent_id,const char* comp_name, const char* comp_type, int* comp_comm, bool* considered_in_parent_coupling_gen, bool* change_dir,const char* annotation);
//     void ccpl_interface_mod_mp_ccpl_set_normal_time_step_(int *comp_id,int *time_step_in_second,const char* annotation);
//     int ccpl_interface_mod_mp_ccpl_register_h2d_grid_global_online_c1d_m1d_float_(int* comp_id, const char *grid_name,const char*edge_type,const char* coord_unit, const char* cyclic_or_acyclic, int *dim_size1,int *dim_size2,float * min_lon);
// }
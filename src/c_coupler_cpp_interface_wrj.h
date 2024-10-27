#pragma once
/*1*/
void register_dg_demo_component_wrj(int *comp_comm);
/*2*/
void register_component_coupling_configuration_wrj(float* HS_from_swan, float* T_from_swan, float* DIR_from_swan, float *QB_from_swan, float* WLEN_from_swan, float* UBOT_from_swan, float* TMBOT_from_swan, float *H_to_swan, float* U_to_swan,float* V_to_swan,float* test_to_swan);
/*3*/
int register_H2D_grid_global_online_C1D_M1D_float_wrj(int *comp_id, const char* grid_name, const char* edge_type, const char* coord_unit, const char* cyclic_or_acyclic, int* dim_size1, int* dim_size2, float *min_lon, float * max_lon, float * min_lat, float * max_lat, float* center_lon, float* center_lat, const char* annotation);
/*4*/
int register_normal_parallel_decomp_wrj(const char* decomp_name,int * grid_id,int*  num_local_cells,int* local_cells_global_index,int *size_local_cells_global_index,const char* annotation);
/*5*/
int define_single_timer_wrj(int* comp_id, const char* period_unit, int* period_count, int* local_lag_count, int* remote_lag_count,const char* annotation);
/*6*/
int register_export_interface_wrj(const char* interface_name, int* num_field_instances, int* field_instance_IDs, int* timer_ID, const char* annotation);
/*7*/
int register_import_interface_wrj(const char* interface_name,int* num_field_instances,int* field_instance_IDs, int* timer_ID, int* inst_or_aver,const char* annotation);
/*8*/
bool execute_interface_using_name_wrj(int* component_id, const char* interface_name, bool bypass_timer, const char* annotation);
/*9*/
void advance_time_wrj(int *comp_id,const char* annotation);
/*10*/
void  end_coupling_configuration_wrj(int* comp_id, const char* annotation);
/*11*/
void finalize_wrj(bool to_finalize_MPI, const char* annotation);

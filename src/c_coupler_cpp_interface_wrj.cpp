#include "coupling_interface.h"


/*1注册分量模式*/
void register_dg_demo_component_wrj(int *comp_comm)
{   

    int comp_id; 
    int considered_in_parent_coupling_gen = 1;
    int change_dir= 1;

    initialize_ccpl_mgrs_();
    register_root_component_(comp_comm,"dg_demo","ocn","register ocn model dg_demo",&comp_id,&considered_in_parent_coupling_gen,&change_dir,"./dg_demo");
    *dg_demo_comp_id = comp_id;
}


void set_normal_time_step_wrj(int* comp_id, int* time_step_in_second, const char* annotation)
{
    set_component_time_step_(comp_id,time_step_in_second,annotation);
}

/*3*/
int register_H2D_grid_global_online_C1D_M1D_float_wrj(int *comp_id, const char* grid_name, const char* edge_type, const char* coord_unit, const char* cyclic_or_acyclic, int* dim_size1, int* dim_size2,
                                                           float *min_lon, float * max_lon, float * min_lat, float * max_lat, float* center_lon, float* center_lat, const char* annotation)
{
    int grid_id;

    int size_center_lon = *dim_size1;
    int size_center_lat = *dim_size1;
    int size_mask =-1;
    int size_area =-1;
    int size_vertex_lon =-1;
    int size_vertex_lat = -1;

    int *temp_int = new int;
    float *temp_float_1d = new float;
    float *temp_float_2d = new float;

    int *temp_mask = temp_int;
    float *temp_area = temp_float_1d;
    float *temp_vertex_lon = temp_float_2d;
    float *temp_vertex_lat = temp_float_2d;
    
    const char* real4= "real4";
    register_h2d_grid_with_global_data_(comp_id, &grid_id, grid_name, edge_type, coord_unit, cyclic_or_acyclic, 
                                        real4, dim_size1, dim_size2, &size_center_lon, &size_center_lat, &size_mask, &size_area, &size_vertex_lon, &size_vertex_lat,
                                        (char*)min_lon, (char*)max_lon, (char*)min_lat, (char*)max_lat, (char*)center_lon, (char*)center_lat,
                                        temp_mask, (char*)temp_area, (char*)temp_vertex_lon, (char*)temp_vertex_lat, annotation);

    delete  temp_int;
    delete  temp_float_1d;
    delete  temp_float_2d;

    return grid_id;

}

/*4*/
int register_normal_parallel_decomp_wrj(const char* decomp_name,int * grid_id,int*  num_local_cells,int* local_cells_global_index,int *size_local_cells_global_index,const char* annotation) 
{
    int decomp_id;

    register_parallel_decomposition_(&decomp_id, grid_id, num_local_cells, size_local_cells_global_index, local_cells_global_index, decomp_name, annotation);
    
    return decomp_id;
}



/*4*/
int register_model_float_1D_data_wrj(float *data_buf,int *size_data_buf, const char* field_name, int* decomp_id, int* comp_or_grid_id, int* buf_mark,
                                         int* usage_tag, const char* field_unit, const char* annotation)
{                                         
    int field_instance_id;
    long *data_buffer_ptr=(long*)&data_buf;

    register_external_field_instance_(&field_instance_id ,field_name, data_buffer_ptr, size_data_buf, decomp_id, comp_or_grid_id, buf_mark, usage_tag, field_unit,"real4", annotation);
    
    int CCPL_register_model_float_1D_data = field_instance_id;
    return CCPL_register_model_float_1D_data;

}

/*5*/
int define_single_timer_wrj(int* comp_id, const char* period_unit, int* period_count, int* local_lag_count, int* remote_lag_count,const char* annotation) 
{
    int timer_id;

    define_single_timer_(comp_id, &timer_id, period_unit, period_count, local_lag_count, remote_lag_count, annotation);

    int CCPL_define_single_timer = timer_id; 

    return CCPL_define_single_timer;
}

/*6*/
int register_export_interface_wrj(const char* interface_name, int* num_field_instances, int* field_instance_IDs, int* timer_ID, const char* annotation)
{   
    int interface_id;
    int import_or_export=1;
    int size_fields=5;
    

    register_inout_interface_(interface_name, &interface_id, &import_or_export, num_field_instances, field_instance_IDs, timer_ID, field_instance_IDs, annotation,&size_fields);

    int CCPL_register_export_interface = interface_id;
    return CCPL_register_export_interface;
}


/*7*/
int register_import_interface_wrj(const char* interface_name,int* num_field_instances,int* field_instance_IDs, int* timer_ID, int* inst_or_aver, const char* annotation)
{
   
    int interface_id;
    int import_or_export = 0;
    int size_fields=8;
   
    register_inout_interface_(interface_name, &interface_id, &import_or_export, num_field_instances, field_instance_IDs, timer_ID, inst_or_aver, annotation, &size_fields);

   int CCPL_register_import_interface = interface_id;
   return CCPL_register_import_interface;
}



/*8*/
bool execute_interface_using_name_wrj(int* component_id, const char* interface_name, bool bypass_timer, const char* annotation)
{
    int temp_field_update_status[4096];
    int local_bypass_timer;
    int num_dst_fields;
    int size_field_update_status= 4000;
    bool CCPL_execute_interface_using_name;

    if (bypass_timer)
    {
        local_bypass_timer = 1;
    }
    else
    {
        local_bypass_timer = 0;
    }
    
    execute_inout_interface_with_name_(component_id, interface_name, &local_bypass_timer, temp_field_update_status, &size_field_update_status, &num_dst_fields, annotation);
    
    CCPL_execute_interface_using_name = true;
    return CCPL_execute_interface_using_name;
}

/*9*/
void advance_time_wrj(int *comp_id,const char* annotation)
{

    advance_component_time_(comp_id, annotation);
}

/*10*/
void  end_coupling_configuration_wrj(int* comp_id, const char* annotation)
{
    
    ccpl_end_registration_(comp_id, annotation);
}

/*11*/
void finalize_wrj(bool to_finalize_MPI, const char* annotation)
{
    int local_to_finalize_MPI;
    if (to_finalize_MPI)
    {
        local_to_finalize_MPI = 1;
    }else
    {
        local_to_finalize_MPI = 0;
    }
    
    finalize_ccpl_(&local_to_finalize_MPI,annotation);   
}

/*2注册耦合变量实例*/
void register_component_coupling_configuration_wrj(float* HS_from_swan, float* T_from_swan, float* DIR_from_swan, float *QB_from_swan, float *WLEN_from_swan, float *UBOT_from_swan, float *TMBOT_from_swan, float *H_to_swan, float* U_to_swan,float* V_to_swan,float* test_to_swan)//需要注意指针地址
{

    int field_id_HS, field_id_DIR, field_id_RTP, field_id_QB, field_id_WLEN, field_id_UBOT, field_id_TMBOT;
    int field_id_H, field_id_U, field_id_V,field_id_test;
    int fields_id_s[8],fields_id_r[8];
    int timer_id;
    int dim_size2=0;
    float min_lon = 0.0;
    float max_lon = 360.0;
    float min_lat = 0.0;
    float max_lat = 90;
    
    float *lon = new float[*meshunion->mesh2d_p->Nv2d];

    float *lat = new float[*meshunion->mesh2d_p->Nv2d];
    int *local_grid_cell_index=new int[*meshunion->mesh2d_p->Nv2d];

    int buf_mark = 0;
    int usage_tag = 3;
    int local_lag_count=0;
    int remote_lag_count=0;


    for (size_t i = 0; i < *meshunion->mesh2d_p->Nv2d; i++)
    {
        lon[i] = (float)((float)360.0/(float)(i+1));
        lat[i] = (float)((float)90.0/(float)(i+1));
        local_grid_cell_index[i] = i + 1;
    }
    
	//注册时间步长，水平网格与并行剖分

    set_normal_time_step_wrj(dg_demo_comp_id, time_step, "setting the time step for dg_demo");
    *grid_h2d_id = register_H2D_grid_global_online_C1D_M1D_float_wrj(dg_demo_comp_id, "dg_demo_H2D_grid", "LON_LAT", "degrees", "acyclic", meshunion->mesh2d_p->Nv2d, &dim_size2, &min_lon, &max_lon, &min_lat, &max_lat, lon, lat, "register dg_demo H2D grid");
    *decomp_id = register_normal_parallel_decomp_wrj("decomp_dg_demo_grid", grid_h2d_id, meshunion->mesh2d_p->Nv2d, local_grid_cell_index, meshunion->mesh2d_p->Nv2d, "allocate decomp for dg_demo grid");


// ------------register field instances to C-Coupler2--------------
    field_id_HS = register_model_float_1D_data_wrj(HS_from_swan, meshunion->mesh2d_p->Nv2d, "HS", decomp_id, grid_h2d_id, &buf_mark, &usage_tag, "m", "register field instance of Significant wave height");
    field_id_DIR = register_model_float_1D_data_wrj(DIR_from_swan, meshunion->mesh2d_p->Nv2d, "DIR", decomp_id, grid_h2d_id, &buf_mark, &usage_tag, "d", "register field instance of Mean wave direction");
    field_id_RTP = register_model_float_1D_data_wrj(T_from_swan, meshunion->mesh2d_p->Nv2d, "RTP", decomp_id, grid_h2d_id, &buf_mark, &usage_tag, "s", "register field instance of Relative peak wave period");  
    field_id_QB = register_model_float_1D_data_wrj(QB_from_swan, meshunion->mesh2d_p->Nv2d, "QB", decomp_id, grid_h2d_id, &buf_mark, &usage_tag, "1", "register field instance of Fraction of breakers in expression of Battjes and Janssen");
	field_id_WLEN = register_model_float_1D_data_wrj(WLEN_from_swan, meshunion->mesh2d_p->Nv2d, "WLEN", decomp_id, grid_h2d_id, &buf_mark, &usage_tag, "m", "register field instance of Wave Length");
	field_id_UBOT = register_model_float_1D_data_wrj(UBOT_from_swan, meshunion->mesh2d_p->Nv2d, "UBOT", decomp_id, grid_h2d_id, &buf_mark, &usage_tag, "m/s", "register field instance of Bottom wave point velocity");
	field_id_TMBOT = register_model_float_1D_data_wrj(TMBOT_from_swan, meshunion->mesh2d_p->Nv2d, "TMBOT", decomp_id, grid_h2d_id, &buf_mark, &usage_tag, "s", "register field instance of Bottom wave point period");

    field_id_H = register_model_float_1D_data_wrj(H_to_swan, meshunion->mesh2d_p->Nv2d, "H", decomp_id, grid_h2d_id, &buf_mark, &usage_tag, "m", "register field instance of Water Level");                                   
    field_id_U = register_model_float_1D_data_wrj(U_to_swan, meshunion->mesh2d_p->Nv2d, "U", decomp_id, grid_h2d_id, &buf_mark, &usage_tag, "m/s", "register field instance of The Velocity in the x direction");
    field_id_V = register_model_float_1D_data_wrj(V_to_swan, meshunion->mesh2d_p->Nv2d, "V", decomp_id, grid_h2d_id, &buf_mark, &usage_tag, "m/s", "register field instance of The Velocity in the y direction"); 
    field_id_test = register_model_float_1D_data_wrj(test_to_swan, meshunion->mesh2d_p->Nv2d, "test", decomp_id, grid_h2d_id, &buf_mark, &usage_tag, "1", "register field instance of tset"); 
// --------register coupling frequency to C-Coupler2-------------
    timer_id = define_single_timer_wrj(dg_demo_comp_id, "seconds", coupling_freq, &local_lag_count, &remote_lag_count,"define a single timer for swan_demo");

// --------register export & import interface to C-Coupler2------
     fields_id_s[0] = field_id_H;
     fields_id_s[1] = field_id_U;
     fields_id_s[2] = field_id_V;
     fields_id_s[3] = field_id_test;
    int num_field_instances=4;

    int export_interface_id, import_interface_id;
    export_interface_id = register_export_interface_wrj("send_data_to_swan", &num_field_instances, fields_id_s, &timer_id, "register interface for sending data to swan");



    fields_id_r[0] = field_id_HS;
    fields_id_r[1] = field_id_DIR;
    fields_id_r[2] = field_id_RTP;
    fields_id_r[3] = field_id_QB;
	fields_id_r[4] = field_id_WLEN;
	fields_id_r[5] = field_id_UBOT;
	fields_id_r[6] = field_id_TMBOT;
    num_field_instances = 7;
    int inst_or_aver = 0;
    
    import_interface_id = register_import_interface_wrj("receive_data_from_swan", &num_field_instances, fields_id_r, &timer_id, &inst_or_aver, "register interface for receiving data from swan");

    std::cout<<"before end_coupling_configuration_wrj  "<< *dg_demo_comp_id <<std::endl;
    end_coupling_configuration_wrj(dg_demo_comp_id, "component atm_demo ends configuration");
    
    std::cout<<"after end_coupling_configuration_wrj"<<std::endl;

	//	std::cout<<"before end_coupling_configuration_wrj"<<std::endl;
}

// bool *bypass_timer = new bool;
// *bypass_timer = false;
// bool interface_status;
// interface_status = execute_interface_using_name_wrj(dg_demo_comp_id,"send_data_to_swan", bypass_timer, "execute interface for sending data to swan");
// interface_status = execute_interface_using_name_wrj(dg_demo_comp_id,"receive_data_to_swan", bypass_timer, "execute interface for receiving data from swan");
// advance_time_wrj(dg_demo_comp_id, "atm_demo advances time for one step");


// (int* component_id, const char* interface_name, bool* bypass_timer, const char* annotation)



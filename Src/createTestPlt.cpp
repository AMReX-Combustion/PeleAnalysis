#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <string>
#include <map>

using namespace amrex;

// Enum for field types
enum class FieldType {
    Constant,
    PlaneStep,
    PlaneSmooth,
    DoublePlaneStep,
    DoublePlaneSmooth,
    CircleStep,
    CircleSmooth,
    RingStep,
    RingSmooth,
    Sine
};

// Structure to hold field configuration
struct FieldConfig {
    char name[16];
    FieldType type;
    
    // Common parameters
    Real value_inside = 1.0;
    Real value_outside = 0.0;
    
    // Plane parameters
    int plane_axis = 0;      // 0=x, 1=y, 2=z
    Real plane_position = 0.5;
    Real plane_position_second = 0.6;
    Real smooth_width = 0.1;
    
    // Circle/Sphere parameters
    GpuArray<Real, AMREX_SPACEDIM> center = {AMREX_D_DECL(0.5, 0.5, 0.5)};
    Real radius = 0.25;

    // Ring/Spherical shell parameters
    Real radius_inner = 0.25;
    Real radius_outer = 0.25;
    
    // Sine parameters
    GpuArray<Real, AMREX_SPACEDIM> frequency = {AMREX_D_DECL(1.0, 1.0, 1.0)};
    GpuArray<Real, AMREX_SPACEDIM> phase = {AMREX_D_DECL(0.0, 0.0, 0.0)};
    Real amplitude = 1.0;
    Real offset = 0.0;
};

// Helper function to parse field type
FieldType parseFieldType(const std::string& type_str) {
    static const std::map<std::string, FieldType> type_map = {
        {"constant", FieldType::Constant},
        {"plane_step", FieldType::PlaneStep},
        {"plane_smooth", FieldType::PlaneSmooth},
        {"double_plane_step", FieldType::DoublePlaneStep},
        {"double_plane_smooth", FieldType::DoublePlaneSmooth},
        {"circle_step", FieldType::CircleStep},
        {"circle_smooth", FieldType::CircleSmooth},
        {"sphere_step", FieldType::CircleStep},    // Same as circle in code
        {"sphere_smooth", FieldType::CircleSmooth},
        {"ring_step", FieldType::RingStep},
        {"ring_smooth", FieldType::RingSmooth},
        {"spherical_shell_step", FieldType::RingStep},    // Same as ring in code
        {"spherical_shell_smooth", FieldType::RingSmooth},
        {"sine", FieldType::Sine}
    };
    
    auto it = type_map.find(type_str);
    if (it != type_map.end()) {
        return it->second;
    }
    
    amrex::Abort("Unknown field type: " + type_str);
    return FieldType::Constant;
}

// Smooth step function
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real smoothstep(Real edge0, Real edge1, Real x) {
    Real t = std::max(0.0, std::min(1.0, (x - edge0) / (edge1 - edge0)));
    
    // smoothstep
    //return t * t * (3.0 - 2.0 * t);

    // smootherstep
    return t * t * t * (t * (6.0 * t - 15.0) + 10.0);

    // tanh
    //return 0.5 * (1 + std::tanh(6 * t - 3));
}

// Function to evaluate field value at a point
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real evaluateField(const FieldConfig& config, 
                   Real x, Real y, Real z,
                   const GpuArray<Real, AMREX_SPACEDIM>& dx,
                   const GpuArray<Real, AMREX_SPACEDIM>& prob_lo) {
    
    Real result = 0.0;
    
    switch (config.type) {
        case FieldType::Constant:
            result = config.value_inside;
            break;
            
        case FieldType::PlaneStep: {
            Real coord;
            if (config.plane_axis == 0) coord = x;
            else if (config.plane_axis == 1) coord = y;
            else coord = z;
            
            result = (coord < config.plane_position) ? 
                     config.value_inside : config.value_outside;
            break;
        }
        
        case FieldType::PlaneSmooth: {
            Real coord;
            if (config.plane_axis == 0) coord = x;
            else if (config.plane_axis == 1) coord = y;
            else coord = z;
            
            Real edge0 = config.plane_position - config.smooth_width / 2.0;
            Real edge1 = config.plane_position + config.smooth_width / 2.0;
            Real blend = smoothstep(edge0, edge1, coord);
            result = config.value_inside * (1.0 - blend) + config.value_outside * blend;
            break;
        }

        case FieldType::DoublePlaneStep: {
            Real coord;
            if (config.plane_axis == 0) coord = x;
            else if (config.plane_axis == 1) coord = y;
            else coord = z;
            
            if ((coord < config.plane_position) || (coord > config.plane_position_second)) {
                result = config.value_outside;
            } else {
                result = config.value_inside;
            }
            break;
        }
        
        case FieldType::DoublePlaneSmooth: {
            Real coord;
            if (config.plane_axis == 0) coord = x;
            else if (config.plane_axis == 1) coord = y;
            else coord = z;
            
            // First plane transition
            Real edge0_first = config.plane_position - config.smooth_width / 2.0;
            Real edge1_first = config.plane_position + config.smooth_width / 2.0;
            Real blend_first = smoothstep(edge0_first, edge1_first, coord);
            
            // Second plane transition
            Real edge0_second = config.plane_position_second - config.smooth_width / 2.0;
            Real edge1_second = config.plane_position_second + config.smooth_width / 2.0;
            Real blend_second = smoothstep(edge0_second, edge1_second, coord);
            
            // Combine: outside if coord < first OR coord > second, inside if first < coord < second
            // Assumes plane_position < plane_position_second
            // blend_first: 0 when coord < first, 1 when coord > first
            // blend_second: 0 when coord < second, 1 when coord > second
            Real blend = blend_first * (1.0 - blend_second);
            
            result = config.value_inside * blend + config.value_outside * (1.0 - blend);
            break;
        }
        
        case FieldType::CircleStep: {
#if AMREX_SPACEDIM == 2
            Real dx_val = x - config.center[0];
            Real dy_val = y - config.center[1];
            Real r = std::sqrt(dx_val * dx_val + dy_val * dy_val);
#else
            Real dx_val = x - config.center[0];
            Real dy_val = y - config.center[1];
            Real dz_val = z - config.center[2];
            Real r = std::sqrt(dx_val * dx_val + dy_val * dy_val + dz_val * dz_val);
#endif
            result = (r < config.radius) ? 
                     config.value_inside : config.value_outside;
            break;
        }
        
        case FieldType::CircleSmooth: {
#if AMREX_SPACEDIM == 2
            Real dx_val = x - config.center[0];
            Real dy_val = y - config.center[1];
            Real r = std::sqrt(dx_val * dx_val + dy_val * dy_val);
#else
            Real dx_val = x - config.center[0];
            Real dy_val = y - config.center[1];
            Real dz_val = z - config.center[2];
            Real r = std::sqrt(dx_val * dx_val + dy_val * dy_val + dz_val * dz_val);
#endif
            Real edge0 = config.radius - config.smooth_width / 2.0;
            Real edge1 = config.radius + config.smooth_width / 2.0;
            Real blend = smoothstep(edge0, edge1, r);
            result = config.value_inside * (1.0 - blend) + config.value_outside * blend;
            break;
        }

        case FieldType::RingStep: {
#if AMREX_SPACEDIM == 2
            Real dx_val = x - config.center[0];
            Real dy_val = y - config.center[1];
            Real r = std::sqrt(dx_val * dx_val + dy_val * dy_val);
#else
            Real dx_val = x - config.center[0];
            Real dy_val = y - config.center[1];
            Real dz_val = z - config.center[2];
            Real r = std::sqrt(dx_val * dx_val + dy_val * dy_val + dz_val * dz_val);
#endif

            if ((r < config.radius_inner) || (r > config.radius_outer)) {
                result = config.value_outside;
            } else {
                result = config.value_inside;
            }
            break;
        }
        
        case FieldType::RingSmooth: {
#if AMREX_SPACEDIM == 2
            Real dx_val = x - config.center[0];
            Real dy_val = y - config.center[1];
            Real r = std::sqrt(dx_val * dx_val + dy_val * dy_val);
#else
            Real dx_val = x - config.center[0];
            Real dy_val = y - config.center[1];
            Real dz_val = z - config.center[2];
            Real r = std::sqrt(dx_val * dx_val + dy_val * dy_val + dz_val * dz_val);
#endif
            // Smooth transition at inner radius (outside -> inside)
            Real edge0_inner = config.radius_inner - config.smooth_width / 2.0;
            Real edge1_inner = config.radius_inner + config.smooth_width / 2.0;
            Real blend_inner = smoothstep(edge0_inner, edge1_inner, r);
            
            // Smooth transition at outer radius (inside -> outside)
            Real edge0_outer = config.radius_outer - config.smooth_width / 2.0;
            Real edge1_outer = config.radius_outer + config.smooth_width / 2.0;
            Real blend_outer = smoothstep(edge0_outer, edge1_outer, r);
            
            // Combine: outside if r < inner OR r > outer, inside if inner < r < outer
            // blend_inner: 0 when r < inner, 1 when r > inner
            // blend_outer: 0 when r < outer, 1 when r > outer
            Real blend = blend_inner * (1.0 - blend_outer);
            
            result = config.value_inside * blend + config.value_outside * (1.0 - blend);
            break;
        }

        
        case FieldType::Sine: {
            Real arg_x = 2.0 * M_PI * config.frequency[0] * x + config.phase[0];
            Real arg_y = 2.0 * M_PI * config.frequency[1] * y + config.phase[1];

#if AMREX_SPACEDIM == 3
            Real arg_z = 2.0 * M_PI * config.frequency[2] * z + config.phase[2];
#else
            Real arg_z = 0.0;
#endif
            result = config.offset + config.amplitude * 
                     std::sin(arg_x) * std::cos(arg_y) * std::sin(arg_z);
            break;
        }
    }
    
    return result;
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        // Read parameters
        ParmParse pp;
        
        // Geometry parms
        Vector<int> n_cell(AMREX_SPACEDIM);
        pp.getarr("amr.n_cell", n_cell, 0, AMREX_SPACEDIM);
        
        // Domain parms
        RealBox real_box;
        Vector<Real> pp_prob_x(AMREX_SPACEDIM, 0.0);
        pp.getarr("geometry.prob_lo", pp_prob_x, 0, AMREX_SPACEDIM);
        GpuArray<Real,AMREX_SPACEDIM> prob_lo = {AMREX_D_DECL(pp_prob_x[0],pp_prob_x[1],pp_prob_x[2])};
        pp.getarr("geometry.prob_hi", pp_prob_x, 0, AMREX_SPACEDIM);
        GpuArray<Real,AMREX_SPACEDIM> prob_hi = {AMREX_D_DECL(pp_prob_x[0],pp_prob_x[1],pp_prob_x[2])};
        
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            real_box.setLo(i, prob_lo[i]);
            real_box.setHi(i, prob_hi[i]);
        }
        
        // Periodicity
        IntVect pp_is_per;
        pp.getarr("geometry.is_periodic",pp_is_per);
        Array<int,AMREX_SPACEDIM> is_per = {AMREX_D_DECL(pp_is_per[0],pp_is_per[1],pp_is_per[2])};
        
        // Coordinate system
        int coord = 0;
        pp.query("geometry.coord_sys", coord);

        // Only cartesian supported
        AMREX_ALWAYS_ASSERT(coord==0);
        
        // Create geometry
        IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
        IntVect domain_hi(AMREX_D_DECL(n_cell[0]-1, 
                                        n_cell[1]-1, 
                                        n_cell[2]-1));
        Box domain(domain_lo, domain_hi);
        
        Geometry geom(domain, real_box, coord, is_per);
        
        // Get cell size
        const auto dx = geom.CellSizeArray();
        
        // Create BoxArray and DistributionMapping 
        BoxArray ba(domain);
        
        // Distribute domain
        int max_grid_size = 32;
        pp.query("amr.max_grid_size", max_grid_size);
        ba.maxSize(max_grid_size);
        int ngrow = 0; // Ghost cells
        pp.query("amr.ngrow", ngrow);
        
        DistributionMapping dm(ba);
        
        // Parse field configurations
        int ncomp = pp.countval("field.names");
        Vector<std::string> field_names(ncomp);
        Vector<FieldConfig> field_configs(ncomp);
        
        if (ncomp > 0) {
            pp.getarr("field.names", field_names);
            Print() << "Fields:" << std::endl;
            
            for (int i = 0; i < ncomp; ++i) {
                ParmParse ppf(field_names[i]);
                FieldConfig& config = field_configs[i];
                std::strncpy(config.name, field_names[i].c_str(), 15);
                config.name[15] = '\0';
                
                std::string type_str;
                ppf.get("type", type_str);
                config.type = parseFieldType(type_str);
                
                Print() << "  Field " << i << ": " << field_names[i] 
                        << "  (type: " << type_str << ")" << std::endl;
                
                // Read parameters based on field type
                switch (config.type) {
                    case FieldType::Constant:{
                        ppf.query("value", config.value_inside);
                        break;
                    }
                        
                    case FieldType::PlaneStep:
                    case FieldType::PlaneSmooth:{
                        ppf.query("axis", config.plane_axis);
                        ppf.query("position", config.plane_position);
                        ppf.query("value_left", config.value_inside);
                        ppf.query("value_right", config.value_outside);
                        if (config.type == FieldType::PlaneSmooth) {
                            ppf.query("smooth_width", config.smooth_width);
                        }
                        break;
                    }

                    case FieldType::DoublePlaneStep:
                    case FieldType::DoublePlaneSmooth:{
                        ppf.query("axis", config.plane_axis);
                        ppf.query("position", config.plane_position);
                        ppf.query("position_second", config.plane_position_second);
                        ppf.query("value_inside", config.value_inside);
                        ppf.query("value_outside", config.value_outside);
                        if (config.type == FieldType::DoublePlaneSmooth) {
                            ppf.query("smooth_width", config.smooth_width);
                        }
                        break;
                    }
                        
                    case FieldType::CircleStep:
                    case FieldType::CircleSmooth: {
                        Vector<Real> pp_center(AMREX_SPACEDIM, 0.0);
                        ppf.queryarr("center", pp_center);
                        for (int dim=0; dim<pp_center.size(); ++dim){
                            config.center[dim] = pp_center[dim];
                        }
                        ppf.query("radius", config.radius);
                        ppf.query("value_inside", config.value_inside);
                        ppf.query("value_outside", config.value_outside);
                        if (config.type == FieldType::CircleSmooth) {
                            ppf.query("smooth_width", config.smooth_width);
                        }
                        break;
                    }

                    case FieldType::RingStep:
                    case FieldType::RingSmooth: {
                        Vector<Real> pp_center(AMREX_SPACEDIM, 0.0);
                        ppf.queryarr("center", pp_center);
                        for (int dim=0; dim<pp_center.size(); ++dim){
                            config.center[dim] = pp_center[dim];
                        }
                        ppf.query("radius_inner", config.radius_inner);
                        ppf.query("radius_outer", config.radius_outer);
                        ppf.query("value_inside", config.value_inside);
                        ppf.query("value_outside", config.value_outside);
                        if (config.type == FieldType::RingSmooth) {
                            ppf.query("smooth_width", config.smooth_width);
                        }
                        break;
                    }
                        
                    case FieldType::Sine:{
                        Vector<Real> pp_frequency(AMREX_SPACEDIM, 0.0);
                        ppf.queryarr("frequency", pp_frequency);
                        for (int dim=0; dim<pp_frequency.size(); ++dim){
                            config.frequency[dim] = pp_frequency[dim];
                        }
                        Vector<Real> pp_phase(AMREX_SPACEDIM, 0.0);
                        ppf.queryarr("phase", pp_phase);
                        for (int dim=0; dim<pp_phase.size(); ++dim){
                            config.phase[dim] = pp_phase[dim];
                        }                        
                        ppf.query("amplitude", config.amplitude);
                        ppf.query("offset", config.offset);
                        break;
                    }
                }
            }
        } else {
            Abort("No fields specified...");
        }
        
        // Create MultiFab
        MultiFab mf(ba, dm, ncomp, ngrow);
        
        // Fill MultiFab with field data
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const Box& box = mfi.validbox();
            auto const& fab = mf.array(mfi);
            
            // Copy configs to device-accessible memory
            Gpu::DeviceVector<FieldConfig> d_configs(ncomp);
            Gpu::copyAsync(Gpu::hostToDevice, field_configs.begin(), 
                          field_configs.end(), d_configs.begin());
            Gpu::streamSynchronize();
            
            FieldConfig* p_configs = d_configs.data();
            
            amrex::ParallelFor(box, ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                // Get physical coordinates at cell center
                Real x = prob_lo[0] + (i + 0.5) * dx[0];
                Real y = prob_lo[1] + (j + 0.5) * dx[1];
#if AMREX_SPACEDIM == 3
                Real z = prob_lo[2] + (k + 0.5) * dx[2];
#else
                Real z = 0.0;
#endif
                
                fab(i,j,k,n) = evaluateField(p_configs[n], x, y, z, dx, prob_lo);
            });
        }
        
        // Variable names for plotfile
        Vector<std::string> varnames(ncomp);
        for (int n = 0; n < ncomp; ++n) {
            varnames[n] = field_names[n];
        }
        
        // Write plot file
        std::string plotfile_name = "pltTestFile";
        pp.query("plotfile_name", plotfile_name);
        
        WriteSingleLevelPlotfile(plotfile_name, mf, varnames, geom, 0.0, 0);
        
        amrex::Print() << "Plotfile created: " << plotfile_name << "\n";
    }
    amrex::Finalize();
    return 0;
}
#include "Adaptor_legacy.h"

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkStringArray.h>

namespace {
vtkCPProcessor *Processor = nullptr;
vtkImageData *VTKGrid = nullptr;
const char *InputName = "input";

int _start_x;
int _nx;
int _dx;
int _start_y;
int _ny;
int _dy;
int _start_z;
int _nz;
int _dz;
const Collective *_sim_params{};

//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
void UpdateVTKAttributes(vtkCPInputDataDescription *idd, EMfields3D *EMf) {
  // I am not sure whether we need to do this check
  if (idd->IsFieldNeeded("B", vtkDataObject::POINT) == true) {
    // Create a VTK object representing magnetic field array

    // Get a reference to the grid's point data object.
    vtkPointData *vtk_point_data = VTKGrid->GetPointData();

    // We need to create a new VTK array object and attach it to the point data,
    // if it hasn't been done yet.

    int ns = _sim_params->getNs();
    const string rhotag[]={"rhoe0", "rhoi1", "rhoe2", "rhoi3", "rhoe4", "rhoi5", "rhoe6", "rhoi7"};
    const string Vtag[]={"Ve0", "Vi1", "Ve2", "Vi3", "Ve4", "Vi5", "Ve6", "Vi7"};
    const string Tcart_tag[]={"Tcart_e0", "Tcart_i1", "Tcart_e2", "Tcart_i3", "Tcart_e4", "Tcart_i5", "Tcart_e6", "Tcart_i7"};
    const string Tperpar_tag[]={"Tperpar_e0", "Tperpar_i1", "Tperpar_e2", "Tperpar_i3", "Tperpar_e4", "Tperpar_i5", "Tperpar_e6", "Tperpar_i7"};


    if (vtk_point_data->GetNumberOfArrays() == 0) {
      //B
      vtkNew<vtkDoubleArray> field_array_B;
      field_array_B->SetName("B");
      field_array_B->SetNumberOfComponents(3);
      field_array_B->SetNumberOfTuples(VTKGrid->GetNumberOfPoints());
      vtk_point_data->AddArray(field_array_B);

      //E
      vtkNew<vtkDoubleArray> field_array_E;
      field_array_E->SetName("E");
      field_array_E->SetNumberOfComponents(3);
      field_array_E->SetNumberOfTuples(VTKGrid->GetNumberOfPoints());
      vtk_point_data->AddArray(field_array_E);

      //add a VTKdoubleArray for each species (everycycle overwrite the vtkDoubleArray)  
      for(int si=0; si<ns; si++){
        //rho
        vtkNew<vtkDoubleArray> rho_array;
        rho_array->SetName(rhotag[si].c_str());
        rho_array->SetNumberOfComponents(1);
        rho_array->SetNumberOfTuples(VTKGrid->GetNumberOfPoints());
        vtk_point_data->AddArray(rho_array);
        //V
        vtkNew<vtkDoubleArray> field_array_V;
        field_array_V->SetName(Vtag[si].c_str());
        field_array_V->SetNumberOfComponents(3);
        field_array_V->SetNumberOfTuples(VTKGrid->GetNumberOfPoints());
        vtk_point_data->AddArray(field_array_V);
        //T_cart
        vtkNew<vtkDoubleArray> field_array_Tcart;
        field_array_Tcart->SetName(Tcart_tag[si].c_str());
        field_array_Tcart->SetNumberOfComponents(6);
        field_array_Tcart->SetNumberOfTuples(VTKGrid->GetNumberOfPoints());
        vtk_point_data->AddArray(field_array_Tcart);
        //T_perpar
        vtkNew<vtkDoubleArray> field_array_Tperpar;
        field_array_Tperpar->SetName(Tperpar_tag[si].c_str());
        field_array_Tperpar->SetNumberOfComponents(6);
        field_array_Tperpar->SetNumberOfTuples(VTKGrid->GetNumberOfPoints());
        vtk_point_data->AddArray(field_array_Tperpar);
      }      
    }
    
    vtkDoubleArray *field_array_B =
        vtkDoubleArray::SafeDownCast(vtk_point_data->GetArray("B"));
    
    vtkDoubleArray *field_array_E =
        vtkDoubleArray::SafeDownCast(vtk_point_data->GetArray("E"));

    vtkDoubleArray *rhons[ns];
    vtkDoubleArray *Vns[ns];
    vtkDoubleArray *Tcart_ns[ns];
    vtkDoubleArray *Tperpar_ns[ns];
    for(int si = 0; si<ns; ++si){
      rhons[si]         = vtkDoubleArray::SafeDownCast(vtk_point_data->GetArray(rhotag[si].c_str()));
      Vns[si]           = vtkDoubleArray::SafeDownCast(vtk_point_data->GetArray(Vtag[si].c_str()));
      Tcart_ns[si]      = vtkDoubleArray::SafeDownCast(vtk_point_data->GetArray(Tcart_tag[si].c_str()));
      Tperpar_ns[si]    = vtkDoubleArray::SafeDownCast(vtk_point_data->GetArray(Tperpar_tag[si].c_str()));
    }
    
    

    // Feed the data into VTK array. Since we don't know the memory layout of
    // our B field data, we feed it point-by-point, in a very slow way

    // Array of grid's dimensions
    int *dims = VTKGrid->GetDimensions();

    auto Bx = EMf->getBxTot();
    auto By = EMf->getByTot();
    auto Bz = EMf->getBzTot();

    auto Ex = EMf->getEx();
    auto Ey = EMf->getEy();
    auto Ez = EMf->getEz();

    auto Jxs = EMf->getJxs();
    auto Jys = EMf->getJys();
    auto Jzs = EMf->getJzs();

    auto rho = EMf->getRHOns();




    // Cycle over all VTK grid's points, get their indices and copy the data.
    // We want to have only one cycle over point's ID to efficiently use
    // multi-threading.
    for(int si=0; si<ns; si++){
      EMf->calcT_si(si);

      auto Tcart = EMf -> getTcart();
      auto Tperpar = EMf -> getTperpar();

      for (vtkIdType p = 0; p < VTKGrid->GetNumberOfPoints(); ++p) {
        // Get cells's indices i, j , k
        const size_t k = p / (dims[0] * dims[1]);
        const size_t j = (p - k * dims[0] * dims[1]) / dims[0];
        const size_t i = p - k * dims[0] * dims[1] - j * dims[0];

        //not repeating for each species
        if(si==0){
          // CAUTION!!! K should be always zero in the 2D case?
          field_array_B->SetComponent(p, 0, Bx[i+1][j+1][k+1]);
          field_array_B->SetComponent(p, 1, By[i+1][j+1][k+1]);
          field_array_B->SetComponent(p, 2, Bz[i+1][j+1][k+1]);

          field_array_E->SetComponent(p, 0, Ex[i+1][j+1][k+1]);
          field_array_E->SetComponent(p, 1, Ey[i+1][j+1][k+1]);
          field_array_E->SetComponent(p, 2, Ez[i+1][j+1][k+1]);
        }
        //printf("settig coordinates: %d %d %d spicies: %d\n", i,j,k,si);
       
        rhons[si]-> SetValue(p, rho[si][i+1][j+1][k+1]*4*3.1415926535897);

        Vns[si]  -> SetComponent(p, 0, Jxs[si][i+1][j+1][k+1]/rho[si][i+1][j+1][k+1]*4*3.1415926535897);
        Vns[si]  -> SetComponent(p, 1, Jys[si][i+1][j+1][k+1]/rho[si][i+1][j+1][k+1]*4*3.1415926535897);
        Vns[si]  -> SetComponent(p, 2, Jzs[si][i+1][j+1][k+1]/rho[si][i+1][j+1][k+1]*4*3.1415926535897);

        Tcart_ns[si]  -> SetComponent(p,0, Tcart[0][i+1][j+1][k+1]);
        Tcart_ns[si]  -> SetComponent(p,1, Tcart[1][i+1][j+1][k+1]);
        Tcart_ns[si]  -> SetComponent(p,2, Tcart[2][i+1][j+1][k+1]);
        Tcart_ns[si]  -> SetComponent(p,3, Tcart[3][i+1][j+1][k+1]);
        Tcart_ns[si]  -> SetComponent(p,4, Tcart[4][i+1][j+1][k+1]);
        Tcart_ns[si]  -> SetComponent(p,5, Tcart[5][i+1][j+1][k+1]);

        Tperpar_ns[si]  -> SetComponent(p,0, Tperpar[0][i+1][j+1][k+1]);
        Tperpar_ns[si]  -> SetComponent(p,1, Tperpar[1][i+1][j+1][k+1]);
        Tperpar_ns[si]  -> SetComponent(p,2, Tperpar[2][i+1][j+1][k+1]);
        Tperpar_ns[si]  -> SetComponent(p,3, Tperpar[3][i+1][j+1][k+1]);
        Tperpar_ns[si]  -> SetComponent(p,4, Tperpar[4][i+1][j+1][k+1]);
        Tperpar_ns[si]  -> SetComponent(p,5, Tperpar[5][i+1][j+1][k+1]);

        //if (si == 0 & i== 3 & j == 3 & k == 3) printf("adaptor Tcart: %f\n",Tperpar[0][i][j][k]);
        //if (Tcart[0][i+1][j+1][k+1] != 0 ) printf("coordinates: %d %d %d spicies: %d \t temp: %f\n", i,j,k,si,Tcart[0][i+1][j+1][k+1]);
      }
    }

    /// Fast way, if memry layout is correct.
    // velocityData->SetArray(const_cast<double*>(velocity.data()),
    // static_cast<vtkIdType>(velocity.size()), 1);
  }
  //  if (idd->IsFieldNeeded("collision", vtkDataObject::POINT) == true)
  //  {
  //    if (VTKGrid->GetPointData()->GetArray("collision") == nullptr)
  //    {
  //      // velocity array
  //      vtkNew<vtkIntArray> collisionData;
  //      collisionData->SetName("collision");
  //      collisionData->SetNumberOfComponents(1);
  //      collisionData->SetNumberOfTuples(static_cast<vtkIdType>(collisions.size()));
  //      VTKGrid->GetPointData()->AddArray(collisionData);
  //    }
  //    vtkIntArray* collisionData =
  //      vtkIntArray::SafeDownCast(VTKGrid->GetPointData()->GetArray("collision"));
  //
  //    collisionData->SetArray(const_cast<int*>(collisions.data()),
  //    static_cast<vtkIdType>(collisions.size()), 1);
  //  }
}

//----------------------------------------------------------------------------
void BuildVTKDataStructures(vtkCPInputDataDescription *idd, EMfields3D *EMf) {
  // feed data to grid
  UpdateVTKAttributes(idd, EMf);
}
} // namespace

namespace Adaptor_legacy {

//----------------------------------------------------------------------------
void Initialize(const Collective *sim_params, const int start_x,
                const int start_y, const int start_z, const int nx,
                const int ny, const int nz, const double dx, const double dy,
                const double dz) {
  if (Processor == NULL) {
    Processor = vtkCPProcessor::New();
    Processor->Initialize();
  } else {
    Processor->RemoveAllPipelines();
  }
  vtkNew<vtkCPPythonScriptPipeline> pipeline;
  pipeline->Initialize(sim_params->getParaviewScriptPath().c_str());
  Processor->AddPipeline(pipeline);

  _start_x = start_x;
  _nx = nx;
  _dx = dx;
  _start_y = start_y;
  _ny = ny;
  _dy = dy;
  _start_z = start_z;
  _nz = nz;
  _dz = dz;

  _sim_params = sim_params;

  if (VTKGrid == NULL) {
    // The grid structure isn't changing so we only build it
    // the first time it's needed. If we needed the memory
    // we could delete it and rebuild as necessary.
    VTKGrid = vtkImageData::New();
    printf("start: (%d, %d, %d); end: (%d,%d,%d)\n", start_x, start_y, start_z, start_x + nx- 3, start_y + ny- 3, start_z + nz - 3);
    VTKGrid->SetExtent(start_x, start_x + nx - 3, start_y , start_y + ny - 3,
                       start_z, start_z + nz - 3);
    VTKGrid->SetSpacing(dx, dy, dz);
  }
}

//----------------------------------------------------------------------------
void Finalize() {
  if (Processor) {
    Processor->Delete();
    Processor = NULL;
  }
  if (VTKGrid) {
    VTKGrid->Delete();
    VTKGrid = NULL;
  }
}

//----------------------------------------------------------------------------
void CoProcess(double time, unsigned int timeStep, EMfields3D *EMf) {
  vtkNew<vtkCPDataDescription> dataDescription;
  dataDescription->AddInput(InputName);
  dataDescription->SetTimeData(time, timeStep);

  vtkNew<vtkStringArray> fd0{};
  fd0->SetName("CaseName");
  fd0->SetNumberOfComponents(1);
  fd0->InsertNextValue(_sim_params->getCase().c_str());
  VTKGrid->GetFieldData()->AddArray(fd0);

  vtkNew<vtkIntArray> fd1{};
  fd1->SetName("TimeStep");
  fd1->SetNumberOfComponents(1);
  fd1->InsertNextValue(timeStep);
  VTKGrid->GetFieldData()->AddArray(fd1);

  std::vector<std::pair<std::string, double>> params{
      {"B0x", _sim_params->getB0x()},
      {"B0y", _sim_params->getB0y()},
      {"B0z", _sim_params->getB0z()},
      {"ns", _sim_params->getNs()}};

  for (const auto &pair : params) {
    vtkNew<vtkDoubleArray> fd{};
    fd->SetNumberOfComponents(1);
    fd->SetName(pair.first.c_str());
    fd->InsertNextValue(pair.second);
    VTKGrid->GetFieldData()->AddArray(fd);
  }

 if (Processor->RequestDataDescription(dataDescription) != 0) {

    vtkCPInputDataDescription *idd =
        dataDescription->GetInputDescriptionByName(InputName);
    BuildVTKDataStructures(idd, EMf);
    idd->SetGrid(VTKGrid);
    Processor->CoProcess(dataDescription);
  }
}
} // namespace Adaptor

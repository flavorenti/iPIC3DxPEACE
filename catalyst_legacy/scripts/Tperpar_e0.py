
#--------------------------------------------------------------

# Global timestep output options
timeStepToStartOutputAt=0
forceOutputAtFirstCall=False

# Global screenshot output options
imageFileNamePadding=3
rescale_lookuptable=False

# Whether or not to request specific arrays from the adaptor.
requestSpecificArrays=False

# a root directory under which all Catalyst output goes
rootDirectory='/home/lquerci/NICE/build_iPIC3D/data/images'

# makes a cinema D index table
make_cinema_table=False

#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# paraview version 5.7.0
#--------------------------------------------------------------

from paraview.simple import *
from paraview import coprocessing

# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.7.0

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      # trace generated using paraview version 5.7.0
      #
      # To ensure correct image size when batch processing, please search 
      # for and uncomment the line `# renderView*.ViewSize = [*,*]`

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # get the material library
      materialLibrary1 = GetMaterialLibrary()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1216, 729]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [3.9999959468841553, 3.9999959468841553, 0.9999974966049194]
      renderView1.StereoType = 'Crystal Eyes'
      renderView1.CameraPosition = [2.2629547077843464, -0.9567793057646335, 22.56484691861354]
      renderView1.CameraFocalPoint = [3.9999959468841566, 3.9999959468841575, 0.99999749660492]
      renderView1.CameraViewUp = [0.003574389112650775, 0.9745170552558612, 0.2242849365381373]
      renderView1.CameraFocalDisk = 1.0
      renderView1.CameraParallelScale = 5.744556566297824
      renderView1.Background = [0.32, 0.34, 0.43]
      renderView1.BackEnd = 'OSPRay raycaster'
      renderView1.OSPRayMaterialLibrary = materialLibrary1

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='image_%t.png', freq=1, fittoscreen=0, magnification=1, width=1216, height=729, cinema={}, compression=-1)
      renderView1.ViewTime = datadescription.GetTime()

      SetActiveView(None)

      # ----------------------------------------------------------------
      # setup view layouts
      # ----------------------------------------------------------------

      # create new layout object 'Layout #1'
      layout1 = CreateLayout(name='Layout #1')
      layout1.AssignView(0, renderView1)

      # ----------------------------------------------------------------
      # restore active view
      SetActiveView(renderView1)
      # ----------------------------------------------------------------

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'Legacy VTK Reader'
      # create a producer from a simulation input
      exosphere_Tperpare0_0vtk = coprocessor.CreateProducer(datadescription, 'input')

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from exosphere_Tperpare0_0vtk
      exosphere_Tperpare0_0vtkDisplay = Show(exosphere_Tperpare0_0vtk, renderView1)

      # get color transfer function/color map for 'Tperpar_e0'
      tperpar_eLUT = GetColorTransferFunction('Tperpar_e0')
      tperpar_eLUT.RGBPoints = [-4.969829669326797e-12, 0.231373, 0.298039, 0.752941, 1.167813867851441e-05, 0.865003, 0.865003, 0.865003, 2.335628232685849e-05, 0.705882, 0.0156863, 0.14902]
      tperpar_eLUT.ScalarRangeInitialized = 1.0
      tperpar_eLUT.VectorMode = 'Component'

      # get opacity transfer function/opacity map for 'Tperpar_e0'
      tperpar_ePWF = GetOpacityTransferFunction('Tperpar_e0')
      tperpar_ePWF.Points = [-4.969829669326797e-12, 0.0, 0.5, 0.0, 2.335628232685849e-05, 1.0, 0.5, 0.0]
      tperpar_ePWF.ScalarRangeInitialized = 1

      # trace defaults for the display properties.
      exosphere_Tperpare0_0vtkDisplay.Representation = 'Slice'
      exosphere_Tperpare0_0vtkDisplay.ColorArrayName = ['POINTS', 'Tperpar_e0']
      exosphere_Tperpare0_0vtkDisplay.LookupTable = tperpar_eLUT
      exosphere_Tperpare0_0vtkDisplay.OSPRayScaleArray = 'Tperpar_e0'
      exosphere_Tperpare0_0vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      exosphere_Tperpare0_0vtkDisplay.SelectOrientationVectors = 'None'
      exosphere_Tperpare0_0vtkDisplay.ScaleFactor = 0.7999992000000001
      exosphere_Tperpare0_0vtkDisplay.SelectScaleArray = 'None'
      exosphere_Tperpare0_0vtkDisplay.GlyphType = 'Arrow'
      exosphere_Tperpare0_0vtkDisplay.GlyphTableIndexArray = 'None'
      exosphere_Tperpare0_0vtkDisplay.GaussianRadius = 0.03999996
      exosphere_Tperpare0_0vtkDisplay.SetScaleArray = ['POINTS', 'Tperpar_e0']
      exosphere_Tperpare0_0vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      exosphere_Tperpare0_0vtkDisplay.OpacityArray = ['POINTS', 'Tperpar_e0']
      exosphere_Tperpare0_0vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
      exosphere_Tperpare0_0vtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
      exosphere_Tperpare0_0vtkDisplay.PolarAxes = 'PolarAxesRepresentation'
      exosphere_Tperpare0_0vtkDisplay.ScalarOpacityUnitDistance = 0.29423594422452637
      exosphere_Tperpare0_0vtkDisplay.ScalarOpacityFunction = tperpar_ePWF
      exosphere_Tperpare0_0vtkDisplay.Slice = 7

      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
      exosphere_Tperpare0_0vtkDisplay.ScaleTransferFunction.Points = [-4.969829669326797e-12, 0.0, 0.5, 0.0, 2.335628232685849e-05, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
      exosphere_Tperpare0_0vtkDisplay.OpacityTransferFunction.Points = [-4.969829669326797e-12, 0.0, 0.5, 0.0, 2.335628232685849e-05, 1.0, 0.5, 0.0]

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for tperpar_eLUT in view renderView1
      tperpar_eLUTColorBar = GetScalarBar(tperpar_eLUT, renderView1)
      tperpar_eLUTColorBar.Title = 'Tperpar_e0'
      tperpar_eLUTColorBar.ComponentTitle = 'XX'

      # set color bar visibility
      tperpar_eLUTColorBar.Visibility = 1

      # show color legend
      exosphere_Tperpare0_0vtkDisplay.SetScalarBarVisibility(renderView1, True)

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(exosphere_Tperpare0_0vtk)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [1]}
  coprocessor.SetUpdateFrequencies(freqs)
  if requestSpecificArrays:
    arrays = []
    coprocessor.SetRequestedArrays('input', arrays)
  coprocessor.SetInitialOutputOptions(timeStepToStartOutputAt,forceOutputAtFirstCall)

  if rootDirectory:
      coprocessor.SetRootDirectory(rootDirectory)

  if make_cinema_table:
      coprocessor.EnableCinemaDTable()

  return coprocessor


#--------------------------------------------------------------
# Global variable that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView and the update frequency
coprocessor.EnableLiveVisualization(False, 1)

# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=rescale_lookuptable,
        image_quality=0, padding_amount=imageFileNamePadding)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)

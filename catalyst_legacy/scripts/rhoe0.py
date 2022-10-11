
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
      renderView1.CameraPosition = [3.468647663212199, 25.572389073024155, 6.194146857011161]
      renderView1.CameraFocalPoint = [3.9999959468841553, 3.9999959468841553, 0.9999974966049194]
      renderView1.CameraViewUp = [0.0406181700205616, 0.23484007148403344, -0.9711850004451031]
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
      exosphere_rhoe0_0vtk = coprocessor.CreateProducer(datadescription, 'input')

      # create a new 'Slice'
      slice1 = Slice(Input=exosphere_rhoe0_0vtk)
      slice1.SliceType = 'Plane'
      slice1.SliceOffsetValues = [0.0]

      # init the 'Plane' selected for 'SliceType'
      slice1.SliceType.Origin = [3.9999960000000003, 3.9999960000000003, 0.9999975000000001]
      slice1.SliceType.Normal = [0.0, 1.0, 0.0]

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from slice1
      slice1Display = Show(slice1, renderView1)

      # get color transfer function/color map for 'rhoe0'
      rhoeLUT = GetColorTransferFunction('rhoe0')
      rhoeLUT.RGBPoints = [-1.1581283807754517, 0.231373, 0.298039, 0.752941, -0.5790641903877258, 0.865003, 0.865003, 0.865003, 0.0, 0.705882, 0.0156863, 0.14902]
      rhoeLUT.ScalarRangeInitialized = 1.0

      # trace defaults for the display properties.
      slice1Display.Representation = 'Surface'
      slice1Display.ColorArrayName = ['POINTS', 'rhoe0']
      slice1Display.LookupTable = rhoeLUT
      slice1Display.OSPRayScaleArray = 'rhoe0'
      slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      slice1Display.SelectOrientationVectors = 'None'
      slice1Display.ScaleFactor = 0.7999991893768311
      slice1Display.SelectScaleArray = 'rhoe0'
      slice1Display.GlyphType = 'Arrow'
      slice1Display.GlyphTableIndexArray = 'rhoe0'
      slice1Display.GaussianRadius = 0.03999995946884155
      slice1Display.SetScaleArray = ['POINTS', 'rhoe0']
      slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
      slice1Display.OpacityArray = ['POINTS', 'rhoe0']
      slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
      slice1Display.DataAxesGrid = 'GridAxesRepresentation'
      slice1Display.PolarAxes = 'PolarAxesRepresentation'

      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
      slice1Display.ScaleTransferFunction.Points = [-1.0671108961105347, 0.0, 0.5, 0.0, 0.0, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
      slice1Display.OpacityTransferFunction.Points = [-1.0671108961105347, 0.0, 0.5, 0.0, 0.0, 1.0, 0.5, 0.0]

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for rhoeLUT in view renderView1
      rhoeLUTColorBar = GetScalarBar(rhoeLUT, renderView1)
      rhoeLUTColorBar.Title = 'rhoe0'
      rhoeLUTColorBar.ComponentTitle = ''

      # set color bar visibility
      rhoeLUTColorBar.Visibility = 1

      # show color legend
      slice1Display.SetScalarBarVisibility(renderView1, True)

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get opacity transfer function/opacity map for 'rhoe0'
      rhoePWF = GetOpacityTransferFunction('rhoe0')
      rhoePWF.Points = [-1.1581283807754517, 0.0, 0.5, 0.0, 0.0, 1.0, 0.5, 0.0]
      rhoePWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(slice1)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [1, 1]}
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


#--------------------------------------------------------------

# Global timestep output options
timeStepToStartOutputAt=0
forceOutputAtFirstCall=False

# Global screenshot output options
imageFileNamePadding=4
rescale_lookuptable=False

# Whether or not to request specific arrays from the adaptor.
requestSpecificArrays=False

# a root directory under which all Catalyst output goes
rootDirectory='./data/images'

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
      renderView1.ViewSize = [1280, 720]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [3.9999959468841553, 3.9999959468841553, 0.9999974966049194]
      renderView1.StereoType = 'Crystal Eyes'
      renderView1.CameraPosition = [3.9999959468841553, 3.9999959468841553, 23.195259688466994]
      renderView1.CameraFocalPoint = [3.9999959468841553, 3.9999959468841553, 0.9999974966049194]
      renderView1.CameraFocalDisk = 1.0
      renderView1.CameraParallelScale = 5.744556566297824
      renderView1.Background = [0.32, 0.34, 0.43]
      renderView1.BackEnd = 'OSPRay raycaster'
      renderView1.OSPRayMaterialLibrary = materialLibrary1

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='image_%t.png', freq=1, fittoscreen=0, magnification=1, width=1280, height=720, cinema={}, compression=-1)
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
      exosphere_Je_ = coprocessor.CreateProducer(datadescription, 'input')

      # create a new 'Stream Tracer'
      streamTracer2 = StreamTracer(Input=exosphere_Je_,
          SeedType='Point Source')
      streamTracer2.Vectors = ['POINTS', 'Ve0']
      streamTracer2.MaximumStreamlineLength = 7.999992000000001

      # init the 'Point Source' selected for 'SeedType'
      streamTracer2.SeedType.Center = [1.0, 3.9999960000000003, 0.9999975000000001]
      streamTracer2.SeedType.NumberOfPoints = 60
      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from streamTracer2
      streamTracer2Display = Show(streamTracer2, renderView1)

      # trace defaults for the display properties.
      streamTracer2Display.Representation = 'Surface'
      streamTracer2Display.ColorArrayName = [None, '']
      streamTracer2Display.OSPRayScaleArray = 'AngularVelocity'
      streamTracer2Display.OSPRayScaleFunction = 'PiecewiseFunction'
      streamTracer2Display.SelectOrientationVectors = 'Normals'
      streamTracer2Display.ScaleFactor = 0.7883861541748047
      streamTracer2Display.SelectScaleArray = 'AngularVelocity'
      streamTracer2Display.GlyphType = 'Arrow'
      streamTracer2Display.GlyphTableIndexArray = 'AngularVelocity'
      streamTracer2Display.GaussianRadius = 0.039419307708740234
      streamTracer2Display.SetScaleArray = ['POINTS', 'AngularVelocity']
      streamTracer2Display.ScaleTransferFunction = 'PiecewiseFunction'
      streamTracer2Display.OpacityArray = ['POINTS', 'AngularVelocity']
      streamTracer2Display.OpacityTransferFunction = 'PiecewiseFunction'
      streamTracer2Display.DataAxesGrid = 'GridAxesRepresentation'
      streamTracer2Display.PolarAxes = 'PolarAxesRepresentation'

      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
      streamTracer2Display.ScaleTransferFunction.Points = [-0.010548283550484048, 0.0, 0.5, 0.0, 0.010489982733140888, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
      streamTracer2Display.OpacityTransferFunction.Points = [-0.010548283550484048, 0.0, 0.5, 0.0, 0.010489982733140888, 1.0, 0.5, 0.0]

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(streamTracer2)
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

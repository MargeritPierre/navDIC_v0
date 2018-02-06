# navDIC_v0

An Open Source Matlab DIC App. *Beta version*.

Written in [Laboratoire NAVIER](https://navier.enpc.fr/?lang=en), Ecole des Ponts ParisTech, Champs-sur-Marne, France.

The main motivations behind the writing of this code are:
- acquire images **simultaneously** with other data sources *(Force, Temperatureor any sensor)*
- define the complete DIC setup in one app: cameras, inputs, seeds, previews.
- allow for real-time DIC processing during experiments
- provide an open framework, with classes of seeds, displacement and strain computation methods that can easily be added/modified

## Requirements

To make the app run, the following Matlab Toolboxes are needed:
- [Image Processing Toolbox](https://www.mathworks.com/products/image.html): to perform the DIC, image format conversion, etc.. (mandatory)
- [Image Acquisition Toolbox](https://www.mathworks.com/products/imaq.html): to control cameras
- [Data Acquisition Toolbox](https://www.mathworks.com/products/daq.html): to control external inputs/outputs

Depending on the used cameras, particular Matlab adaptators packages could be needed too.

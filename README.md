# navDIC_v0

An Open Source Matlab DIC App. *Beta version*.

Written in [Laboratoire NAVIER](https://navier.enpc.fr/?lang=en), Ecole des Ponts ParisTech, Champs-sur-Marne, France.

The main motivations behind the writing of this code are:
- acquire images **simultaneously** with other data sources (Force, Temperature or any sensor).
- define the complete DIC setup in one app: cameras, inputs, seeds, previews.
- allow for real-time DIC processing during experiments.
- provide an open framework, with classes of seeds, displacement and strain computation methods that can easily be added/modified.


## Requirements

To make the app run, the following Matlab Toolboxes are needed:
- [Image Processing Toolbox](https://www.mathworks.com/products/image.html): to perform the DIC, image format conversion, etc.. (mandatory)
- [Image Acquisition Toolbox](https://www.mathworks.com/products/imaq.html): to control cameras
- [Data Acquisition Toolbox](https://www.mathworks.com/products/daq.html): to control external inputs/outputs

Depending on the used cameras, particular [Matlab adaptators packages](https://www.mathworks.com/help/imaq/installing-the-support-packages-for-image-acquisition-toolbox-adaptors.html) could be needed too. The app has been tested and works with the following adaptors:
- [dcam](https://www.mathworks.com/hardware-support/dcam.html) : IEEE 1394 port
- [GigE](https://www.mathworks.com/hardware-support/gige.html) : RJ45 port


## Installation

As any other Matlab Library, just copy/paste the navDIC directory in your matlab root path.

# amPEPpy model

`amPEP.onnx` is a portable conversion of the official pretrained random-forest model
from amPEPpy commit `85aab3428b328d9fe4744052258746d8f4ba7bf6`. It was converted with
scikit-learn 1.4.0 and verified against the source model before deployment. The app
uses the ONNX artifact so inference is independent of the Python/scikit-learn version.

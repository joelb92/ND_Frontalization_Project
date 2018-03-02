This is another approach for frontalizing facial images, which doesn't need any 3D models. It is based on the old idea of view morphing originally introduced in "S. M. Seitz and C. R. Dyer: View Morphing, Proc. SIGGRAPH 96, 1996, 21-30." While frontalization with view morphing is not really that novel, I have not seen any fully automatic implementation of the procedure, as reliable landmark detectors were not available until recently.

Our procedure first mirrors the (input) face image to produce two different view points of the same face. Facial landmarks (68 of them) are then detected on each image and used to produce a frontal view of the face. The frontalization results appear to be visually more pleasing than those produced with the 3D model, especially for  more difficult images.

To run a demo of the procedure replace the path to mexopencv with your own path and run demo_viewmorph.m. There are a lot of setting for the approach and can be tested by changing the configuration in viewmorph.m. One particularly interesting setting is to use SDM for the landmarks on the face and use our CMR approach for the facial outline. For some images this produces better results.    



**Version**: 1.1.0

================  ================  =================  ====================================================================================
Attribute         Type              Array Type         Description                                                                         
================  ================  =================  ====================================================================================
``PYP_SPEC``      str                                  PypeIt spectrograph name                                                            
``bpmtilts``      `numpy.ndarray`_  `numpy.integer`_   Bad pixel mask for tilt solutions. Keys are taken from SlitTraceSetBitmask          
``coeffs``        `numpy.ndarray`_  `numpy.floating`_  2D coefficents for the fit on the initial slits.  One set per slit/order (3D array).
``func2d``        str                                  Function used for the 2D fit                                                        
``nslit``         int                                  Total number of slits.  This can include masked slits                               
``spat_flexure``  float                                Flexure shift from the input TiltImage                                              
``spat_id``       `numpy.ndarray`_  `numpy.integer`_   Slit spat_id                                                                        
``spat_order``    `numpy.ndarray`_  `numpy.integer`_   Order for spatial fit (nslit)                                                       
``spec_order``    `numpy.ndarray`_  `numpy.integer`_   Order for spectral fit (nslit)                                                      
================  ================  =================  ====================================================================================

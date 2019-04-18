# exportCTF
The function *exportCTF* exports EBSD data from *ebsd* objects created with the crystallographic texture analysis tool **MTEX** to *Channel Text Files* (*ctf*) for further processing with interactive applications such as **Oxford Instruments Channel5** or **Atex**. The function was tested with data that was imported into MTEX from Oxford Instruments *cpr/crc* files, but should work with any kind of imported data. 

### Function arguments:
*function exportCTF(ebsd,fName,varargin)*

- *ebsd*:		  MTEX ebsd-object
- *fName*:		Name of output file including file ending ‘*.ctf*’ 		(Example ‘*ebsdDataOut.ctf*’)
- *varargin*:	Variable argument in – optional arguments as parameter-value pairs:
  - *‘params’*: The parameter *params* may be used to define optional microscopy acquisition parameters. By default the microscopy acquisition parameters are set to zeros. *params* may be combined with these values:
    - A structure containing the file information from an imported cpr/crc EBSD file (Example: *…,’params’,cprStruct,…*)
    -	String ‘manual’ which will prompt the user to manually supply the parameters (Example: *…,’params’,’manual’,…*)
  - *‘flip’*:	The parameter *flip* may be paired with Boolean values *0* and *1* and indicates whether the *ebsd* spatial data (not the orientations) should be rotated by 180°. By default this parameter is set to 0 (no rotation).
	(Example: *…,’flip’, 1,…*)

### Contact

If you encounter problems or have suggestions for improvement of this function feel free to contact me via *contactnospam@fniessen.com* (remove the nospam to make this email address work).

*Frank Niessen, 18/04/2019*

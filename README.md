# DTASelect2MzId
Data file format converter, from DTASelect-filter.txt files (output from DTASelect software) to mzIdentML version 1.1.0 and 1.2.0

Usage: DTASelect2MzId -i [input file or folder]  
  
************  
Input file is missing  
************  
 -d,--decoy <arg>                 [OPTIONAL] decoy regular expression.  
                                  Ignores matching entries. Example: 'Reverse'.  
 -i,--input <arg>                 path to the input file (or folder)  
 -n,--file_name <arg>             [OPTIONAL] To use a input file name  
                                  different from the default 'DTASelect-filter.txt'  
 -ns,--no_spectra                 [OPTIONAL] If provided, no MS2 files  
                                  will be readed in order to match PSMs with spectra in the output file  
 -r,--recursive                   [OPTIONAL] In case of using a folder as  
                                  input '-i', it will search recursively for all the DTASelect-filter.txt  
                                  files.  
 -rs,--referenceToSpectra <arg>   [OPTIONAL] Reference to spectra.  
                                  Possible values: 'MZXML', 'MS2'. Default: 'MS2'  
 -sky,--skyline <arg>             [OPTIONAL] Whether the generated  
                                  mzIdentML will be compatible with skyline for importing results. If  
                                  'true', the file will contain additional features to make it compatible.  
                                  Possible values: 'true' or 'false'. Default: 'true'  
 -u,--unique_output_file          [OPTIONAL] A single mzIdentML file will  
                                  be created collapsing all the information of all the input files.  
                                  Otherwise, a mzIdentML file will be created for each input file.  
 -v,--version <arg>               [OPTIONAL] Version of the output  
                                  mzIdentML '-v', Possible values: '1.1.1', '1.2.0'. Default: 1.1.1  
Contact Salvador Martinez-Bartolome at salvador@scripps.edu for more help  

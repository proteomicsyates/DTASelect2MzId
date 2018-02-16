# DTASelect2MzId
Data file format converter, from DTASelect-filter.txt files (output from DTASelect software) to mzIdentML version 1.1.0 and 1.2.0

Download:
Latest version available at: [http://sealion.scripps.edu/dtaselect2mzid/]  

Usage: 
```
DTASelect2MzId -i [input file or folder]  
 
 -d,--decoy <arg>                 [OPTIONAL] decoy regular expression. Ignores matching entries. Example: 'Reverse'.
   
 -i,--input <arg>                 path to the input file (or folder)   
   
 -n,--file_name <arg>             [OPTIONAL] To use a input file name different from the default 'DTASelect-filter.txt'  
    
 -ns,--no_spectra                 [OPTIONAL] If provided, no MS2 files will be readed in order to match PSMs with spectra in the output file  
                                  
 -r,--recursive                   [OPTIONAL] In case of using a folder as input '-i', it will search recursively for all the DTASelect-filter.txt files.
 
 -rs,--referenceToSpectra <arg>   [OPTIONAL] Reference to spectra. Possible values: 'MZXML', 'MS2'. Default: 'MS2'  
 
 -sky,--skyline <arg>             [OPTIONAL] Whether the generated mzIdentML will be compatible with skyline for importing results. If  
                                  'true', the file will contain additional features to make it compatible. Possible values: 'true' or 'false'. Default: 'true'  
                                  
 -u,--unique_output_file          [OPTIONAL] A single mzIdentML file will be created collapsing all the information of all the input files. Otherwise, a mzIdentML file will be created for each input file.  
   
 -v,--version <arg>               [OPTIONAL] Version of the output mzIdentML '-v', Possible values: '1.1.1', '1.2.0'. Default: 1.1.1  
```
  
   
---
  
   
*Note about **referenceToSpectra** parameter:*  
   - The use of the parameter *referenceToSpectra* will make the converter to read the associated spectra files (mzXML or ms2) in order to properly include the reference to the spectra of the peptide spectrum matches.
   - However, the converter will not include the actual identified/matched fragment ions in the output mzIdentML file.
   - As an example, if the converter reads the spectrum line from the DTASelect file as:
```
*	042117_PB1_trypsin.18815.18815.3	1.7635	0.3047	100.0	2275.1868	2275.1765	4.5	1986159.2	2	3.9090989	12.0	31.5	1	G.GGWSGSHAFILVM(15.9949)AALTTRAGR.K 
```  

it will look for the file named as 
```
042117_PB1_trypsin.extension
```
where 'extension' is the parameter value of **-rs** input parameter.
  - If everything goes well, then you will be able to open the output mzIdentML in a software sich as [PRIDE Inspector](https://github.com/PRIDE-Toolsuite/pride-inspector) and then be able to select the associated spectra files, which will load the matches spectra in the screen with their corresponding annotated matched fragment ions.
    
---    
  
 #### How to integrate search engine input parameters in the mzIdentML?
  You should have a ***search.xml*** file in the same forlder as the input files. If it is found, elements such as:
   - *Enzymes*, 
   - *FragmentTolerance*, 
   - *ParentTolerance*, 
   - *SearchType* or
   - *ModificationParams*  
   
 will be added to the  *SpectrumIdentificationProtocol* element of the resulting mzIdentML file.

Contact Salvador Martinez-Bartolome at salvador at scripps.edu for more help  

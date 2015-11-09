OMD-TOIL
===
Copyright (c) 2015, Katherine D. McMahon and Lackeys  
URL: [https://mcmahonlab.wisc.edu/](https://mcmahonlab.wisc.edu/)  
URL: [https://github.com/McMahonLab/](https://github.com/McMahonLab/)  
All rights reserved.

Overview
--
During August 20-21, 2015, the McMahon Lab carried out a 24 hour sampling of Lake Mendota. This repo describes the data and analysis associated with this expedition, known as OMD-TOIL (Operation Mendota Drain: Transcriptomics of Inland Lakes).

Samples were collected approximately every two hours from 6am on August 20 to 6am on August 21. A subset of samples were subject to RNA extraction, spiking with an internal standard, and sent to the [UW Biotech Center](https://www.biotech.wisc.edu/) for sequencing. Resulting sequences were subject to rRNA removal.

Repo Structure
--
    project
    |- README     # the top level description of content
    |- data       # raw and primary data, are not changed once created
    |- metadata   # metadata associated with the samples, not changed once created
    |- protocols  # protocols associated with sample collection and processing
    |- scripts    # protocols associated with processing and analysis of sequence data

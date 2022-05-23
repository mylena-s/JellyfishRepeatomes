# JellyfishRepeatomes
Supplementary files for Chapter 2 of my master's dissertation "Jellyfish repeatomes uncover hidden satellite diversity and provide insights into the C-value enigma"


﻿The files and directories included in this folder correspond to "Supplementary file S1" of the dissertation
 and include the following:

|directory/file			                    	| description
|-------------------------------------------------------|----------------------------------------------------------
|/RepeatExplorer/			                |directory including datasets resulting from RE2 pipeline                    
|                                         		|and the scripts used to process them and carry on                      
|                                        		|statistics.
|                                        		|
|/RepeatExplorer/AssemblyAnalysis.py     		|Python script used to manipulate files resulting from  
|                                        		|repeatmasker repeat annotation in assemblies, 
|                                        		|to identify clusters after clusterscan,
|                                        		|and to classify them into proximal or distal to
|                                        		|gaps and ends.
|							|
|/RepeatExplorer/ScyphoSat_table 		      	|table containing all information regarding satellite
|				                       	|monomers obtained with RE2, including manual curation   
|                               	  	     	|information, classification and sequence information.
|							|
|/RepeatExplorer/assemblySat_plots.py 			|Python script to plot the distribution of intragenomic
|					               	|divergence of satellites in assemblies and reads.
|							|
|/RepeatExplorer/Genus_ComparisonStats.py		|Python script to compare satDNA features distribution in 
|                                         		|samples of the same genus.
|							|
|/RepeatExplorer/SampleStatistics.py	    		|Python script to compare satDNA features across all  
|                                        		|species and samples.
|                                      			|      
|/RepeatExplorer/duplicateMonomer.py     		|Python script to duplicate and rename satDNA monomers.
|							|
|/RepeatExplorer/filterBlastn.py         		|Python script to classify satDNA blast hits as variants, 
|                                        		|families and superfamilies.
|							|
|/dnaPipeTE/		            		        |directory including datasets resulting from dnaPipeTE and 
|                                        		|the scripts used to process them and carry on statistics 
|		  					|
|/dnaPipeTE/SupplementaryTableS2extended 		|general repetitive content information including
|					                |detailed abundance information about all TE Superfamilies
|							|
|/dnaPipeTE/dnaPipeTE_spearmanr.R	      		|R script used for statistical analysis between genome      
|                                        		|size and dnaPipeTE results		
|							|
|/dnaPipeTE/filter_RTcompleteness.py	    		|Python script used for filtering protein domains with  
|                                        		|coverage below certain threshold

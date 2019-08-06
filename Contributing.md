# Contributing

  * All contributions are under the same licence (MIT) as the rest of the repository. 
  * Assembly code should go to the `src/` directory followed by the primary assembly at the second level and assembly name at the third directory level.
  * General code should go to the `src/scripts` directory.
  * Setup scripts should
     * be named `prepare.sh` and located in the assembly source code directory
     * install assemblies in a `$primaryAssembly/$assemblyName` directory following the general structure for assemblies (see Readme.md)
     * use MD5 checksums to check all downloaded files
     * Bash is fine, also e.g. Make, Snakemake, etc.
     * All tools should be in Conda. The Conda environment should be updated to include new tools/libraries used by new code.
     * Write reusable code (document inputs and output; explain why you made certain decisions; etc.)
  

     
  
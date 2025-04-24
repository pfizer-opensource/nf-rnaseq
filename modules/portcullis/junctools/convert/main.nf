/*
Nextflow module to run Portcullis junctools convert

Requires: none

Params: none

Input:
  -- input_file (path): path to the input file to convert
  -- input_type (string): format of the input file (see: https://portcullis.readthedocs.io/en/latest/junctools.html#convert)
  -- output_type (string): format of the output file (see: https://portcullis.readthedocs.io/en/latest/junctools.html#convert)

Output:
  -- path to the converted output file, emitted as: `converted`
  -- _version.yaml containing the version string, emitted as: `version`   
     
Customization:
    ext.args in module.config to add optional arguments
    publishDir in module.config
*/

nextflow.enable.dsl = 2

process JUNCTOOLS_CONVERT {
    tag "$input_file"

    label 'process_6cpu_12gb_6h'

    input:
    path input_file
    val input_type
    val output_type
  
    output:
    path '[!_]*', emit: converted
    path "_version.yaml", emit: version

    script:
    def outname = input_file.getName() + ".${output_type}"
    def args = task.ext.args ? task.ext.args.join(' \\\n      ').trim() : ''
    """
    junctools convert \\
      -if $input_type \\
      -of $output_type \\
      -o $outname \\
      $args \\
      $input_file

    cat <<-END_VERSION >_version.yaml
    "${task.process}":
        junctools: \$(junctools --version 2>&1 | sed -nre 's/^[^0-9]*(([0-9]+\\.)*[0-9]+).*/\\1/p')
    END_VERSION
    """
}

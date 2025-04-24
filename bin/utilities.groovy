// Define helper function to create metadata, file path tuples
def create_fastq_entry(LinkedHashMap row) {
    def meta = [:]
    meta.id = row.sample_name
    meta.stranded = row.stranded
    meta.orientation = row.orientation

    def array = []
    fastq_files = file("${params.fastq_dir}/${row.sample_name}*.{fastq,fq}{,.gz}")
    for(fastq in fastq_files) {
        def is_gzip = "file -L ${fastq}".toString().execute() | 'grep -cim1 gzip'.execute()
        is_gzip.waitFor()
        if(!(is_gzip.text.toBoolean())) {
            exit 1, "Sample ${fastq.getName()} is not compressed. To conserve disk space compress all input fastq files with gzip."
        }
    }
    if( fastq_files.isEmpty() ) {
        if(params.allow_missing_samples) {
            log.warn("Sample ${row.sample_name} has no matching files")
            return null
            } else {
                exit 1, "Sample ${row.sample_name} has no matching files"
            }
    } else {
        switch( fastq_files.size ) {
            case 1: 
                meta.is_paired_end = false;
                break;
            case 2: 
                meta.is_paired_end = true;
                break;
            default: 
                exit 1, "Sample ${row.sample_name} has more than 2 matching samples.";
                break;
        }
        // verify that the correct number of files were detected based on SampleSheet designation of single- vs. paired-end
        if (meta.is_paired_end && row.read_type =='single') {
            exit 1, "Sample ${row.sample_name} is read_type 'single' in SampleSheet but detected 2 fastq files with matching sample name."
        }
        if (!meta.is_paired_end && row.read_type =='paired') {
            exit 1, "Sample ${row.sample_name} is read_type 'paired' in SampleSheet but only detected 1 fastq file with matching sample name."
        }
        array = [ meta, fastq_files.sort() ]
    }
    return array
}

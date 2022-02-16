package com.hartwig.hmftools.teal.tellength

import com.beust.jcommander.Parameter
import com.beust.jcommander.ParameterException

enum class SampleType
{
    tumor, ref
}

// NOTE: We use nullable string for required argument types, this is to avoid (default: <empty string>) getting
// printed in the usage by JCommander
data class TelLengthParams
(
    @Parameter(names = ["-sample"], required = true, description = "ID of sample")
    var sampleId: String? = null,

    @Parameter(names = ["-type"], required = true, description = "Type of sample")
    var sampleType: SampleType? = null,

    @Parameter(names = ["-output_file"], required = true, description = "Path to output tsv file")
    var outputFile: String? = null,

    @Parameter(names = ["-telbam_file"], required = true, description = "Path to input telbam file")
    var telbamFile: String? = null,

    @Parameter(names = ["-germline_telomere_length"],
        description = "Germline telomere length. Only required for determining tumor length in a mixed sample")
    var germlineTelomereLength: Double? = null,

    @Parameter(names = ["-purity"],
        description = "Purity of the sample. Only required for determining tumor length in a mixed sample")
    var purity: Double = 1.0,

    @Parameter(names = ["-ploidy"],
        description = "Ploidy of the sample. Only required for determining tumor length in a mixed sample")
    var ploidy: Double = 2.0,

    @Parameter(names = ["-duplicate_proportion"],
        description = "Proportion of reads that are marked duplicates in the sample BAM")
    var duplicatePercent: Double = 0.0,

    @Parameter(names = ["-mean_reads_per_kb"], required = true, description = "Mean reads per KB")
    var meanReadsPerKb: Int = 0,

    @Parameter(names = ["-gc50_reads_per_kb"], description = "GC 50 reads per KB. Defaults to mean reads per KB if not provided")
    var gc50ReadsPerKb: Int? = null
)
{
    fun validate()
    {
        if (duplicatePercent < 0.0 || duplicatePercent >= 1.0)
        {
            throw ParameterException("invalid arg: $duplicatePercent for duplicatePercent, must be between [0, 1)")
        }
    }
}
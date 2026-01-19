package com.hartwig.hmftools.teal.tellength

import com.beust.jcommander.Parameter
import com.beust.jcommander.ParameterException
import com.hartwig.hmftools.common.sequencing.SequencingType

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

    @Parameter(names = ["-" + SequencingType.SEQUENCING_TYPE_CFG], description = "Sequencing types")
    var sequencingType: SequencingType = SequencingType.ILLUMINA,

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

    @Parameter(names = ["-mean_read_depth"], required = true, description = "Mean read depth")
    var meanReadDepth: Double = 0.0,

    @Parameter(names = ["-gc50_read_depth"], description = "GC 50 read depth. Defaults to mean reads per KB if not provided")
    var gc50ReadDepth: Double? = null
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
package com.hartwig.hmftools.teal.breakend

import com.beust.jcommander.IStringConverter
import com.beust.jcommander.Parameter
import com.beust.jcommander.ParameterException
import com.hartwig.hmftools.common.genome.bed.NamedBedFile
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import com.hartwig.hmftools.common.genome.region.GenomeRegion
import com.hartwig.hmftools.common.genome.region.GenomeRegions

private const val DEFAULT_BREAK_POINT_MARK_DUP_DISTANCE = 60
private const val DEFAULT_SPLIT_TELOMERE_MATCH_THRESHOLD = 0.9

class BreakEndParams
{
    // the distance which we deem the same breakpoint and mark as duplicates
    var markDuplicateDistance: Int = DEFAULT_BREAK_POINT_MARK_DUP_DISTANCE
    var telomereMatchThreshold: Double = DEFAULT_SPLIT_TELOMERE_MATCH_THRESHOLD
    // alignedSegmentTelomereMatchThreshold is used to determine if the aligned part
    // of the split read is also telomeric
    // I find using 0.8 can still get all the break ends gripss found.
    // using 0.7 and we start excluding some
    var alignedSegmentTelomereMatchThreshold: Double = 0.8

    // filter for poly base
    var maxAlignedPolyBaseThreshold: Double = 0.85

    @Parameter(names = ["-sample_id"], required = true, description = "ID of tumor sample")
    lateinit var sampleId: String

    @Parameter(names = ["-reference_telbam"], required = true, description = "Path to reference telbam.bam file")
    lateinit var germlineTelbamFile: String

    @Parameter(names = ["-tumor_telbam"], required = true, description = "Path to tumor telbam.bam file")
    lateinit var tumorTelbamFile: String

    @Parameter(names = ["-output"], required = true, description = "Output breakend Gzipped TSV (tsv.gz) file")
    lateinit var outputFile: String

    @Parameter(names = ["-" + RefGenomeVersion.REF_GENOME_VERSION],
                description = RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC)
    var refGenomeVersionStr: String = "37"
        set(value) {
            field = value
            refGenomeVersion = RefGenomeVersion.from(value)
            excludedGenomeRegions = loadExcludedRegionBed(refGenomeVersion)
        }

    var refGenomeVersion = RefGenomeVersion.V37
    var excludedGenomeRegions: List<GenomeRegion> = emptyList()

    @Parameter(names = ["-genome_regions"],
        listConverter = IncludedGenomeRegionsConverter::class,
        description = "only these genome regions will be processed, in the form of 1:10000-2000;X:1400-2000")
    var includedGenomeRegions: List<GenomeRegion>? = null

    companion object
    {
        private fun loadExcludedRegionBed(refGenomeVersion: RefGenomeVersion): List<GenomeRegion>
        {
            val resourcePath = when (refGenomeVersion)
            {
                RefGenomeVersion.V37 -> "blacklistedTelomereRegions.37.bed"
                RefGenomeVersion.V38 -> "blacklistedTelomereRegions.38.bed"
                else -> null
            }
            if (resourcePath != null)
            {
                val bedStream: java.io.InputStream = BreakEndParams::class.java.classLoader.getResourceAsStream(resourcePath)!!
                // write the resource out to a temp file and read it back
                val tempFile = java.io.File.createTempFile("teal-region-bed", null)
                tempFile.deleteOnExit()

                // cannot use readAllBytes since we are still on java8
                val bytes = ByteArray(bedStream.available())
                java.io.DataInputStream(bedStream).readFully(bytes)
                tempFile.writeBytes(bytes)
                return NamedBedFile.readBedFile(tempFile.path)
            }
            return emptyList()
        }

        class IncludedGenomeRegionsConverter : IStringConverter<List<GenomeRegion>?>
        {
            override fun convert(value: String): List<GenomeRegion>?
            {
                return parseIncludedGenomeRegions(value)
            }
        }

        //
        // Parse the argument string for genome regions in the form of "1:10000-2000;X:1400-2000"
        private fun parseIncludedGenomeRegions(incGenomeRegionArg: String?) : List<GenomeRegion>?
        {
            if (incGenomeRegionArg == null)
            {
                return null
            }

            val genomeRegions = ArrayList<GenomeRegion>()

            val pattern = Regex("""(.+):(\d+)-(\d+)""")
            for (s in incGenomeRegionArg.split(";"))
            {
                val matchResult = pattern.matchEntire(s) ?:
                throw ParameterException("Invalid value $incGenomeRegionArg for option '-genome_regions'")

                val (chromosome, start, end) = matchResult.destructured
                genomeRegions.add(GenomeRegions.create(chromosome, start.toInt(), end.toInt()))
            }

            return genomeRegions
        }
    }
}
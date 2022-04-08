package com.hartwig.hmftools.teal.tellength

import com.hartwig.hmftools.teal.ReadGroup
import com.hartwig.hmftools.teal.TealUtils
import org.apache.logging.log4j.LogManager
import tech.tablesaw.api.DoubleColumn
import tech.tablesaw.api.IntColumn
import tech.tablesaw.api.StringColumn
import tech.tablesaw.api.Table
import tech.tablesaw.columns.numbers.NumberColumnFormatter
import tech.tablesaw.io.csv.CsvWriteOptions
import java.text.NumberFormat
import java.util.*

// analyse the telo bam and come up with our own score
// of F1, F2, F4
class TelomereLengthCalc(
    val purity: Double,
    val ploidy: Double,
    val duplicateProportion: Double,
    val meanReadsPerKb: Double,
    val gc50ReadsPerKb: Double,
    val germlineTelomereLength: Double?)
{
    private val logger = LogManager.getLogger(javaClass)

    enum class FragmentType
    {
        UNKNOWN,
        NOT_TELOMERE,
        F1, // both sides telomeric
        F2, // one side telomeric C rich
        F4  // one side telomeric G rich
    }

    private val mFragmentTypeCount: MutableMap<FragmentType, Int> = EnumMap(FragmentType::class.java)
    fun getFragmentTypeCount(fragType: FragmentType): Int
    {
        return mFragmentTypeCount.getOrDefault(fragType, 0)
    }

    fun processReadGroups(readGroups: Collection<ReadGroup>)
    {
        for (readGroup in readGroups)
        {
            onReadGroup(readGroup)
        }

        logger.info(
            "frag type counts: F1({}) F2({}) F4({})",
            getFragmentTypeCount(FragmentType.F1),
            getFragmentTypeCount(FragmentType.F2),
            getFragmentTypeCount(FragmentType.F4))

        logger.info("gc bias adj: {}", calcGcBiasAdj())
        logger.info("sample mix length: {}", calcTelomereLength())
        logger.info("tumor only length: {}", calcTumorTelomereLength())
    }

    fun onReadGroup(readGroup: ReadGroup)
    {
        // we want to classify this read group
        val fragType = classifyFragment(readGroup)
        mFragmentTypeCount[fragType] = mFragmentTypeCount.getOrDefault(fragType, 0) + 1
    }

    fun classifyFragment(readGroup: ReadGroup): FragmentType
    {
        if (readGroup.Reads.size < 2)
        {
            return FragmentType.UNKNOWN
        }

        // we want to check
        val read1 = readGroup.firstOfPair
        val read2 = readGroup.secondOfPair

        if (read1 == null || read2 == null)
        {
            return FragmentType.UNKNOWN
        }

        var seq1 = read1.readString!!
        var seq2 = read2.readString!!

        // qual_str1 = read1['BaseQualities']
        // qual_str2 = read2['BaseQualities']
        if (read1.readNegativeStrandFlag)
        {
            seq1 = TealUtils.reverseComplementSequence(seq1)
        }
        if (read2.readNegativeStrandFlag)
        {
            seq2 = TealUtils.reverseComplementSequence(seq2)
        }

        // try each one in tern
        var fragType = classifyFragment(seq1, seq2)
        if (fragType == FragmentType.NOT_TELOMERE)
        {
            fragType = classifyFragment(seq2, seq1)
        }
        return fragType
    }

    private fun classifyFragment(readPairG: String, readPairC: String): FragmentType
    {
        val isGTelomeric = TealUtils.isLikelyGTelomeric(readPairG)
        val isCTelomeric = TealUtils.isLikelyCTelomeric(readPairC)

        // first we want to find the TTAGGG motif, then fill it up backwards
        if (isGTelomeric && isCTelomeric)
        {
            return FragmentType.F1
        }
        if (isGTelomeric)
        {
            // only g is telomeric
            return FragmentType.F4
        }
        return if (isCTelomeric)
        {
            // only C term is telomeric, could be F2a or F2b
            FragmentType.F2
        }
        else
        {
            // not telomeric
            FragmentType.NOT_TELOMERE
        }
    }

    fun calcGcBiasAdj(): Double
    {
        val gc50bias = gc50ReadsPerKb / meanReadsPerKb

        if (gc50bias < 0.625) return 1.2427
        if (gc50bias < 0.675) return 1.0595
        if (gc50bias < 0.725) return 1.0263
        if (gc50bias < 0.775) return 0.9965
        if (gc50bias < 0.825) return 0.9789
        if (gc50bias < 0.875) return 0.9721
        if (gc50bias < 0.925) return 0.9941780963
        if (gc50bias < 0.975) return 0.9926282328
        if (gc50bias < 1.025) return 0.9754729509
        if (gc50bias < 1.075) return 0.9953500389
        if (gc50bias < 1.125) return 1.0497
        else                 return 1.0497
    }

    // this is the number of reads that have telomeric content - number of interstitial reads
    fun calcTotalTelomericReads() : Int
    {
        return 2 * getFragmentTypeCount(FragmentType.F1) + getFragmentTypeCount(FragmentType.F2) - getFragmentTypeCount(FragmentType.F4)
    }

    fun calcTelomereLength() : Double
    {
        val totalTelomericReads = calcTotalTelomericReads()
        val meanTelomereLength = totalTelomericReads * (1 - duplicateProportion) * 1000 / meanReadsPerKb / calcGcBiasAdj() / 46
        return meanTelomereLength
    }

    fun calcTumorTelomereLength() : Double
    {
        val sampleLength = calcTelomereLength()

        if (germlineTelomereLength == null)
        {
            // if we don't have germline sample then we just return sample length
            return sampleLength
        }

        val tumorLength = (sampleLength * (purity * ploidy + 2 * (1 - purity)) - germlineTelomereLength * (1 - purity) * 2) / purity / ploidy
        return tumorLength
    }

    fun writeTelLengthTsv(lengthTsvPath: String, sampleId: String, sampleType: SampleType)
    {
        val teloLengthTable = Table.create("TealLength")
        teloLengthTable.addColumns(
            StringColumn.create("sampleId"),
            StringColumn.create("type"),
            createDoubleColumn("rawTelomereLength", 2),
            createDoubleColumn("finalTelomereLength", 2),
            IntColumn.create("fullFragments"),
            IntColumn.create("cRichPartialFragments"),
            IntColumn.create("gRichPartialFragments"),
            IntColumn.create("totalTelomericReads"),
            DoubleColumn.create("purity"),
            DoubleColumn.create("ploidy"),
            DoubleColumn.create("duplicateProportion"),
            IntColumn.create("meanReadsPerKb"),
            IntColumn.create("gc50ReadsPerKb"))

        val row = teloLengthTable.appendRow()
        row.setString("sampleId", sampleId)
        row.setString("type", sampleType.name)
        row.setDouble("rawTelomereLength", calcTelomereLength())
        row.setDouble("finalTelomereLength", if (sampleType == SampleType.ref) calcTelomereLength() else calcTumorTelomereLength())
        row.setInt("fullFragments", getFragmentTypeCount(FragmentType.F1))
        row.setInt("cRichPartialFragments", getFragmentTypeCount(FragmentType.F2))
        row.setInt("gRichPartialFragments", getFragmentTypeCount(FragmentType.F4))
        row.setInt("totalTelomericReads", calcTotalTelomericReads())
        row.setDouble("purity", purity)
        row.setDouble("ploidy", ploidy)
        row.setDouble("duplicateProportion", duplicateProportion)
        row.setDouble("meanReadsPerKb", meanReadsPerKb)
        row.setDouble("gc50ReadsPerKb", gc50ReadsPerKb)

        try
        {
            val writeOptions = CsvWriteOptions.builder(lengthTsvPath).separator('\t').usePrintFormatters(true).build()
            teloLengthTable.write().csv(writeOptions)
        }
        catch (e: java.io.IOException)
        {
            throw IllegalStateException("Could not save to csv file: $lengthTsvPath: ${e.message}")
        }
    }

    companion object
    {
        fun createDoubleColumn(name: String, decimals: Int) : DoubleColumn
        {
            val col = DoubleColumn.create(name)

            val format = NumberFormat.getInstance(Locale.ENGLISH)
            format.isGroupingUsed = false
            format.maximumFractionDigits = decimals
            col.setPrintFormatter(NumberColumnFormatter(format, "0"))
            return col
        }
    }
}
package com.hartwig.hmftools.teal.tellength

import com.hartwig.hmftools.common.sequencing.SequencingType
import com.hartwig.hmftools.common.teal.ImmutableTelomereLength
import com.hartwig.hmftools.common.teal.TelomereLengthFile
import com.hartwig.hmftools.teal.ReadGroup
import org.apache.logging.log4j.Level
import org.apache.logging.log4j.LogManager
import java.util.*

// analyse the telo bam and come up with our own score
// of F1, F2, F4
class TelomereLengthCalc(
    val sequencingType: SequencingType,
    val purity: Double,
    val ploidy: Double,
    val duplicateProportion: Double,
    val meanReadDepth: Double,
    val gc50ReadDepth: Double,
    val germlineTelomereLength: Double?)
{
    private val logger = LogManager.getLogger(javaClass)

    // following fields used to calculate mean read length
    private var readLengthSum: Long = 0
    private var numReads: Long = 0
    private val fragmentClassifier = FragmentClassifier(sequencingType)

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
        logger.printf(Level.INFO, "sample mix length: %.2f", calcTelomereLength())
        logger.printf(Level.INFO,"tumor only length: %.2f", calcTumorTelomereLength())
    }

    fun onReadGroup(readGroup: ReadGroup)
    {
        // we want to classify this read group
        val fragType = fragmentClassifier.classifyFragment(readGroup)
        mFragmentTypeCount[fragType] = mFragmentTypeCount.getOrDefault(fragType, 0) + 1

        readGroup.allReads.forEach { r -> readLengthSum += r.readLength }
        numReads += readGroup.allReads.size
    }

    fun calcGcBiasAdj(): Double
    {
        val gc50bias = gc50ReadDepth / meanReadDepth

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
        return (2 * getFragmentTypeCount(FragmentType.F1) + getFragmentTypeCount(FragmentType.F2) - getFragmentTypeCount(FragmentType.F4))
            .coerceAtLeast(0)
    }

    fun calcTelomereLength() : Double
    {
        // protect against division by 0
        if (numReads == 0L)
            return 0.0

        val meanReadLength = readLengthSum / numReads.toDouble()
        val totalTelomericReads = calcTotalTelomericReads()
        val meanTelomereLength = totalTelomericReads * (1 - duplicateProportion) * meanReadLength / meanReadDepth / calcGcBiasAdj() / 46
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
        val telomereLength = ImmutableTelomereLength.builder()
            .type(sampleType.name)
            .rawTelomereLength(calcTelomereLength())
            .finalTelomereLength(if (sampleType == SampleType.ref) calcTelomereLength() else calcTumorTelomereLength())
            .fullFragments(getFragmentTypeCount(FragmentType.F1))
            .cRichPartialFragments(getFragmentTypeCount(FragmentType.F2))
            .gRichPartialFragments(getFragmentTypeCount(FragmentType.F4))
            .totalTelomericReads(calcTotalTelomericReads())
            .purity(purity)
            .ploidy(ploidy)
            .duplicateProportion(duplicateProportion)
            .meanReadDepth(meanReadDepth)
            .gc50ReadDepth(gc50ReadDepth)
            .build()

        TelomereLengthFile.write(lengthTsvPath, sampleId, telomereLength)
    }
}
package com.hartwig.hmftools.gripss

import com.hartwig.hmftools.bedpe.Breakend
import com.hartwig.hmftools.bedpe.Breakpoint
import com.hartwig.hmftools.extensions.*
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder

const val SHORT_EVENT_SIZE = 1000

class StructuralVariantContext(private val context: VariantContext, normalOrdinal: Int = 0, tumorOrdinal: Int = 1) {
    private companion object {
        private val polyG = "G".repeat(16);
        private val polyC = "C".repeat(16);
    }

    val contig = context.contig!!

    val imprecise = context.imprecise()
    val precise = !imprecise
    val variantType = context.toVariantType()
    val orientation = variantType.startOrientation
    val isSingle = variantType is Single
    val isShortDup = variantType is Duplication && variantType.length < SHORT_EVENT_SIZE
    val isShortDel = variantType is Deletion && variantType.length < SHORT_EVENT_SIZE
    val isShortIns = variantType is Insertion && variantType.length < SHORT_EVENT_SIZE
    val isShort = isShortDup || isShortDel || isShortIns
    val vcfId: String = context.id
    val mateId: String? = context.mate()
    val confidenceInterval = context.confidenceInterval()
    val start = context.start
    private val remoteConfidenceInterval = context.confidenceInterval()

    val startBreakend: Breakend = Breakend(contig, start + confidenceInterval.first, start + confidenceInterval.second, orientation)
    val endBreakend: Breakend? = (variantType as? Paired)?.let { Breakend(it.otherChromosome, it.otherPosition + remoteConfidenceInterval.first, it.otherPosition + remoteConfidenceInterval.second, it.endOrientation) }
    val breakpoint: Breakpoint? = endBreakend?.let { Breakpoint(startBreakend, it) }
    val minStart = startBreakend.start

    val maxStart = startBreakend.end
    val insertSequenceLength = variantType.insertSequence.length
    val duplicationLength = (variantType as? Duplication)?.let { it.length + 1 } ?: 0
    val qual = context.phredScaledQual


    private val normalGenotype = context.getGenotype(normalOrdinal);
    private val tumorGenotype = context.getGenotype(tumorOrdinal);
    private val tumorAF = tumorGenotype.allelicFrequency(isSingle, isShort)


    fun confidenceIntervalsOverlap(other: StructuralVariantContext): Boolean {
        return contig == other.contig && other.minStart <= maxStart && other.maxStart >= minStart
    }

    fun context(localLink: String, remoteLink: String, altPath: String?, isHotspot: Boolean, filters: Set<String>): VariantContext {
        val builder = VariantContextBuilder(context).filters()

        builder.attribute(TAF, tumorAF)
                .attribute(LOCAL_LINKED_BY, localLink)
                .attribute(REMOTE_LINKED_BY, remoteLink)
                .attribute(HOTSPOT, isHotspot)

        altPath?.let { x -> builder.attribute(ALT_PATH, x) }
        filters.forEach { x -> builder.filter(x) }
        if (filters.isEmpty()) {
            builder.filter(PASS)
        }

        return builder.make()
    }

    fun assemblies(): List<String> = context.assemblies();

    fun isHardFilter(config: GripssFilterConfig) = normalSupportFilter(config.maxNormalSupport)

    fun softFilters(config: GripssFilterConfig): Set<String> {
        val result = mutableSetOf<String>()

        if (normalCoverageFilter(config.minNormalCoverage)) {
            result.add(MIN_NORMAL_COVERAGE)
        }

        if (allelicFrequencyFilter(config.minTumorAF)) {
            result.add(MIN_TUMOR_AF)
        }

        if (strandBiasFilter(config.maxShortStrandBias)) {
            result.add(SHORT_STRAND_BIAS)
        }

        if (qualFilter(config.minQualBreakEnd, config.minQualBreakPoint)) {
            result.add(MIN_QUAL)
        }

        if (impreciseFilter()) {
            result.add(IMPRECISE)
        }

        if (polyGCFilter()) {
            result.add(MAX_POLY_G_LENGTH)
        }

        if (homologyLengthFilter(config.maxHomLength)) {
            result.add(MAX_HOM_LENGTH)
        }

        if (homologyLengthFilterShortInversion(config.maxHomLengthShortInversion)) {
            result.add(MAX_HOM_LENGTH_SHORT_INV)
        }

        if (inexactHomologyLengthFilter(config.maxInexactHomLength)) {
            result.add(MAX_INEXACT_HOM_LENGTH)
        }

        if (inexactHomologyLengthShortDelFilter(config.maxInexactHomLengthShortDel)) {
            result.add(MAX_INEXACT_HOM_LENGTH_SHORT_DEL)
        }

        if (shortSplitReadTumorFilter()) {
            result.add(SHORT_SR_SUPPORT)
        }

        if (shortSplitReadNormalFilter()) {
            result.add(SHORT_SR_NORMAL)
        }

        if (longDPSupportFilter()) {
            result.add(LONG_DP_SUPPORT)
        }

        if (breakendAssemblyReadPairsFilter()) {
            result.add(BREAK_END_ASSEMBLY_READ_PAIR)
        }

        if (minSizeFilter(config.minSize)) {
            result.add(MIN_SIZE)
        }

        return result
    }

    fun qualFilter(minQualBreakEnd: Int, minQualBreakPoint: Int): Boolean {
        val minQual = if (isSingle) minQualBreakEnd else minQualBreakPoint
        return context.phredScaledQual < minQual.toDouble()
    }

    fun polyGCFilter(): Boolean {
        return isSingle && variantType.insertSequence.contains(polyG) || variantType.insertSequence.contains(polyC)
    }

    fun inexactHomologyLengthFilter(maxInexactHomLength: Int): Boolean {
        return !isSingle && !isShortDup && context.inexactHomologyLength() > maxInexactHomLength
    }

    fun inexactHomologyLengthShortDelFilter(maxInexactHomLength: Int, minDelLength: Int = 100, maxDelLength: Int = 800): Boolean {
        return variantType is Deletion && variantType.length >= minDelLength && variantType.length <= maxDelLength && context.inexactHomologyLength() > maxInexactHomLength;
    }

    fun breakendAssemblyReadPairsFilter(): Boolean {
        return isSingle && context.breakendAssemblyReadPairs() == 0
    }

    fun minSizeFilter(minSize: Int): Boolean {
        return when (variantType) {
            is Deletion -> variantType.length + variantType.insertSequence.length - 1 < minSize
            is Insertion -> variantType.length + variantType.insertSequence.length + 1 < minSize
            is Duplication -> variantType.length + variantType.insertSequence.length < minSize
            else -> false
        }
    }


    fun homologyLengthFilter(maxHomLength: Int): Boolean {
        return !isSingle && context.homologyLength() > maxHomLength
    }

    fun homologyLengthFilterShortInversion(maxHomLength: Int, maxInversionLength: Int = 40): Boolean {
        return variantType is Inversion && variantType.length <= maxInversionLength && context.homologyLength() > maxHomLength
    }

    fun shortSplitReadTumorFilter(): Boolean {
        return isShort && tumorGenotype.splitRead() == 0
    }

    fun shortSplitReadNormalFilter(): Boolean {
        return isShort && normalGenotype.splitRead() > 0
    }

    fun longDPSupportFilter(): Boolean {
        return !isSingle && !isShort
                && normalGenotype.readPairs() == 0
                && normalGenotype.assemblyReadPairs() == 0
                && tumorGenotype.readPairs() == 0
                && tumorGenotype.assemblyReadPairs() == 0
    }

    fun impreciseFilter(): Boolean {
        return context.imprecise();
    }

    fun strandBiasFilter(maxShortStrandBias: Double): Boolean {
        return isShort && context.strandBias() > maxShortStrandBias
    }

    fun allelicFrequencyFilter(minTumorAf: Double): Boolean {
        return tumorAF < minTumorAf
    }

    fun normalCoverageFilter(minNormalCoverage: Int): Boolean {

        val supportingFragments = normalGenotype.fragmentSupport(isSingle)
        val ref = normalGenotype.refSupportRead()
        val refPair = normalGenotype.refSupportReadPair()

        return supportingFragments + ref + refPair < minNormalCoverage
    }

    fun normalSupportFilter(maxNormalSupport: Double): Boolean {
        val normalSupport = normalGenotype.fragmentSupport(isSingle)
        val tumorSupport = tumorGenotype.fragmentSupport(isSingle)

        return normalSupport > maxNormalSupport * tumorSupport;
    }

    override fun toString(): String {
        return "${context.id} ${context.contig}:${context.start} QUAL:${context.phredScaledQual} Orientation:${variantType.startOrientation} ${context.alleles[0].displayString}  > ${context.alleles[1].displayString}"
    }
}
package com.hartwig.hmftools.gripss

import com.hartwig.hmftools.extensions.*
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder

const val SHORT_EVENT_SIZE = 1000

class StructuralVariantContext(private val context: VariantContext, normalOrdinal: Int = 0, tumorOrdinal: Int = 1) {
    private companion object Constants {
        private val polyG = "G".repeat(16);
        private val polyC = "C".repeat(16);
    }

    val variantType = context.toVariantType()
    val isSingle = variantType is Single
    val isShortDup = variantType is Duplication && variantType.length < SHORT_EVENT_SIZE
    val isShortDel = variantType is Deletion && variantType.length < SHORT_EVENT_SIZE
    val isShortIns = variantType is Insertion && variantType.length < SHORT_EVENT_SIZE
    val isShort = isShortDup || isShortDel || isShortIns

    private val normalGenotype = context.getGenotype(normalOrdinal);
    private val tumorGenotype = context.getGenotype(tumorOrdinal);
    private val mateId: String? = context.mateId();
    private val tumorAF = tumorGenotype.allelicFrequency(isSingle, isShort)

    fun context(config: GripssFilterConfig): VariantContext {
        val builder = VariantContextBuilder(context).filters()

        if (normalCoverageFilter(config.minNormalCoverage)) {
            builder.filter(MIN_NORMAL_COVERAGE)
        }

        if (allelicFrequencyFilter(config.minTumorAF)) {
            builder.filter(MIN_TUMOR_AF)
        }

        if (strandBiasFilter(config.maxShortStrandBias)) {
            builder.filter(SHORT_STRAND_BIAS)
        }

        if (qualFilter(config.minQualBreakEnd, config.minQualBreakPoint)) {
            builder.filter(MIN_QUAL)
        }

        if (impreciseFilter()) {
            builder.filter(IMPRECISE)
        }

        if (polyGCFilter()) {
            builder.filter(MAX_POLY_G_LENGTH)
        }

        if (homologyLengthFilter(config.maxHomLength)) {
            builder.filter(MAX_HOM_LENGTH)
        }

        if (homologyLengthFilterShortInversion(config.maxHomLengthShortInversion)) {
            builder.filter(MAX_HOM_LENGTH_SHORT_INV)
        }

        if (inexactHomologyLengthFilter(config.maxInexactHomLength)) {
            builder.filter(MAX_INEXACT_HOM_LENGTH)
        }

        if (inexactHomologyLengthShortDelFilter(config.maxInexactHomLengthShortDel)) {
            builder.filter(MAX_INEXACT_HOM_LENGTH_SHORT_DEL)
        }

        if (shortSplitReadTumorFilter()) {
            builder.filter(SHORT_SR_SUPPORT)
        }

        if (shortSplitReadNormalFilter()) {
            builder.filter(SHORT_SR_NORMAL)
        }

        if (longDPSupportFilter()) {
            builder.filter(LONG_DP_SUPPORT)
        }

        return builder.attribute(TAF, tumorAF).make()
    }

    fun isHardFilter(config: GripssFilterConfig) = normalSupportFilter(config.maxNormalSupport)

    fun qualFilter(minQualBreakEnd: Int, minQualBreakPoint: Int): Boolean {
        val minQual = if (isSingle) minQualBreakEnd else minQualBreakPoint
        return context.phredScaledQual < minQual.toDouble()
    }

    fun polyGCFilter(): Boolean {
        return isSingle && variantType.insertSequence.contains(polyG) || variantType.insertSequence.contains(polyC)
    }

    fun inexactHomologyLengthShortDelFilter(maxInexactHomLength: Int, minDelLength: Int = 100, maxDelLength: Int = 800): Boolean {
        return variantType is Deletion && variantType.length >= minDelLength && variantType.length <= maxDelLength && context.inexactHomologyLength() > maxInexactHomLength;
    }

    fun homologyLengthFilter(maxHomLength: Int): Boolean {
        return !isSingle && context.homologyLength() > maxHomLength
    }

    fun homologyLengthFilterShortInversion(maxHomLength: Int, maxInversionLength: Int = 40): Boolean {
        return variantType is Inversion && variantType.length <= maxInversionLength && context.homologyLength() > maxHomLength
    }

    fun inexactHomologyLengthFilter(maxInexactHomLength: Int): Boolean {
        return !isSingle && !isShortDup && context.inexactHomologyLength() > maxInexactHomLength
    }

    fun shortSplitReadTumorFilter(): Boolean {
        return isShort && tumorGenotype.splitRead() == 0
    }

    fun shortSplitReadNormalFilter(): Boolean {
        return isShort && normalGenotype.splitRead() > 0
    }

    fun longDPSupportFilter(): Boolean {
        return !isShort
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

        return supportingFragments + ref + refPair > minNormalCoverage
    }

    fun normalSupportFilter(maxNormalSupport: Double): Boolean {
        val normalSupport = normalGenotype.fragmentSupport(isSingle)
        val tumorSupport = tumorGenotype.fragmentSupport(isSingle)

        return normalSupport > maxNormalSupport * tumorSupport;
    }


    private fun VariantContext.mateId(): String? {
        if (this.hasAttribute("MATE_ID")) {
            return this.getAttributeAsString("MATE_ID", "")
        }

        if (this.hasAttribute("PAR_ID")) {
            return this.getAttributeAsString("PAR_ID", "")
        }

        return null
    }
}
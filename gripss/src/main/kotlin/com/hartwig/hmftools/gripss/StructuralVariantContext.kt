package com.hartwig.hmftools.gripss

import com.hartwig.hmftools.extensions.*
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import kotlin.math.abs

const val SHORT_EVENT_SIZE = 1000

class StructuralVariantContext(private val context: VariantContext, normalOrdinal: Int = 0, tumorOrdinal: Int = 1) {

    private val POLY_G = "G".repeat(16);
    private val POLY_C = "C".repeat(16);

    val breakJunction: BreakJunction = BreakJunction.create(context.alleles[0].displayString, context.alleles[1].displayString)
    val isShortDelDup = breakJunction is BreakPoint
            && context.contig == breakJunction.chromosome
            && breakJunction.startOrientation != breakJunction.endOrientation
            && abs(context.start - breakJunction.position) < SHORT_EVENT_SIZE

    val isShortDup = breakJunction is BreakPoint && isShortDelDup &&
            (context.start <= breakJunction.position && breakJunction.startOrientation == (-1).toByte()
                    || context.start >= breakJunction.position && breakJunction.startOrientation == 1.toByte())


    private val normalGenotype = context.getGenotype(normalOrdinal);
    private val tumorGenotype = context.getGenotype(tumorOrdinal);
    private val mateId: String? = context.mateId();
    private val isSingleBreakend = breakJunction is BreakEnd
    private val tumorAF = tumorGenotype.allelicFrequency(isSingleBreakend, isShortDelDup)

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

        if (inexactHomologyLengthFilter(config.maxInexactHomLength)) {
            builder.filter(MAX_INEXACT_HOM_LENGTH)
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
        val minQual = if (isSingleBreakend) minQualBreakEnd else minQualBreakPoint
        return context.phredScaledQual < minQual.toDouble()
    }

    fun polyGCFilter(): Boolean {
        return breakJunction is BreakEnd && breakJunction.insertSequence.contains(POLY_G) || breakJunction.insertSequence.contains(POLY_C)
    }

    fun homologyLengthFilter(maxHomLength: Int): Boolean {
        return breakJunction is BreakPoint && context.homologyLength() > maxHomLength
    }

    fun inexactHomologyLengthFilter(maxInexactHomLength: Int): Boolean {
        return breakJunction is BreakPoint && !isShortDup && context.inexactHomologyLength() > maxInexactHomLength
    }

    fun shortSplitReadTumorFilter(): Boolean {
        return isShortDelDup && tumorGenotype.splitRead() == 0
    }

    fun shortSplitReadNormalFilter(): Boolean {
        return isShortDelDup && normalGenotype.splitRead() > 0
    }

    fun longDPSupportFilter(): Boolean {
        return !isShortDelDup
                && normalGenotype.readPairs() == 0
                && normalGenotype.assemblyReadPairs() == 0
                && tumorGenotype.readPairs() == 0
                && tumorGenotype.assemblyReadPairs() == 0
    }

    fun impreciseFilter(): Boolean {
        return context.imprecise();
    }

    fun strandBiasFilter(maxShortStrandBias: Double): Boolean {
        return isShortDelDup && context.strandBias() > maxShortStrandBias
    }

    fun allelicFrequencyFilter(minTumorAf: Double): Boolean {
        return tumorAF < minTumorAf
    }

    fun normalCoverageFilter(minNormalCoverage: Int): Boolean {

        val supportingFragments = normalGenotype.fragmentSupport(isSingleBreakend)
        val ref = normalGenotype.refSupportRead()
        val refPair = normalGenotype.refSupportReadPair()

        return supportingFragments + ref + refPair > minNormalCoverage
    }

    fun normalSupportFilter(maxNormalSupport: Double): Boolean {
        val normalSupport = normalGenotype.fragmentSupport(isSingleBreakend)
        val tumorSupport = tumorGenotype.fragmentSupport(isSingleBreakend)

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
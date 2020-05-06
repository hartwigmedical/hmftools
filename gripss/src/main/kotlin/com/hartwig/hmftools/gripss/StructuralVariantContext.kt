package com.hartwig.hmftools.gripss

import com.hartwig.hmftools.extensions.fragmentSupport
import com.hartwig.hmftools.extensions.refSupportRead
import com.hartwig.hmftools.extensions.refSupportReadPair
import com.hartwig.hmftools.extensions.requireAttributeAsInt
import htsjdk.variant.variantcontext.Genotype
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder

class StructuralVariantContext(normalOrdinal: Int, tumorOrdinal: Int, context: VariantContext) {
    constructor(context: VariantContext) : this(0, 1, context)

    private var context = VariantContextBuilder(context).filters().make()
    private val mateId: String? = context.mateId();
    private val isSingleBreakend = mateId == null;
    private val normalGenotype = context.getGenotype(normalOrdinal);
    private val tumorGenotype = context.getGenotype(tumorOrdinal);

    fun context(): VariantContext = context

    fun isHardFilter(config: GripssFilterConfig) = normalSupportFilter(config.maxNormalSupport)

    fun filter(config: GripssFilterConfig) {
        if (normalSupportFilter(config.maxNormalSupport)) {
            context.commonInfo.addFilter(MAX_NORMAL_SUPPORT);
        }

        if (normalCoverageFilter(config.minNormalCoverage)) {
            context.commonInfo.addFilter(MIN_NORMAL_COVERAGE);
        }
    }

    fun isPassing(): Boolean = context.filters.isEmpty() || context.filters.size == 1 && context.filters.contains("PASS")

// gridss_bp_af = function(gr, vcf, ordinal) {
//  return(.gridss_af(gr, vcf, ordinal, !is_short_deldup(gr), includeBreakpointSupport=TRUE, includeBreakendSupport=FALSE))
//}
//gridss_be_af = function(gr, vcf, ordinal) {
//  return(.gridss_af(gr, vcf, ordinal, includeRefPair=rep(TRUE, length(gr)), includeBreakpointSupport=FALSE, includeBreakendSupport=TRUE))
//}
//.gridss_af = function(gr, vcf, ordinal, includeRefPair, no_coverage_af=0, includeBreakpointSupport=TRUE, includeBreakendSupport=FALSE) {
//  assertthat::are_equal(length(gr), length(includeRefPair))
//  genotype = geno(vcf[names(gr)])
//  g = lapply(names(genotype), function(field) { if (is.numeric(genotype[[field]])) { .genosum(genotype[[field]], ordinal) } else { genotype[[field]] } })
//  names(g) <- names(genotype)
//  support = rep(0, length(gr))
//  if (includeBreakpointSupport) {
//    support = support + g$VF
//  }
//  if (includeBreakendSupport) {
//    support = support + g$BVF
//  }
//  vf_af = support / (support + g$REF + ifelse(!includeRefPair, 0, g$REFPAIR))
//  vf_af[is.nan(vf_af)] = no_coverage_af
//  return(vf_af)
//}

    fun allelicFrequency(): Double {

        val fragmentSupport = tumorGenotype.fragmentSupport(isSingleBreakend)


        return 0.0;
    }


    fun normalCoverageFilter(minNormalCoverage: Int): Boolean {

        val supportingFragments = normalGenotype.fragmentSupport(isSingleBreakend)
        val ref = normalGenotype.refSupportRead()
        val refPair = normalGenotype.refSupportReadPair()

        return supportingFragments + ref + refPair > minNormalCoverage
    }

    fun normalSupportFilter(maxNormalSupport: Double): Boolean {
        val fragment = if (isSingleBreakend) { x: Genotype -> x.requireAttributeAsInt("BVF") } else { x: Genotype -> x.requireAttributeAsInt("VF") };
        val normalSupport = fragment(normalGenotype)
        val tumorSupport = fragment(tumorGenotype)

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
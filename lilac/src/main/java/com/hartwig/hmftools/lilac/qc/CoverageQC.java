package com.hartwig.hmftools.lilac.qc;

import com.hartwig.hmftools.lilac.coverage.HlaComplexCoverage;

public class CoverageQC {

    public final int ATypes;
    public final int BTypes;
    public final int CTypes;
    public final int TotalFragments;
    public final int UniqueFragments;
    public final int SharedFragments;
    public final int WildcardFragments;

    public CoverageQC(final int ATypes, final int BTypes, final int CTypes, final int totalFragments, final int uniqueFragments,
            final int sharedFragments, final int wildcardFragments) {
        this.ATypes = ATypes;
        this.BTypes = BTypes;
        this.CTypes = CTypes;
        TotalFragments = totalFragments;
        UniqueFragments = uniqueFragments;
        SharedFragments = sharedFragments;
        WildcardFragments = wildcardFragments;
    }

    // TODO
    public double getPercentWildcard() { return 0; }

    public static CoverageQC create(int totalFragments, HlaComplexCoverage winner) {
        return null;

        /*
        val alleles = winner.AlleleCoverage.map { it.allele }
        val aTypes = alleles.filter { it.gene == "A" }.distinct().size
        val bTypes = alleles.filter { it.gene == "B" }.distinct().size
        val cTypes = alleles.filter { it.gene == "C" }.distinct().size

        if (aTypes == 0 || bTypes == 0 || cTypes == 0) {
            logger.warn("    UNMATCHED_TYPE - $aTypes A alleles, $bTypes B alleles, $cTypes C alleles")
        }

        if (winner.WildCoverage > 0) {
            logger.warn("    WILDCARD_MATCH - winning solution contains wildcards")
        }

        return CoverageQC(aTypes, bTypes, cTypes, totalFragments, winner.UniqueCoverage, winner.SharedCoverage, winner.WildCoverage)

         */
    }

    /*
    private val percentFormatter = DecimalFormat("00.0%")
    private fun Double.percent(): String {
        return percentFormatter.format(this)
    }

    fun header(): List<String> {
        return listOf("aTypes", "bTypes", "cTypes", "unusedFragments", "fittedFragments", "percentUnique", "percentShared", "percentWildcard")
        }

        fun body(): List<String> {
        return listOf(aTypes.toString(), bTypes.toString(), cTypes.toString(), unusedFragments.toString(), fittedFragments.toString(), percentUnique.percent(), percentShared.percent(), percentWildcard.percent())
        }

        override fun toString(): String {
        return "CoverageQC(aTypes=$aTypes, bTypes=$bTypes, cTypes=$cTypes, unusedFragments=$unusedFragments, uniqueFragments=$uniqueFragments, sharedFragments=$sharedFragments, wildcardFragments=$wildcardFragments, fittedFragments=$fittedFragments, percentUnique=${percentUnique.percent()}, percentShared=${percentShared.percent()}, percentWildcard=${percentWildcard.percent()})"
        }

        }

     */

}
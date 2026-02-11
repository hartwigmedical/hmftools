package com.hartwig.hmftools.compar.lilac;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_DISC_ALIGN_FRAGS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_DISC_INDELS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_HLA_Y;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_QC_STATUS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_TOTAL_FRAGS;
import static com.hartwig.hmftools.common.hla.LilacQcData.FLD_FIT_FRAGS;
import static com.hartwig.hmftools.compar.common.CategoryType.LILAC;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.hla.LilacQcData;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.driver.DriverData;

public class LilacData implements ComparableItem
{
    public final LilacQcData QcData;
    public final List<LilacAllele> Alleles;

    protected static final String FLD_ALLELES = "Alleles";
    protected static final String FLD_VARIANTS = "SomaticVariants";

    private static final String ALLELE_DELIM = ":";

    public LilacData(final LilacQcData qcData, final List<LilacAllele> alleles)
    {
        QcData = qcData;
        Alleles = alleles;
    }

    @Override
    public CategoryType category() { return LILAC; }

    @Override
    public String key()
    {
        return QcData.genes();
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%s", QcData.status()));

        StringJoiner alleles = new StringJoiner(ALLELE_DELIM);
        Alleles.forEach(x -> alleles.add(x.allele()));
        values.add(format("%s", alleles));

        values.add(format("%d", somaticVariantCount()));
        return values;
    }

    @Override
    public boolean reportable() { return true; }

    @Override
    public boolean isPass() {
        return true;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final LilacData otherData = (LilacData)other;

        if(!QcData.genes().equals(otherData.QcData.genes()))
            return false;

        return true;
    }

    protected int somaticVariantCount()
    {
        return (int)Alleles.stream().mapToDouble(x -> x.somaticVariantCount()).sum();
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches)
    {
        final LilacData otherData = (LilacData)other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_QC_STATUS, QcData.status(), otherData.QcData.status());
        checkDiff(diffs, FLD_TOTAL_FRAGS, QcData.totalFragments(), otherData.QcData.totalFragments(), thresholds);
        checkDiff(diffs, FLD_FIT_FRAGS, QcData.fittedFragments(), otherData.QcData.fittedFragments(), thresholds);
        checkDiff(diffs, FLD_DISC_ALIGN_FRAGS, QcData.discardedAlignmentFragments(), otherData.QcData.discardedAlignmentFragments(), thresholds);
        checkDiff(diffs, FLD_DISC_INDELS, QcData.discardedIndels(), otherData.QcData.discardedIndels(), thresholds);
        checkDiff(diffs, FLD_HLA_Y, QcData.hlaYAllele(), otherData.QcData.hlaYAllele());

        List<LilacAllele> origDiffs = Alleles.stream().filter(x -> !hasAllele(x, otherData.Alleles)).collect(Collectors.toList());
        List<LilacAllele> newDiffs = otherData.Alleles.stream().filter(x -> !hasAllele(x, Alleles)).collect(Collectors.toList());

        if(!origDiffs.isEmpty() || !newDiffs.isEmpty())
        {
            StringJoiner origDiffsSj = new StringJoiner(ALLELE_DELIM);
            origDiffs.forEach(x -> origDiffsSj.add(x.allele()));
            StringJoiner newDiffsSj = new StringJoiner(ALLELE_DELIM);
            newDiffs.forEach(x -> newDiffsSj.add(x.allele()));

            diffs.add(String.format("%s(%s/%s)", FLD_ALLELES, origDiffsSj, newDiffsSj));
        }

        // matches alleles in order when an allele is homozygous
        List<LilacAllele> newAllelesToMatch = Lists.newArrayList(otherData.Alleles);
        for(LilacAllele refAllele : Alleles)
        {
            LilacAllele matchingNewAllele =
                    newAllelesToMatch.stream().filter(x -> x.allele().equals(refAllele.allele())).findFirst().orElse(null);
            if(matchingNewAllele != null)
            {
                List<String> temporaryDiffs = Lists.newArrayList();
                checkDiff(temporaryDiffs, LilacAllele.FLD_MISSENSE, refAllele.somaticMissense(), matchingNewAllele.somaticMissense(), thresholds);
                checkDiff(temporaryDiffs, LilacAllele.FLD_NFS, refAllele.somaticNonsenseOrFrameshift(),
                        matchingNewAllele.somaticNonsenseOrFrameshift(), thresholds);
                checkDiff(temporaryDiffs, LilacAllele.FLD_SPLICE, refAllele.somaticSplice(), matchingNewAllele.somaticSplice(), thresholds);
                checkDiff(temporaryDiffs, LilacAllele.FLD_INDEL, refAllele.somaticInframeIndel(), matchingNewAllele.somaticInframeIndel(), thresholds);
                checkDiff(temporaryDiffs, LilacAllele.FLD_TUMOR_CN, refAllele.tumorCopyNumber(), matchingNewAllele.tumorCopyNumber(), thresholds);
                if(matchLevel == MatchLevel.DETAILED)
                {
                    checkDiff(temporaryDiffs, LilacAllele.FLD_REF_TOTAL, refAllele.refFragments(), matchingNewAllele.refFragments(), thresholds);
                    checkDiff(temporaryDiffs, LilacAllele.FLD_TUMOR_TOTAL, refAllele.tumorFragments(), matchingNewAllele.tumorFragments(), thresholds);
                    checkDiff(temporaryDiffs, LilacAllele.FLD_SYNON, refAllele.somaticSynonymous(), matchingNewAllele.somaticSynonymous(), thresholds);

                }

                temporaryDiffs.stream().map(d -> refAllele.allele() + ALLELE_DELIM + d).forEach(d -> diffs.add(d));
                newAllelesToMatch.remove(matchingNewAllele);
            }
        }

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }

    private static boolean hasAllele(final LilacAllele allele, final List<LilacAllele> alleles)
    {
        return alleles.stream().anyMatch(x -> x.allele().equals(allele.allele()));
    }
}

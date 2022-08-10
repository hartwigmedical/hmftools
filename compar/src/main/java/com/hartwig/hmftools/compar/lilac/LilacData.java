package com.hartwig.hmftools.compar.lilac;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.Category.LILAC;
import static com.hartwig.hmftools.compar.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.MismatchType.VALUE;

import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.hla.LilacSummaryData;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.DiffThresholds;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;

public class LilacData implements ComparableItem
{
    public final LilacSummaryData Data;

    protected static final String FLD_STATUS = "Status";
    protected static final String FLD_ALLELES = "Alleles";
    protected static final String FLD_VARIANTS = "SomaticVariants";
    private static final String ALLELE_DELIM = ":";

    public LilacData(final LilacSummaryData data)
    {
        Data = data;
    }

    @Override
    public Category category() { return LILAC; }

    @Override
    public String key()
    {
        return "";
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%s", Data.qc()));

        StringJoiner alleles = new StringJoiner(ALLELE_DELIM);
        Data.alleles().forEach(x -> alleles.add(x.allele()));
        values.add(format("%s", alleles.toString()));

        values.add(format("%d", Data.somaticVariantCount()));
        return values;
    }

    @Override
    public boolean reportable() { return true; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        // a single record for each sample
        return true;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds)
    {
        final LilacData otherData = (LilacData)other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_STATUS, Data.qc(), otherData.Data.qc());
        checkDiff(diffs, FLD_VARIANTS, Data.somaticVariantCount(), otherData.Data.somaticVariantCount());

        List<LilacAllele> origDiffs = Data.alleles().stream().filter(x -> !hasAllele(x, otherData.Data.alleles())).collect(Collectors.toList());
        List<LilacAllele> newDiffs = otherData.Data.alleles().stream().filter(x -> !hasAllele(x, Data.alleles())).collect(Collectors.toList());

        if(!origDiffs.isEmpty() || !newDiffs.isEmpty())
        {
            StringJoiner origDiffsSj = new StringJoiner(ALLELE_DELIM);
            origDiffs.forEach(x -> origDiffsSj.add(x.allele()));
            StringJoiner newDiffsSj = new StringJoiner(ALLELE_DELIM);
            newDiffs.forEach(x -> newDiffsSj.add(x.allele()));

            diffs.add(String.format("%s(%s/%s)", FLD_ALLELES, origDiffsSj, newDiffsSj.toString()));
        }

        return !diffs.isEmpty() ? new Mismatch(this, other, VALUE, diffs) : null;
    }

    private static boolean hasAllele(final LilacAllele allele, final List<LilacAllele> alleles)
    {
        return alleles.stream().anyMatch(x -> x.allele().equals(allele.allele()));
    }
}

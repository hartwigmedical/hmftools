package com.hartwig.hmftools.compar.peach;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.PEACH;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class PeachData implements ComparableItem
{
    public final PeachGenotype Genotype;

    protected static final String FLD_ALLELE_COUNT = "AlleleCount";
    protected static final String FLD_FUNCTION = "Function";
    protected static final String FLD_DRUGS = "Drugs";
    protected static final String FLD_PRESCRIPTION_URLS = "PrescriptionUrls";

    public PeachData(final PeachGenotype genotype)
    {
        Genotype = genotype;
    }

    @Override
    public CategoryType category()
    {
        return PEACH;
    }

    @Override
    public String key()
    {
        return Genotype.gene() + " " + Genotype.allele();
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%d", Genotype.alleleCount()));
        values.add(format("%s", Genotype.function()));
        values.add(format("%s", Genotype.linkedDrugs()));
        values.add(format("%s", Genotype.urlPrescriptionInfo()));
        return values;
    }

    @Override
    public boolean reportable()
    {
        return true;
    }

    @Override
    public boolean isPass() {
        return true;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final PeachData otherData = (PeachData) other;
        if(!Genotype.gene().equals(otherData.Genotype.gene()))
        {
            return false;
        }
        return Genotype.allele().equals(otherData.Genotype.allele());
    }

    @Override
    public Mismatch findMismatch(
            final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds, final boolean includeMatches)
    {
        final PeachData otherData = (PeachData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_ALLELE_COUNT, Genotype.alleleCount(), otherData.Genotype.alleleCount());
        checkDiff(diffs, FLD_FUNCTION, Genotype.function(), otherData.Genotype.function());
        checkDiff(diffs, FLD_DRUGS, Genotype.linkedDrugs(), otherData.Genotype.linkedDrugs());
        checkDiff(diffs, FLD_PRESCRIPTION_URLS, Genotype.urlPrescriptionInfo(), otherData.Genotype.urlPrescriptionInfo());

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }
}

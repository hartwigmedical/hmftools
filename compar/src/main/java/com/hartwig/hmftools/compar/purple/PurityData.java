package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.Category.PURITY;
import static com.hartwig.hmftools.compar.CommonUtils.checkDiff;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;

public class PurityData implements ComparableItem
{
    public final PurityContext Purity;

    public PurityData(final PurityContext purityContext) {
        Purity = purityContext;
    }

    public Category category() {
        return PURITY;
    }

    @Override
    public String key()
    {
        return "";
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(String.format("Pass(%s)", Purity.qc().pass()));
        values.add(String.format("Purity(%.2f)", Purity.bestFit().purity()));
        values.add(String.format("Ploidy(%.2f)", Purity.bestFit().ploidy()));
        return values;
    }

    @Override
    public boolean reportable() {
        return true;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        // a single record for each sample
        return true;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel)
    {
        final PurityData otherPurity = (PurityData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, "fitMethod", Purity.method().toString(), otherPurity.Purity.method().toString());
        checkDiff(diffs, "qcPass", Purity.qc().pass(), otherPurity.Purity.qc().pass());

        if (matchLevel == REPORTABLE)
            return null;

        checkDiff(diffs, "purity", Purity.bestFit().purity(), otherPurity.Purity.bestFit().purity());
        checkDiff(diffs, "ploidy", Purity.bestFit().ploidy(), otherPurity.Purity.bestFit().ploidy());

        return null;
    }
}

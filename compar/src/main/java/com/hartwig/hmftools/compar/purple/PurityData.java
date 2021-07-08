package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.Category.DRIVER;
import static com.hartwig.hmftools.compar.Category.PURITY;
import static com.hartwig.hmftools.compar.CommonUtils.ITEM_DELIM;
import static com.hartwig.hmftools.compar.CommonUtils.checkDiff;
import static com.hartwig.hmftools.compar.CommonUtils.diffValue;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;

import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriver;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.MatchLevel;

import org.apache.commons.compress.utils.Lists;

public class PurityData implements ComparableItem
{
    public final PurityContext Purity;

    public PurityData(final PurityContext purityContext)
    {
        Purity = purityContext;
    }

    public Category category() { return PURITY; }

    public boolean reportable() { return true; }

    public boolean matches(final ComparableItem other)
    {
        // a single record for each sample
        return true;
    }

    public List<String> findDifferences(final ComparableItem other, final MatchLevel matchLevel)
    {
        final PurityData otherPurity = (PurityData)other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, "fitMethod", Purity.method().toString(), otherPurity.Purity.method().toString());
        checkDiff(diffs, "qcPass", Purity.qc().pass(), otherPurity.Purity.qc().pass());

        if(matchLevel == REPORTABLE)
            return diffs;

        checkDiff(diffs, "purity", Purity.bestFit().purity(), otherPurity.Purity.bestFit().purity());
        checkDiff(diffs, "ploidy", Purity.bestFit().ploidy(), otherPurity.Purity.bestFit().ploidy());

        return diffs;
    }

    public String description()
    {
        return String.format("%s_%.2f_%.2f", Purity.qc().pass(), Purity.bestFit().purity(), Purity.bestFit().ploidy());
    }
}

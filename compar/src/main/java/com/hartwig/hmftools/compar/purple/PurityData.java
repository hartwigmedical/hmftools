package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.Category.PURITY;
import static com.hartwig.hmftools.compar.CommonUtils.checkDiff;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityQCContext;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.MatchLevel;

public class PurityData implements ComparableItem {
    public final PurityQCContext Purity;

    public PurityData(final PurityQCContext purityQCContext) {
        Purity = purityQCContext;
    }

    public Category category() {
        return PURITY;
    }

    public boolean reportable() {
        return true;
    }

    public boolean matches(final ComparableItem other) {
        // a single record for each sample
        return true;
    }

    public List<String> findDifferences(final ComparableItem other, final MatchLevel matchLevel) {
        final PurityData otherPurity = (PurityData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, "fitMethod", Purity.purityContext().method().toString(), otherPurity.Purity.purityContext().method().toString());
        checkDiff(diffs, "qcPass", Purity.qc().pass(), otherPurity.Purity.qc().pass());

        if (matchLevel == REPORTABLE) {
            return diffs;
        }

        checkDiff(diffs, "purity", Purity.purityContext().bestFit().purity(), otherPurity.Purity.purityContext().bestFit().purity());
        checkDiff(diffs, "ploidy", Purity.purityContext().bestFit().ploidy(), otherPurity.Purity.purityContext().bestFit().ploidy());

        return diffs;
    }

    public String description() {
        return String.format("%s_%.2f_%.2f",
                Purity.qc().pass(),
                Purity.purityContext().bestFit().purity(),
                Purity.purityContext().bestFit().ploidy());
    }
}

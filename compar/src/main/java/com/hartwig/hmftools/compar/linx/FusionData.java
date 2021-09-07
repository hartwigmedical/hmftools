package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.Category.FUSION;
import static com.hartwig.hmftools.compar.CommonUtils.checkDiff;
import static com.hartwig.hmftools.compar.CommonUtils.diffValue;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;

import java.util.List;

import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.MatchLevel;

import org.apache.commons.compress.utils.Lists;

public class FusionData implements ComparableItem
{
    public final LinxFusion Fusion;

    public FusionData(final LinxFusion fusion)
    {
        Fusion = fusion;
    }

    public Category category() { return FUSION; }

    public boolean reportable() { return Fusion.reported(); }

    public boolean matches(final ComparableItem other)
    {
        final FusionData otherFusion = (FusionData)other;

        return otherFusion.Fusion.name().equals(Fusion.name());
    }

    public List<String> findDifferences(final ComparableItem other, final MatchLevel matchLevel)
    {
        final FusionData otherFusion = (FusionData)other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, "reported", Fusion.reported(), otherFusion.Fusion.reported());
        checkDiff(diffs, "reportedType", Fusion.reportedType(), otherFusion.Fusion.reportedType());

        if(matchLevel == REPORTABLE)
            return diffs;

        checkDiff(diffs, "phased", Fusion.phased().toString(), otherFusion.Fusion.phased().toString());
        checkDiff(diffs, "likelihood", Fusion.likelihood().toString(), otherFusion.Fusion.likelihood().toString());
        checkDiff(diffs, "fusedExonsUp", Fusion.fusedExonUp(), otherFusion.Fusion.fusedExonUp());
        checkDiff(diffs, "fusedExonsDown", Fusion.fusedExonDown(), otherFusion.Fusion.fusedExonDown());

        /*
            public abstract int chainLength();
            public abstract int chainLinks();
            public abstract boolean chainTerminated();
            public abstract String domainsKept();
            public abstract String domainsLost();
            public abstract int skippedExonsUp();
            public abstract int skippedExonsDown();

         */

        return diffs;
    }

    public String description()
    {
        return String.format("%s", Fusion.name());
    }
}

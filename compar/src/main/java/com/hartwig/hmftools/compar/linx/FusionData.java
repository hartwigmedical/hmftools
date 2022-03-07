package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.Category.FUSION;
import static com.hartwig.hmftools.compar.CommonUtils.checkDiff;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;

public class FusionData implements ComparableItem
{
    public final LinxFusion Fusion;
    public final String GeneMappedName;

    public FusionData(final LinxFusion fusion, final String geneMappedName)
    {
        Fusion = fusion;
        GeneMappedName = geneMappedName;
    }

    @Override
    public Category category() { return FUSION; }

    @Override
    public String key()
    {
        return String.format("%s_%s", Fusion.name(), Fusion.reportedType());
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
//        values.add(String.format("Qual(%.0f)", Variant.qual()));
//        values.add(String.format("Tier(%s)", Variant.tier().toString()));
//        values.add(String.format("TotalReadCount(%d)", Variant.totalReadCount()));
//        values.add(String.format("AlleleReadCount(%d)", Variant.alleleReadCount()));
        return values;
    }

    @Override
    public boolean reportable() { return Fusion.reported(); }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final FusionData otherFusion = (FusionData)other;

        return otherFusion.GeneMappedName.equals(GeneMappedName);
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel)
    {
        return null;
    }

    @Override
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
}

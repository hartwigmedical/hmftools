package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.Category.FUSION;
import static com.hartwig.hmftools.compar.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.compar.MismatchType.VALUE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.DiffThresholds;
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
        values.add(String.format("reportedType=%s", Fusion.reportedType()));
        values.add(String.format("phased=%s", Fusion.phased()));
        values.add(String.format("likelihood=%s", Fusion.likelihood()));
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
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds)
    {
        final FusionData otherFusion = (FusionData)other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, "reported", Fusion.reported(), otherFusion.Fusion.reported());
        checkDiff(diffs, "reportedType", Fusion.reportedType(), otherFusion.Fusion.reportedType());

        checkDiff(diffs, "phased", Fusion.phased().toString(), otherFusion.Fusion.phased().toString());
        checkDiff(diffs, "likelihood", Fusion.likelihood().toString(), otherFusion.Fusion.likelihood().toString());
        checkDiff(diffs, "fusedExonsUp", Fusion.fusedExonUp(), otherFusion.Fusion.fusedExonUp());
        checkDiff(diffs, "fusedExonsDown", Fusion.fusedExonDown(), otherFusion.Fusion.fusedExonDown());
        checkDiff(diffs, "chainLinks", Fusion.chainLinks(), otherFusion.Fusion.chainLinks());
        checkDiff(diffs, "chainTerminated", Fusion.chainTerminated(), otherFusion.Fusion.chainTerminated());
        checkDiff(diffs, "domainsKept", Fusion.domainsKept(), otherFusion.Fusion.domainsKept());
        checkDiff(diffs, "domainsLost", Fusion.domainsLost(), otherFusion.Fusion.domainsLost());

        return !diffs.isEmpty() ? new Mismatch(this, other, VALUE, diffs) : null;
    }
}

package com.hartwig.hmftools.compar.linx;

import static com.hartwig.hmftools.compar.Category.FUSION;
import static com.hartwig.hmftools.compar.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.DiffFunctions.checkDiff;
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

    protected static final String FLD_REPORTED_TYPE = "ReportedType";
    protected static final String FLD_PHASED = "Phased";
    protected static final String FLD_LIKELIHOOD = "likelihood";
    protected static final String FLD_EXON_UP = "fusedExonUp";
    protected static final String FLD_EXON_DOWN = "fusedExonDown";
    protected static final String FLD_CHAIN_LINKS = "chainLinks";
    protected static final String FLD_CHAIN_TERM = "chainTerminated";
    protected static final String FLD_DOMAINS_KEPT = "domainsKept";
    protected static final String FLD_DOMAINS_LOST = "domainsLost";

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
        values.add(String.format("%s", Fusion.reported()));
        values.add(String.format("%s", Fusion.reportedType()));
        values.add(String.format("%s", Fusion.phased()));
        values.add(String.format("%s", Fusion.likelihood()));
        values.add(String.format("%d", Fusion.fusedExonUp()));
        values.add(String.format("%d", Fusion.fusedExonDown()));
        values.add(String.format("%d", Fusion.chainLinks()));
        values.add(String.format("%s", Fusion.chainTerminated()));
        values.add(String.format("%s", Fusion.domainsKept()));
        values.add(String.format("%s", Fusion.domainsLost()));
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

        checkDiff(diffs, FLD_REPORTED, Fusion.reported(), otherFusion.Fusion.reported());
        checkDiff(diffs, FLD_REPORTED_TYPE, Fusion.reportedType(), otherFusion.Fusion.reportedType());
            checkDiff(diffs, FLD_PHASED, Fusion.phased().toString(), otherFusion.Fusion.phased().toString());
        checkDiff(diffs, FLD_LIKELIHOOD, Fusion.likelihood().toString(), otherFusion.Fusion.likelihood().toString());
        checkDiff(diffs, FLD_EXON_UP, Fusion.fusedExonUp(), otherFusion.Fusion.fusedExonUp());
        checkDiff(diffs, FLD_EXON_DOWN, Fusion.fusedExonDown(), otherFusion.Fusion.fusedExonDown());
        checkDiff(diffs, FLD_CHAIN_LINKS, Fusion.chainLinks(), otherFusion.Fusion.chainLinks());
        checkDiff(diffs, FLD_CHAIN_TERM, Fusion.chainTerminated(), otherFusion.Fusion.chainTerminated());
        checkDiff(diffs, FLD_DOMAINS_KEPT, Fusion.domainsKept(), otherFusion.Fusion.domainsKept());
        checkDiff(diffs, FLD_DOMAINS_LOST, Fusion.domainsLost(), otherFusion.Fusion.domainsLost());

        return !diffs.isEmpty() ? new Mismatch(this, other, VALUE, diffs) : null;
    }
}

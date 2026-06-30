package com.hartwig.hmftools.compar.linx;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.FUSION;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.FieldConfig;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class FusionData implements ComparableItem
{
    public final LinxFusion Fusion;
    public final String GeneMappedName;
    public final BreakendData BreakendFive;
    public final BreakendData BreakendThree;

    protected static final String FLD_REPORTED_TYPE = "ReportedType";
    protected static final String FLD_PHASED = "Phased";
    protected static final String FLD_LIKELIHOOD = "Likelihood";
    protected static final String FLD_TRANSCRIPT_UP = "FusedTranscriptUp";
    protected static final String FLD_EXON_UP = "FusedExonUp";
    protected static final String FLD_TRANSCRIPT_DOWN = "FusedTranscriptDown";
    protected static final String FLD_EXON_DOWN = "FusedExonDown";
    protected static final String FLD_CHAIN_LINKS = "ChainLinks";
    protected static final String FLD_CHAIN_TERM = "ChainTerminated";
    protected static final String FLD_DOMAINS_KEPT = "DomainsKept";
    protected static final String FLD_DOMAINS_LOST = "DomainsLost";
    protected static final String FLD_JUNCTION_COPY_NUMBER = "JunctionCopyNumber";

    public FusionData(final LinxFusion fusion, final String geneMappedName, final BreakendData breakendFive, final BreakendData breakendThree)
    {
        Fusion = fusion;
        GeneMappedName = geneMappedName;
        BreakendFive = breakendFive;
        BreakendThree = breakendThree;
    }

    @Override
    public CategoryType category() { return FUSION; }

    @Override
    public String key()
    {
        return format("%s_%s", Fusion.name(), Fusion.reportedType());
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%s", Fusion.reported()));
        values.add(format("%s", Fusion.reportedType()));
        values.add(format("%s", Fusion.phased()));
        values.add(format("%s", Fusion.likelihood()));
        values.add(format("%d", Fusion.fusedExonUp()));
        values.add(format("%d", Fusion.fusedExonDown()));
        values.add(format("%d", Fusion.chainLinks()));
        values.add(format("%s", Fusion.chainTerminated()));
        values.add(format("%s", Fusion.domainsKept()));
        values.add(format("%s", Fusion.domainsLost()));

        return values;
    }

    @Override
    public List<String> extraInfoValues()
    {
        List<String> values = Lists.newArrayList();

        if(BreakendFive != null)
            values.add(format("BreakendUp=%s", BreakendFive.fullStr(true)));

        if(BreakendThree != null)
            values.add(format("BreakendDown=%s", BreakendThree.fullStr(true)));

        return values;
    }

    @Override
    public boolean reportable() {
        return Fusion.reported();
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final FusionData otherFusion = (FusionData)other;

        return otherFusion.GeneMappedName.equals(GeneMappedName);
    }

    @Override
    public Mismatch findMismatch(
            final ComparableItem other, final MatchLevel matchLevel, final FieldConfig fieldConfig, final boolean includeMatches)
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

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }
}

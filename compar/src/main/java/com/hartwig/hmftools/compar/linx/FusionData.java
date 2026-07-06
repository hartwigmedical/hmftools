package com.hartwig.hmftools.compar.linx;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.FUSION;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;

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
}

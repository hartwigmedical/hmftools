package com.hartwig.hmftools.compar.linx;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.FUSION;

import java.util.List;

import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class FusionData extends ComparableItem
{
    public final LinxFusion Fusion;
    public final String GeneMappedName;
    public final BreakendData BreakendFive;
    public final BreakendData BreakendThree;

    public FusionData(
            final LinxFusion fusion, final String geneMappedName, final BreakendData breakendFive, final BreakendData breakendThree,
            final List<FieldInfo> fields)
    {
        Fusion = fusion;
        GeneMappedName = geneMappedName;
        BreakendFive = breakendFive;
        BreakendThree = breakendThree;

        addBoolValue(FusionComparer.Fields.Reported.toString(), fusion.reported(), fields);
        addStringValue(FusionComparer.Fields.ReportedType.toString(), fusion.reportedType(), fields);
        addStringValue(FusionComparer.Fields.Phased.toString(), fusion.phased().toString(), fields);
        addStringValue(FusionComparer.Fields.Likelihood.toString(), fusion.likelihood().toString(), fields);

        addIntValue(FusionComparer.Fields.FusedExonUp.toString(), fusion.fusedExonUp(), fields);
        addIntValue(FusionComparer.Fields.FusedExonDown.toString(), fusion.fusedExonDown(), fields);
        addIntValue(FusionComparer.Fields.ChainLinks.toString(), fusion.chainLinks()    , fields);
        addBoolValue(FusionComparer.Fields.ChainTerminated.toString(), fusion.chainTerminated(), fields);

        addStringValue(FusionComparer.Fields.DomainsKept.toString(), fusion.domainsKept(), fields);
        addStringValue(FusionComparer.Fields.DomainsLost.toString(), fusion.domainsLost(), fields);

        // TODO: add field values so can build for display only
        // addStringValue(FusionComparer.Fields.TumorCopyNumber.toString(), fusion.tumorCopyNumber(), fields);
        // addStringValue(FusionComparer.Fields.TumorCopyNumber.toString(), fusion.tumorCopyNumber(), fields);

    }

    @Override
    public CategoryType category() { return FUSION; }

    @Override
    public String key()
    {
        return format("%s_%s", Fusion.name(), Fusion.reportedType());
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

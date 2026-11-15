package com.hartwig.hmftools.compar.purple;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.purple.ReportedStatus.REPORTED;
import static com.hartwig.hmftools.compar.common.CategoryType.GERMLINE_AMP_DEL;

import java.util.List;

import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class GermlineAmpDelData extends ComparableItem
{
    public final GermlineAmpDel AmpDelData;

    public GermlineAmpDelData(final GermlineAmpDel germlineAmpDel, final List<FieldInfo> fields)
    {
        AmpDelData = germlineAmpDel;

        addStringValue(GermlineAmpDelComparer.Fields.Reported.toString(), germlineAmpDel.Reported.toString(), fields);
        addStringValue(GermlineAmpDelComparer.Fields.GermlineStatus.toString(), germlineAmpDel.NormalStatus.toString(), fields);
        addStringValue(GermlineAmpDelComparer.Fields.TumorStatus.toString(), germlineAmpDel.TumorStatus.toString(), fields);
        addDoubleValue(GermlineAmpDelComparer.Fields.GermlineCopyNumber.toString(), germlineAmpDel.GermlineCopyNumber, fields);
        addDoubleValue(GermlineAmpDelComparer.Fields.TumorCopyNumber.toString(), germlineAmpDel.TumorCopyNumber, fields);
    }

    public CategoryType category() {
        return GERMLINE_AMP_DEL;
    }

    @Override
    public String key()
    {
        return format("%s", AmpDelData.GeneName);
    }

    @Override
    public boolean reportable() {
        return AmpDelData.Reported == REPORTED;
    }

    @Override
    public String geneName() { return AmpDelData.GeneName; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final GermlineAmpDelData otherAmpDel = (GermlineAmpDelData)other;
        return AmpDelData.GeneName.equals(otherAmpDel.AmpDelData.GeneName);
    }
}

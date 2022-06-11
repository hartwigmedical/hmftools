package com.hartwig.hmftools.compar.purple;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.Category.PURITY;
import static com.hartwig.hmftools.compar.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.MismatchType.VALUE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.gene.GermlineDeletion;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.DiffThresholds;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;

public class GermlineDeletionData implements ComparableItem
{
    public final GermlineDeletion Deletion;

    protected static final String FLD_GERMLINE_CN = "germlineCopyNumber";
    protected static final String FLD_TUMOR_CN = "tumorCopyNumber";

    public GermlineDeletionData(final GermlineDeletion germlineDeletion)
    {
        Deletion = germlineDeletion;
    }

    public Category category() {
        return PURITY;
    }

    @Override
    public String key()
    {
        return format("%s", Deletion.GeneName);
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("reported=%s", Deletion.Reported));
        values.add(format("germlineStatus=%s", Deletion.NormalStatus));
        values.add(format("tumorStatus=%s", Deletion.TumorStatus));
        return values;
    }

    @Override
    public boolean reportable() {
        return true;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final GermlineDeletionData otherDeletion = (GermlineDeletionData)other;
        return Deletion.GeneName.equals(otherDeletion.Deletion.GeneName);
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds)
    {
        final GermlineDeletionData otherDeletion = (GermlineDeletionData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, "qcStatus", Deletion.Reported, otherDeletion.Deletion.Reported);
        checkDiff(diffs, "germlineStatus", Deletion.NormalStatus.toString(), otherDeletion.Deletion.NormalStatus.toString());
        checkDiff(diffs, "tumorStatus", Deletion.TumorStatus.toString(), otherDeletion.Deletion.TumorStatus.toString());

        checkDiff(diffs, FLD_GERMLINE_CN, Deletion.GermlineCopyNumber, otherDeletion.Deletion.GermlineCopyNumber, thresholds);
        checkDiff(diffs, FLD_TUMOR_CN, Deletion.TumorCopyNumber, otherDeletion.Deletion.TumorCopyNumber, thresholds);

        return !diffs.isEmpty() ? new Mismatch(this, other, VALUE, diffs) : null;
    }
}

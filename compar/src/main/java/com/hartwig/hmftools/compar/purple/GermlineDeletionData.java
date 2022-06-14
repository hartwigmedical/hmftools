package com.hartwig.hmftools.compar.purple;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.Category.GERMLINE_DELETION;
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

    protected static final String FLD_QC_STATUS = "QcStatus";
    protected static final String FLD_GERMLINE_STATUS = "GermlineStatus";
    protected static final String FLD_TUMOR_STATUS = "TumorStatus";
    protected static final String FLD_GERMLINE_CN = "GermlineCopyNumber";
    protected static final String FLD_TUMOR_CN = "TumorCopyNumber";

    public GermlineDeletionData(final GermlineDeletion germlineDeletion)
    {
        Deletion = germlineDeletion;
    }

    public Category category() {
        return GERMLINE_DELETION;
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
        values.add(format("%s", Deletion.Reported));
        values.add(format("%s", Deletion.NormalStatus));
        values.add(format("%s", Deletion.TumorStatus));
        values.add(format("%s", Deletion.GermlineCopyNumber));
        values.add(format("%s", Deletion.TumorCopyNumber));
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

        checkDiff(diffs, FLD_QC_STATUS, Deletion.Reported, otherDeletion.Deletion.Reported);
        checkDiff(diffs, FLD_GERMLINE_STATUS, Deletion.NormalStatus.toString(), otherDeletion.Deletion.NormalStatus.toString());
        checkDiff(diffs, FLD_TUMOR_STATUS, Deletion.TumorStatus.toString(), otherDeletion.Deletion.TumorStatus.toString());
        checkDiff(diffs, FLD_GERMLINE_CN, Deletion.GermlineCopyNumber, otherDeletion.Deletion.GermlineCopyNumber, thresholds);
        checkDiff(diffs, FLD_TUMOR_CN, Deletion.TumorCopyNumber, otherDeletion.Deletion.TumorCopyNumber, thresholds);

        return !diffs.isEmpty() ? new Mismatch(this, other, VALUE, diffs) : null;
    }
}

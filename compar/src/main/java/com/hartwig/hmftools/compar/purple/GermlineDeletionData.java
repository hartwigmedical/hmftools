package com.hartwig.hmftools.compar.purple;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.purple.ReportedStatus.REPORTED;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.compar.common.Category.GERMLINE_DELETION;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_CHROMOSOME_BAND;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class GermlineDeletionData implements ComparableItem
{
    public final GermlineDeletion Deletion;
    public final String mComparisonChromosome;

    protected static final String FLD_GERMLINE_STATUS = "GermlineStatus";
    protected static final String FLD_TUMOR_STATUS = "TumorStatus";
    protected static final String FLD_GERMLINE_CN = "GermlineCopyNumber";
    protected static final String FLD_TUMOR_CN = "TumorCopyNumber";

    public GermlineDeletionData(final GermlineDeletion germlineDeletion, final String comparisonChromosome)
    {
        Deletion = germlineDeletion;
        mComparisonChromosome = comparisonChromosome;
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
        return Deletion.Reported == REPORTED;
    }

    @Override
    public boolean isPass() {
        return true;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final GermlineDeletionData otherDeletion = (GermlineDeletionData)other;
        return Deletion.GeneName.equals(otherDeletion.Deletion.GeneName);
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches)
    {
        final GermlineDeletionData otherDeletion = (GermlineDeletionData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_REPORTED, Deletion.Reported == REPORTED, otherDeletion.Deletion.Reported == REPORTED);
        checkDiff(diffs, FLD_GERMLINE_STATUS, Deletion.NormalStatus.toString(), otherDeletion.Deletion.NormalStatus.toString());
        checkDiff(diffs, FLD_TUMOR_STATUS, Deletion.TumorStatus.toString(), otherDeletion.Deletion.TumorStatus.toString());
        checkDiff(diffs, FLD_GERMLINE_CN, Deletion.GermlineCopyNumber, otherDeletion.Deletion.GermlineCopyNumber, thresholds);
        checkDiff(diffs, FLD_TUMOR_CN, Deletion.TumorCopyNumber, otherDeletion.Deletion.TumorCopyNumber, thresholds);
        checkDiff(diffs, FLD_CHROMOSOME, mComparisonChromosome, otherDeletion.mComparisonChromosome);
        checkDiff(diffs, FLD_CHROMOSOME_BAND, Deletion.ChromosomeBand, otherDeletion.Deletion.ChromosomeBand);

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }
}

package com.hartwig.hmftools.compar.purple;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.purple.ReportedStatus.REPORTED;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.compar.common.CategoryType.GERMLINE_AMP_DEL;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_CHROMOSOME_BAND;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class GermlineAmpDelData implements ComparableItem
{
    public final GermlineAmpDel AmpDelData;
    public final String mComparisonChromosome;

    protected static final String FLD_GERMLINE_STATUS = "GermlineStatus";
    protected static final String FLD_TUMOR_STATUS = "TumorStatus";
    protected static final String FLD_GERMLINE_CN = "GermlineCopyNumber";
    protected static final String FLD_TUMOR_CN = "TumorCopyNumber";

    public GermlineAmpDelData(final GermlineAmpDel germlineAmpDel, final String comparisonChromosome)
    {
        AmpDelData = germlineAmpDel;
        mComparisonChromosome = comparisonChromosome;
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
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%s", AmpDelData.Reported));
        values.add(format("%s", AmpDelData.NormalStatus));
        values.add(format("%s", AmpDelData.TumorStatus));
        values.add(format("%s", AmpDelData.GermlineCopyNumber));
        values.add(format("%s", AmpDelData.TumorCopyNumber));
        return values;
    }

    @Override
    public boolean reportable() {
        return AmpDelData.Reported == REPORTED;
    }

    @Override
    public boolean isPass() {
        return true;
    }

    @Override
    public String geneName() { return AmpDelData.GeneName; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final GermlineAmpDelData otherAmpDel = (GermlineAmpDelData)other;
        return AmpDelData.GeneName.equals(otherAmpDel.AmpDelData.GeneName);
    }

    @Override
    public Mismatch findMismatch(
            final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds, final boolean includeMatches)
    {
        final GermlineAmpDelData otherDeletion = (GermlineAmpDelData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_REPORTED, AmpDelData.Reported == REPORTED, otherDeletion.AmpDelData.Reported == REPORTED);
        checkDiff(diffs, FLD_GERMLINE_STATUS, AmpDelData.NormalStatus.toString(), otherDeletion.AmpDelData.NormalStatus.toString());
        checkDiff(diffs, FLD_TUMOR_STATUS, AmpDelData.TumorStatus.toString(), otherDeletion.AmpDelData.TumorStatus.toString());
        checkDiff(diffs, FLD_GERMLINE_CN, AmpDelData.GermlineCopyNumber, otherDeletion.AmpDelData.GermlineCopyNumber, thresholds);
        checkDiff(diffs, FLD_TUMOR_CN, AmpDelData.TumorCopyNumber, otherDeletion.AmpDelData.TumorCopyNumber, thresholds);
        checkDiff(diffs, FLD_CHROMOSOME, mComparisonChromosome, otherDeletion.mComparisonChromosome);
        checkDiff(diffs, FLD_CHROMOSOME_BAND, AmpDelData.ChromosomeBand, otherDeletion.AmpDelData.ChromosomeBand);

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }
}

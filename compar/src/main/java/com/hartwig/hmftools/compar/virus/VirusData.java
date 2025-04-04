package com.hartwig.hmftools.compar.virus;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.Category.VIRUS;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class VirusData implements ComparableItem
{
    public final AnnotatedVirus Virus;

    protected static final String FLD_INTEGRATIONS = "Integrations";
    protected static final String FLD_MEAN_COVERAGE = "MeanCoverage";
    protected static final String FLD_DRIVER_LIKELIHOOD = "DriverLikelihood";

    VirusData(final AnnotatedVirus virus)
    {
        Virus = virus;
    }

    @Override
    public Category category() {
        return VIRUS;
    }

    @Override
    public String key()
    {
        return Virus.name();
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%s", Virus.reported()));
        values.add(format("%d", Virus.integrations()));
        values.add(format("%.2f", Virus.meanCoverage()));
        values.add(format("%s", Virus.virusDriverLikelihoodType()));
        return values;
    }

    @Override
    public boolean reportable()
    {
        return Virus.reported();
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final VirusData otherData = (VirusData) other;
        return Virus.name().equals(otherData.Virus.name());
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches)
    {
        final VirusData otherData = (VirusData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_REPORTED, Virus.reported(), otherData.Virus.reported());
        checkDiff(diffs, FLD_INTEGRATIONS, Virus.integrations(), otherData.Virus.integrations(), thresholds);
        checkDiff(diffs, FLD_MEAN_COVERAGE, Virus.meanCoverage(), otherData.Virus.meanCoverage(), thresholds);
        checkDiff(diffs, FLD_DRIVER_LIKELIHOOD, String.valueOf(Virus.virusDriverLikelihoodType()), String.valueOf(otherData.Virus.virusDriverLikelihoodType()));

        return createMismatchFromDiffs(this, other, diffs, includeMatches);
    }
}

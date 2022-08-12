package com.hartwig.hmftools.compar.cuppa;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.Category.CUPPA;
import static com.hartwig.hmftools.compar.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.MismatchType.VALUE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.DiffThresholds;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;

public class CuppaData implements ComparableItem
{
    public final ClassifierData ClassifierResult;
    private final String mKey;

    protected static final String FLD_TOP_CANCER_TYPE = "TopCancerType";
    protected static final String FLD_LIKELIHOOD = "Likelihood";

    public CuppaData(final ClassifierData result)
    {
        ClassifierResult = result;
        mKey = result.DataType;
    }

    @Override
    public Category category() { return CUPPA; }

    @Override
    public String key()
    {
        return mKey;
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%s", ClassifierResult.TopRefCancerType));
        values.add(format("%.3f", ClassifierResult.TopRefCancerValue));
        return values;
    }

    @Override
    public boolean reportable() { return true; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final CuppaData otherCuppaData = (CuppaData)other;
        return otherCuppaData.ClassifierResult.DataType.equals(ClassifierResult.DataType);
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds)
    {
        final CuppaData otherCuppaData = (CuppaData)other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(
                diffs, FLD_TOP_CANCER_TYPE, ClassifierResult.TopRefCancerType, otherCuppaData.ClassifierResult.TopRefCancerType);

        checkDiff(diffs, FLD_LIKELIHOOD, ClassifierResult.TopRefCancerValue, otherCuppaData.ClassifierResult.TopRefCancerValue, thresholds);

        return !diffs.isEmpty() ? new Mismatch(this, other, VALUE, diffs) : null;
    }

    public String toString()
    {
        return format("type(%s) top(cancer=%s value=%.3f)",
            ClassifierResult.DataType, ClassifierResult.TopRefCancerType, ClassifierResult.TopRefCancerValue);
    }
}

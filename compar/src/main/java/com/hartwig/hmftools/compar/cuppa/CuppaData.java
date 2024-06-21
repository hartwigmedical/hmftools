package com.hartwig.hmftools.compar.cuppa;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.Category.CUPPA;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.common.MismatchType.VALUE;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cuppa.CuppaPredictionEntry;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class CuppaData implements ComparableItem
{
    public final CuppaPredictionEntry PredictionEntry;

    protected static final String FLD_CLASSIFIER_NAME = "classifier_name";
    protected static final String FLD_TOP_CANCER_TYPE = "top_cancer_type";
    protected static final String FLD_PROBABILITY = "probability";

    public CuppaData(final CuppaPredictionEntry predictionEntry)
    {
        PredictionEntry = predictionEntry;
    }

    @Override
    public Category category() { return CUPPA; }

    @Override
    public String key()
    {
        return String.format("%s;%s", PredictionEntry.DataType, PredictionEntry.ClassifierName);
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%s", PredictionEntry.CancerType));
        values.add(format("%.3f", PredictionEntry.DataValue));
        return values;
    }

    @Override
    public boolean reportable() { return true; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final CuppaData otherCuppaData = (CuppaData) other;

        // Match by DataType in case we want to compare other DataTypes (e.g. 'feat_contrib' and 'sig_quantile') in the future
        // Currently only support 'prob' DataType
        return otherCuppaData.PredictionEntry.DataType.equals(PredictionEntry.DataType) &
                otherCuppaData.PredictionEntry.ClassifierName.equals(PredictionEntry.ClassifierName);
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds)
    {
        final CuppaData otherCuppaData = (CuppaData)other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_TOP_CANCER_TYPE, PredictionEntry.CancerType, otherCuppaData.PredictionEntry.CancerType);
        checkDiff(diffs, FLD_PROBABILITY, PredictionEntry.DataValue, otherCuppaData.PredictionEntry.DataValue, thresholds);

        return !diffs.isEmpty() ? new Mismatch(this, other, VALUE, diffs) : null;
    }

    public String toString()
    {
        return format("type(%s) top(cancer=%s value=%.3f)",
            PredictionEntry.DataType, PredictionEntry.CancerType, PredictionEntry.DataValue);
    }
}

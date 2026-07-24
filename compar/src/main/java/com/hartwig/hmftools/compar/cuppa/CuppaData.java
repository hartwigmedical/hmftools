package com.hartwig.hmftools.compar.cuppa;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.CUPPA;

import com.hartwig.hmftools.common.cuppa.CuppaPredictionEntry;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;

public class CuppaData implements ComparableItem
{
    public final CuppaPredictionEntry PredictionEntry;

    public CuppaData(final CuppaPredictionEntry predictionEntry)
    {
        PredictionEntry = predictionEntry;
    }

    @Override
    public CategoryType category() { return CUPPA; }

    @Override
    public String key()
    {
        return String.format("%s;%s", PredictionEntry.DataType, PredictionEntry.ClassifierName);
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final CuppaData otherCuppaData = (CuppaData) other;

        // Match by DataType in case we want to compare other DataTypes (e.g. 'feat_contrib' and 'sig_quantile') in the future
        // Currently only support 'prob' DataType
        return otherCuppaData.PredictionEntry.DataType.equals(PredictionEntry.DataType) &&
                otherCuppaData.PredictionEntry.ClassifierName.equals(PredictionEntry.ClassifierName);
    }

    public String toString()
    {
        return format("type(%s) top(cancer=%s value=%.3f)",
            PredictionEntry.DataType, PredictionEntry.CancerType, PredictionEntry.DataValue);
    }
}

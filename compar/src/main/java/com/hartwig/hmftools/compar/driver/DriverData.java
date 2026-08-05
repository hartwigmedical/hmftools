package com.hartwig.hmftools.compar.driver;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.DRIVER;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.field.FieldInfo;
import com.hartwig.hmftools.compar.common.field.FieldValue;

public class DriverData extends ComparableItem
{
    public final DriverCatalog DriverCatalog;

    public final String mComparisonChromosome;
    private final String mKey;
    private final boolean mCheckTranscript;
    private final boolean mIsPass;

    public DriverData(
            final DriverCatalog driverCatalog, final PurplePurity purity, final String comparisonChromosome,
            boolean checkTranscript, boolean isPass, final List<FieldInfo> fields)
    {
        DriverCatalog = driverCatalog;

        mComparisonChromosome = comparisonChromosome;
        mCheckTranscript = checkTranscript;

        String key = format("%s_%s", driverCatalog.driver(), driverCatalog.gene());
        mKey = driverCatalog.isCanonical() ? key : key + "_" + driverCatalog.transcript();
        mIsPass = isPass;

        addStringValue(DriverComparer.Fields.LikelihoodMethod.toString(), driverCatalog.likelihoodMethod().toString(), fields);
        addDoubleValue(DriverComparer.Fields.Likelihood.toString(), driverCatalog.driverLikelihood(), fields);
        addDoubleValue(DriverComparer.Fields.MinCopyNumber.toString(), driverCatalog.minCopyNumber(), fields);
        addDoubleValue(DriverComparer.Fields.MaxCopyNumber.toString(), driverCatalog.maxCopyNumber(), fields);
        addStringValue(DriverComparer.Fields.Chromosome.toString(), comparisonChromosome, fields);
        addStringValue(DriverComparer.Fields.ChromosomeBand.toString(), driverCatalog.chromosomeBand(), fields);

        addDoubleValue(DriverComparer.Fields.Purity.toString(), purity.Purity, fields);
        addDoubleValue(DriverComparer.Fields.Ploidy.toString(), purity.Ploidy, fields);
    }

    @Override
    public CategoryType category() { return DRIVER; }

    @Override
    public String key()
    {
        return mKey;
    }

    @Override
    public boolean isPass()
    {
        return mIsPass;
    }

    @Override
    public boolean reportable() { return DriverCatalog.reportedStatus() == ReportedStatus.REPORTED; }

    @Override
    public String geneName() { return DriverCatalog.gene(); }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final DriverData otherDriver = (DriverData)other;

        if(!DriverCatalog.gene().equals(otherDriver.DriverCatalog.gene()))
            return false;

        if(DriverCatalog.driver() != otherDriver.DriverCatalog.driver())
            return false;

        if(mCheckTranscript || otherDriver.mCheckTranscript)
        {
            if(!DriverCatalog.transcript().equals(otherDriver.DriverCatalog.transcript())
            && DriverCatalog.isCanonical() != otherDriver.DriverCatalog.isCanonical())
            {
                return false;
            }
        }

        return true;
    }

    public Mismatch findMismatch(
            final ItemComparer comparer, final ComparableItem other, final MatchLevel matchLevel, final boolean includeMatches)
    {
        // comparison is custom here so that non-reportable drivers can skip comparisons if synthetically created
        List<String> diffs = Lists.newArrayList();

        // find and compare fields present in both items
        List<FieldInfo> fields = comparer.fieldsList();
        Map<String, FieldValue> oldFieldValues = fieldValues();
        Map<String,FieldValue> newFieldValues = other.fieldValues();

        boolean bothPass = isPass() && other.isPass();

        for(FieldInfo field : fields)
        {
            if(!field.FieldCheck.IsCompared)
                continue;

            if(!bothPass && field.Name.equals(DriverComparer.Fields.Likelihood.toString()))
                continue;

            FieldValue oldValue = oldFieldValues.get(field.Name);
            FieldValue newValue = newFieldValues.get(field.Name);

            if(oldValue == null || newValue == null)
                continue;

            if(oldValue.hasDifference(newValue))
            {
                oldValue.addDiffInfo(oldValue, newValue, diffs);
            }
        }


        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }
}

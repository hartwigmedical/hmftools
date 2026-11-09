package com.hartwig.hmftools.compar.purple;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.GENE_COPY_NUMBER;

import java.util.List;

import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.ExpectedValue;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class GeneCopyNumberData extends ComparableItem
{
    public final GeneCopyNumber CopyNumber;

    public GeneCopyNumberData(final GeneCopyNumber copyNumber, final List<FieldInfo> fields)
    {
        CopyNumber = copyNumber;

        addDoubleValue(GeneCopyNumberComparer.Fields.MinCopyNumber.toString(), copyNumber.MinCopyNumber, fields);
        addDoubleValue(GeneCopyNumberComparer.Fields.MaxCopyNumber.toString(), copyNumber.MaxCopyNumber, fields);
        addStringValue(GeneCopyNumberComparer.Fields.DriverType.toString(), copyNumber.driverType().toString(), fields);

        addIntValue(GeneCopyNumberComparer.Fields.MinRegionStart.toString(), copyNumber.MinRegionStart, fields);
        addIntValue(GeneCopyNumberComparer.Fields.MinRegionEnd.toString(), copyNumber.MinRegionEnd, fields);
    }

    private GeneCopyNumberData(final GeneCopyNumber copyNumber)
    {
        CopyNumber = copyNumber;
    }

    public static GeneCopyNumberData fromExpected(final List<ExpectedValue> expectedValues, final List<FieldInfo> fields)
    {
        String geneName = expectedValues.get(0).Key;
        double minCopyNumber = 0;
        double maxCopyNumber = 0;
        DriverType driverType = DriverType.UNKNOWN;

        for(ExpectedValue expectedValue : expectedValues)
        {
            if(expectedValue.FieldName.equals(GeneCopyNumberComparer.Fields.MinCopyNumber.toString()))
                minCopyNumber = Double.parseDouble(expectedValue.Value);
            else if(expectedValue.FieldName.equals(GeneCopyNumberComparer.Fields.MaxCopyNumber.toString()))
                maxCopyNumber = Double.parseDouble(expectedValue.Value);
            else if(expectedValue.FieldName.equals(GeneCopyNumberComparer.Fields.DriverType.toString()))
                driverType = DriverType.valueOf(expectedValue.Value);
        }

        GeneCopyNumber copyNumber = new GeneCopyNumber(
                "", 0, 0, geneName, "", true, "", maxCopyNumber, minCopyNumber, 0, 0, 0, 0, 0, 0, 0,
                SegmentSupport.NONE, SegmentSupport.NONE, CopyNumberMethod.UNKNOWN, 0);

        copyNumber.setDriverType(driverType);

        GeneCopyNumberData geneCopyNumber = new GeneCopyNumberData(copyNumber);
        geneCopyNumber.addExpectedValues(expectedValues, fields);
        return geneCopyNumber;
    }

    private void addExpectedValues(final List<ExpectedValue> expectedValues, final List<FieldInfo> fields)
    {
        for(ExpectedValue expectedValue : expectedValues)
        {
            if(expectedValue.FieldName.equals(GeneCopyNumberComparer.Fields.MinCopyNumber.toString()))
                addDoubleValue(GeneCopyNumberComparer.Fields.MinCopyNumber.toString(), CopyNumber.MinCopyNumber, fields);
            else if(expectedValue.FieldName.equals(GeneCopyNumberComparer.Fields.MaxCopyNumber.toString()))
                addDoubleValue(GeneCopyNumberComparer.Fields.MaxCopyNumber.toString(), CopyNumber.MaxCopyNumber, fields);
            else if(expectedValue.FieldName.equals(GeneCopyNumberComparer.Fields.DriverType.toString()))
                addStringValue(GeneCopyNumberComparer.Fields.DriverType.toString(), CopyNumber.driverType().toString(), fields);
        }
    }

    public CategoryType category() {
        return GENE_COPY_NUMBER;
    }

    @Override
    public String key()
    {
        return format("%s", CopyNumber.geneName());
    }

    @Override
    public boolean reportable() {
        return CopyNumber.reportedStatus() == ReportedStatus.REPORTED;
    }

    @Override
    public String geneName() { return CopyNumber.GeneName; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final GeneCopyNumberData otherCn = (GeneCopyNumberData)other;

        return CopyNumber.geneName().equals(otherCn.CopyNumber.geneName());
    }
}

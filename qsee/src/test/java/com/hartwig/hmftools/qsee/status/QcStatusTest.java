package com.hartwig.hmftools.qsee.status;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.qsee.feature.NumberFormat;

import org.junit.Test;

public class QcStatusTest
{
    @Test
    public void canFormatPercentThresholdsCorrectly()
    {
        String displayString;

        displayString = QcStatus.formDiplayString(QcStatusType.FAIL, ComparisonOperator.LESS_THAN, 0.10, NumberFormat.NUMBER);
        assertEquals("<0.1", displayString);

        displayString = QcStatus.formDiplayString(QcStatusType.FAIL, ComparisonOperator.LESS_THAN, 0.10, NumberFormat.PERCENT);
        assertEquals("<10%", displayString);

        displayString = QcStatus.formDiplayString(QcStatusType.FAIL, ComparisonOperator.LESS_THAN, 0.101, NumberFormat.PERCENT);
        assertEquals("<10.1%", displayString);
    }

    @Test
    public void canFormatNoThresholdsCorrectly()
    {
        String displayString;

        displayString = QcStatus.formDiplayString(QcStatusType.FAIL, ComparisonOperator.LESS_THAN, Double.NaN, NumberFormat.PERCENT);
        assertEquals("", displayString);

        displayString = QcStatus.formDiplayString(QcStatusType.FAIL, ComparisonOperator.LESS_THAN, Double.NaN, NumberFormat.NUMBER);
        assertEquals("", displayString);
    }
}

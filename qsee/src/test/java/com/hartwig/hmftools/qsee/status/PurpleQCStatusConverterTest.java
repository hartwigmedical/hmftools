package com.hartwig.hmftools.qsee.status;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.purple.PurpleQCStatus;

import org.junit.Test;

public class PurpleQCStatusConverterTest
{
    @Test
    public void canMapPurpleQCStatusToQseeQCStatus()
    {
        ThresholdRegistry qcThresholds = ThresholdRegistry.createWithoutThresholds();
        PurpleQCStatusConverter converter = new PurpleQCStatusConverter(qcThresholds);

        List<QcStatus> qcStatuses = new ArrayList<>();
        for(PurpleQCStatus purpleQCStatus : PurpleQCStatus.values())
        {
            QcStatus qcStatus = converter.toQseeQcStatus(purpleQCStatus);
            qcStatuses.add(qcStatus);
        }

        assertEquals(PurpleQCStatus.values().length, qcStatuses.size());
    }
}

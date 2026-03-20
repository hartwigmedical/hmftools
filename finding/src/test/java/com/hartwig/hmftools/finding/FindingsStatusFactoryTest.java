package com.hartwig.hmftools.finding;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.finding.datamodel.finding.FindingsStatus;
import com.hartwig.hmftools.finding.datamodel.finding.ResultIssue;
import com.hartwig.hmftools.finding.datamodel.finding.ResultStatus;

import org.junit.Test;

public class FindingsStatusFactoryTest
{
    @Test
    public void toFindingStatusOk() {
        FindingsStatus findingsStatus = FindingsStatusFactory.toFindingsStatus(Set.of(PurpleQCStatus.WARN_LOW_PURITY));
        assertEquals(ResultStatus.OK, findingsStatus.status());
        assertTrue(findingsStatus.isOK());
        assertEquals(findingsStatus.warnings(), Set.of(ResultIssue.LOW_PURITY));
        assertEquals(findingsStatus.errors(), Set.of());
    }

    @Test
    public void toFindingStatusNotReliable() {
        FindingsStatus findingsStatus = FindingsStatusFactory.toFindingsStatus(Set.of(PurpleQCStatus.FAIL_CONTAMINATION, PurpleQCStatus.WARN_LOW_PURITY));
        assertEquals(ResultStatus.NOT_RELIABLE, findingsStatus.status());
        assertFalse(findingsStatus.isOK());
        assertEquals(findingsStatus.warnings(), Set.of(ResultIssue.LOW_PURITY));
        assertEquals(findingsStatus.errors(), Set.of(ResultIssue.CONTAMINATION));
    }
}

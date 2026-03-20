package com.hartwig.hmftools.finding;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.finding.ResultIssue;
import com.hartwig.hmftools.finding.datamodel.finding.ResultStatus;

import org.junit.Test;

public class FindingStatusFactoryTest
{
    @Test
    public void toFindingStatusOk() {
        FindingStatus findingStatus = FindingsStatusFactory.toFindingsStatus(Set.of(PurpleQCStatus.WARN_LOW_PURITY));
        assertEquals(ResultStatus.OK, findingStatus.status());
        assertTrue(findingStatus.isOK());
        assertEquals(findingStatus.warnings(), Set.of(ResultIssue.LOW_PURITY));
        assertEquals(findingStatus.errors(), Set.of());
    }

    @Test
    public void toFindingStatusNotReliable() {
        FindingStatus findingStatus = FindingsStatusFactory.toFindingsStatus(Set.of(PurpleQCStatus.FAIL_CONTAMINATION, PurpleQCStatus.WARN_LOW_PURITY));
        assertEquals(ResultStatus.NOT_RELIABLE, findingStatus.status());
        assertFalse(findingStatus.isOK());
        assertEquals(findingStatus.warnings(), Set.of(ResultIssue.LOW_PURITY));
        assertEquals(findingStatus.errors(), Set.of(ResultIssue.CONTAMINATION));
    }
}

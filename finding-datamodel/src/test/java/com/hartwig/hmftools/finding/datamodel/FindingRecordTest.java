package com.hartwig.hmftools.finding.datamodel;

import static org.junit.Assert.assertEquals;

import java.util.Set;
import java.util.TreeSet;

import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatusBuilder;

import org.junit.Test;

public class FindingRecordTest
{
    @Test
    public void testMergeFindingStatus()
    {
        testFindingStatus(createFindingStatus(FindingStatus.Status.OK, Set.of(), Set.of(FindingStatus.Issue.GENDER_MISMATCH)),
                createFindingStatus(FindingStatus.Status.OK, Set.of(), Set.of()),
                createFindingStatus(FindingStatus.Status.OK, Set.of(), Set.of(FindingStatus.Issue.GENDER_MISMATCH)));
        testFindingStatus(createFindingStatus(FindingStatus.Status.OK, Set.of(), Set.of(FindingStatus.Issue.NO_TUMOR)),
                createFindingStatus(FindingStatus.Status.OK, Set.of(), Set.of()),
                createFindingStatus(FindingStatus.Status.NOT_RELIABLE, Set.of(FindingStatus.Issue.NO_TUMOR), Set.of()));
        testFindingStatus(createFindingStatus(FindingStatus.Status.OK, Set.of(), Set.of(FindingStatus.Issue.REF_REQUIRED)),
                createFindingStatus(FindingStatus.Status.OK, Set.of(), Set.of()),
                createFindingStatus(FindingStatus.Status.NOT_AVAILABLE, Set.of(FindingStatus.Issue.REF_REQUIRED), Set.of()));
        testFindingStatus(createFindingStatus(FindingStatus.Status.NOT_RELIABLE, Set.of(FindingStatus.Issue.NO_TUMOR, FindingStatus.Issue.LOW_PURITY), Set.of()),
                createFindingStatus(FindingStatus.Status.NOT_RELIABLE, Set.of(FindingStatus.Issue.NO_TUMOR), Set.of()),
                createFindingStatus(FindingStatus.Status.NOT_RELIABLE, Set.of(FindingStatus.Issue.LOW_PURITY), Set.of()));
        testFindingStatus(createFindingStatus(FindingStatus.Status.NOT_AVAILABLE, Set.of(FindingStatus.Issue.WGS_REQUIRED, FindingStatus.Issue.REF_REQUIRED), Set.of()),
                createFindingStatus(FindingStatus.Status.NOT_AVAILABLE, Set.of(FindingStatus.Issue.WGS_REQUIRED), Set.of()),
                createFindingStatus(FindingStatus.Status.NOT_AVAILABLE, Set.of(FindingStatus.Issue.REF_REQUIRED), Set.of()));
        testFindingStatus(createFindingStatus(FindingStatus.Status.NOT_RELIABLE, Set.of(FindingStatus.Issue.NO_TUMOR, FindingStatus.Issue.REF_REQUIRED), Set.of()),
                createFindingStatus(FindingStatus.Status.NOT_RELIABLE, Set.of(FindingStatus.Issue.NO_TUMOR), Set.of()),
                createFindingStatus(FindingStatus.Status.NOT_AVAILABLE, Set.of(FindingStatus.Issue.REF_REQUIRED), Set.of()));
    }

    private void testFindingStatus(FindingStatus expected, FindingStatus status1, FindingStatus status2)
    {
        assertFindingStatus(expected, FindingRecord.mergeFindingStatus(status1, status2));
        assertFindingStatus(expected, FindingRecord.mergeFindingStatus(status2, status1));
    }

    private void assertFindingStatus(FindingStatus expected, FindingStatus actual)
    {
        assertEquals(expected.status(), actual.status());
        assertEquals(expected.errors(), actual.errors());
        assertEquals(expected.warnings(), actual.warnings());
    }

    private static FindingStatus createFindingStatus(FindingStatus.Status status, Set<FindingStatus.Issue> errors,
            Set<FindingStatus.Issue> warnings)
    {
        return FindingStatusBuilder.builder().status(status).errors(new TreeSet<>(errors)).warnings(new TreeSet<>(warnings)).build();
    }
}

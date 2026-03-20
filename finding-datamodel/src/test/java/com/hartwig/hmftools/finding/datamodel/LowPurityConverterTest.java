package com.hartwig.hmftools.finding.datamodel;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.finding.FindingItem;
import com.hartwig.hmftools.finding.datamodel.finding.FindingList;
import com.hartwig.hmftools.finding.datamodel.finding.FindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatusBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.ResultIssue;
import com.hartwig.hmftools.finding.datamodel.finding.ResultStatus;
import com.hartwig.hmftools.finding.util.LowPurityConverter;

import org.junit.Test;

public class LowPurityConverterTest
{
    @Test
    public void testConvert()
    {
        // TODO: Create this from actual orange record to make sure it's a representative finding record
        FindingRecord original = TestFindingRecordFactory.createMinimalTestFindingRecordBuilder()
                .qc(TestFindingFactory.qcBuilder()
                        .status(Set.of(Qc.QCStatus.PASS, Qc.QCStatus.WARN_LOW_PURITY))
                        .build())
                .hlaAlleles(FindingListBuilder.<HlaAllele>builder()
                        .status(FindingStatusBuilder.builder()
                                .status(ResultStatus.OK)
                                .warnings(new TreeSet<>(Set.of(ResultIssue.LOW_PURITY)))
                                .errors(new TreeSet<>())
                                .build())
                        .findings(List.of(TestFindingFactory.hlaAlleleBuilder().build()))
                        .build())
                .build();
        FindingRecord converted = LowPurityConverter.convert(original);

        assertFindingList(converted.somaticDisruptions());
        assertFindingList(converted.somaticGainDeletions());
        assertFindingList(converted.viruses());
        assertHLA(converted.hlaAlleles());

        assertFindingItem(converted.homologousRecombination());
        assertFindingItem(converted.microsatelliteStability());
        assertFindingItem(converted.tumorMutationalBurden());
        assertFindingItem(converted.tumorMutationalLoad());
    }

    private void assertFindingItem(FindingItem<?> item)
    {
        assertFindingStatus(item.status());
        assertNull(item.finding());
    }

    private void assertFindingList(DriverFindingList<?> list)
    {
        assertFindingStatus(list.status());
        assertTrue(list.findings().isEmpty());
    }

    private void assertFindingStatus(FindingStatus findingStatus) {
        assertEquals(ResultStatus.NOT_RELIABLE, findingStatus.status());
        assertTrue(findingStatus.errors().contains(ResultIssue.LOW_PURITY));
        assertFalse(findingStatus.warnings().contains(ResultIssue.LOW_PURITY));
    }

    private void assertHLA(FindingList<HlaAllele> findingList) {
        FindingStatus findingStatus = findingList.status();
        assertEquals(ResultStatus.OK, findingStatus.status());
        assertFalse(findingStatus.errors().contains(ResultIssue.LOW_PURITY));
        assertTrue(findingStatus.warnings().contains(ResultIssue.LOW_PURITY));
        List<HlaAllele> hlaAlleles = findingList.findings();
        assertFalse(hlaAlleles.isEmpty());
        for(HlaAllele hlaAllele : hlaAlleles)
        {
            assertNull(hlaAllele.tumorCopyNumber());
        }
    }
}

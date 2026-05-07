package com.hartwig.hmftools.finding.util;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.HlaAllele;
import com.hartwig.hmftools.finding.datamodel.TestFindingFactory;
import com.hartwig.hmftools.finding.datamodel.TestFindingRecordFactory;
import com.hartwig.hmftools.finding.datamodel.finding.FindingList;
import com.hartwig.hmftools.finding.datamodel.finding.FindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatusBuilder;

import org.junit.Test;

public class ErrorConverterTest
{
    @Test
    public void testConvertHLANoTumor()
    {
        // TODO: Create this from actual orange record to make sure it's a representative finding record
        FindingRecord original = TestFindingRecordFactory.createMinimalTestFindingRecordBuilder()
                .hlaAlleles(FindingListBuilder.<HlaAllele>builder()
                        .status(FindingStatusBuilder.builder()
                                .status(FindingStatus.Status.NOT_RELIABLE)
                                .errors(new TreeSet<>(Set.of(FindingStatus.Issue.NO_TUMOR)))
                                .warnings(new TreeSet<>())
                                .build())
                        .findings(List.of(TestFindingFactory.hlaAlleleBuilder().build()))
                        .build())
                .build();
        FindingRecord converted = ErrorConverter.convert(original);

        FindingList<HlaAllele> findingList = converted.hlaAlleles();
        FindingStatus findingStatus = findingList.status();
        assertEquals(FindingStatus.Status.OK, findingStatus.status());
        assertFalse(findingStatus.errors().contains(FindingStatus.Issue.NO_TUMOR));
        assertTrue(findingStatus.warnings().contains(FindingStatus.Issue.NO_TUMOR));
        List<HlaAllele> hlaAlleles = findingList.findings();
        assertFalse(hlaAlleles.isEmpty());
        for(HlaAllele hlaAllele : hlaAlleles)
        {
            assertNull(hlaAllele.tumorCopyNumber());
        }
    }

    @Test
    public void testConvertHLANoPass()
    {
        // TODO: Create this from actual orange record to make sure it's a representative finding record
        FindingRecord original = TestFindingRecordFactory.createMinimalTestFindingRecordBuilder()
                .hlaAlleles(FindingListBuilder.<HlaAllele>builder()
                        .status(FindingStatusBuilder.builder()
                                .status(FindingStatus.Status.OK)
                                .errors(new TreeSet<>())
                                .warnings(new TreeSet<>())
                                .build())
                        .findings(List.of(TestFindingFactory.hlaAlleleBuilder().qcStatus(Set.of(HlaAllele.QcStatus.FAIL)).build()))
                        .build())
                .build();
        FindingRecord converted = ErrorConverter.convert(original);

        FindingList<HlaAllele> findingList = converted.hlaAlleles();
        FindingStatus findingStatus = findingList.status();
        assertEquals(FindingStatus.Status.NOT_AVAILABLE, findingStatus.status());
        assertTrue(findingStatus.errors().contains(FindingStatus.Issue.NO_REPORTABLE_VALUE));
        assertTrue(findingStatus.warnings().isEmpty());
        List<HlaAllele> hlaAlleles = findingList.findings();
        assertTrue(hlaAlleles.isEmpty());
    }
}

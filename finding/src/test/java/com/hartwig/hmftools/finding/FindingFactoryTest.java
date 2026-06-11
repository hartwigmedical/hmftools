package com.hartwig.hmftools.finding;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.datamodel.TestOrangeJsonWriter;
import com.hartwig.hmftools.datamodel.isofox.AltSpliceJunctionContext;
import com.hartwig.hmftools.datamodel.isofox.AltSpliceJunctionType;
import com.hartwig.hmftools.datamodel.isofox.ImmutableGeneExpression;
import com.hartwig.hmftools.datamodel.isofox.ImmutableIsofoxRecord;
import com.hartwig.hmftools.datamodel.isofox.ImmutableNovelSpliceJunction;
import com.hartwig.hmftools.datamodel.isofox.ImmutableRnaFusion;
import com.hartwig.hmftools.datamodel.isofox.ImmutableRnaStatistics;
import com.hartwig.hmftools.datamodel.isofox.RnaFusionType;
import com.hartwig.hmftools.datamodel.isofox.RnaQCStatus;
import com.hartwig.hmftools.datamodel.isofox.StructuralVariantType;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeRecord;
import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.NovelSpliceJunction;
import com.hartwig.hmftools.finding.datamodel.RnaFusion;
import com.hartwig.hmftools.finding.datamodel.RnaQc;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;

import org.junit.Ignore;
import org.junit.Test;

public class FindingFactoryTest
{
    @Ignore
    @Test
    public void fromOrangeRecord() throws IOException
    {
        FindingRecordFactory.fromOrangeRecord(TestOrangeJsonWriter.createOrangeRecord(), null);
    }

    @Test
    public void populatesRnaFindingsFromIsofox() throws IOException
    {
        FindingRecord findingRecord = FindingRecordFactory.fromOrangeRecord(ImmutableOrangeRecord.builder()
                .from(TestOrangeJsonWriter.createOrangeRecord())
                .referenceId(null)
                .isofox(ImmutableIsofoxRecord.builder()
                        .summary(ImmutableRnaStatistics.builder()
                                .addQcStatus(RnaQCStatus.PASS)
                                .totalFragments(100)
                                .duplicateFragments(10)
                                .splicedFragmentPerc(0.1)
                                .unsplicedFragmentPerc(0.2)
                                .altFragmentPerc(0.3)
                                .chimericFragmentPerc(0.4)
                                .build())
                        .highExpressionGenes(List.of(ImmutableGeneExpression.builder()
                                .gene("ERBB2")
                                .tpm(20)
                                .medianTpmCohort(5)
                                .percentileCohort(0.95)
                                .medianTpmCancer(6D)
                                .percentileCancer(0.9)
                                .build()))
                        .lowExpressionGenes(List.of(ImmutableGeneExpression.builder()
                                .gene("PTEN")
                                .tpm(1)
                                .medianTpmCohort(8)
                                .percentileCohort(0.05)
                                .build()))
                        .fusions(List.of(ImmutableRnaFusion.builder()
                                .geneStart("EML4")
                                .geneEnd("ALK")
                                .chromosomeStart("2")
                                .chromosomeEnd("2")
                                .positionStart(42522656)
                                .positionEnd(29446394)
                                .junctionTypeStart("EXON")
                                .junctionTypeEnd("EXON")
                                .knownType(RnaFusionType.KNOWN_PAIR)
                                .svType(StructuralVariantType.INV)
                                .splitFragments(12)
                                .realignedFrags(3)
                                .discordantFrags(4)
                                .depthStart(100)
                                .depthEnd(90)
                                .cohortFrequency(2)
                                .build()))
                        .novelSpliceJunctions(List.of(ImmutableNovelSpliceJunction.builder()
                                .gene("MET")
                                .chromosome("7")
                                .junctionStart(116411708)
                                .junctionEnd(116412043)
                                .type(AltSpliceJunctionType.SKIPPED_EXONS)
                                .exonStart(13)
                                .exonEnd(15)
                                .fragmentCount(8)
                                .depthStart(80)
                                .depthEnd(75)
                                .regionStart(AltSpliceJunctionContext.SPLICE_JUNC)
                                .regionEnd(AltSpliceJunctionContext.SPLICE_JUNC)
                                .cohortFrequency(1)
                                .build()))
                        .build())
                .build(), null);

        RnaQc rnaQc = findingRecord.rnaQc();
        assertNotNull(rnaQc);
        assertTrue(rnaQc.errors().isEmpty());
        assertTrue(rnaQc.warnings().isEmpty());
        assertEquals(100, rnaQc.totalFragments());
        assertEquals(0.4, rnaQc.chimericFragmentPercent(), 0);

        assertEquals("ERBB2", findingRecord.highExpressionGenes().findings().get(0).gene());
        assertEquals("PTEN", findingRecord.lowExpressionGenes().findings().get(0).gene());

        RnaFusion fusion = findingRecord.rnaFusions().findings().get(0);
        assertEquals("EML4::ALK", fusion.display());
        assertEquals(RnaFusion.KnownType.KNOWN_PAIR, fusion.knownType());
        assertEquals(RnaFusion.StructuralVariantType.INV, fusion.structuralVariantType());
        assertEquals(3, fusion.realignedFragments());

        NovelSpliceJunction spliceJunction = findingRecord.novelSpliceJunctions().findings().get(0);
        assertEquals("MET", spliceJunction.gene());
        assertEquals(NovelSpliceJunction.Type.SKIPPED_EXONS, spliceJunction.type());
        assertEquals(NovelSpliceJunction.Context.SPLICE_JUNC, spliceJunction.regionStart());
    }

    @Test
    public void setsRnaFindingsToNotReliableWhenIsofoxQcFails() throws IOException
    {
        FindingRecord findingRecord = FindingRecordFactory.fromOrangeRecord(ImmutableOrangeRecord.builder()
                .from(TestOrangeJsonWriter.createOrangeRecord())
                .referenceId(null)
                .isofox(ImmutableIsofoxRecord.builder()
                        .summary(ImmutableRnaStatistics.builder()
                                .addQcStatus(RnaQCStatus.FAIL_LOW_COVERAGE)
                                .totalFragments(100)
                                .duplicateFragments(10)
                                .splicedFragmentPerc(0.1)
                                .unsplicedFragmentPerc(0.2)
                                .altFragmentPerc(0.3)
                                .chimericFragmentPerc(0.4)
                                .build())
                        .highExpressionGenes(List.of())
                        .lowExpressionGenes(List.of())
                        .fusions(List.of())
                        .novelSpliceJunctions(List.of())
                        .build())
                .build(), null);

        assertNotNull(findingRecord.rnaQc());
        assertEquals(Set.of(RnaQc.QcStatus.FAIL_LOW_COVERAGE), findingRecord.rnaQc().errors());
        assertRnaSampleQcNotReliable(findingRecord.highExpressionGenes().status());
        assertRnaSampleQcNotReliable(findingRecord.lowExpressionGenes().status());
        assertRnaSampleQcNotReliable(findingRecord.rnaFusions().status());
        assertRnaSampleQcNotReliable(findingRecord.novelSpliceJunctions().status());
    }

    private static void assertRnaSampleQcNotReliable(FindingStatus status)
    {
        assertEquals(FindingStatus.Status.NOT_RELIABLE, status.status());
        assertEquals(Set.of(FindingStatus.Issue.RNA_SAMPLE_QUALITY_CONTROL), status.errors());
    }
}

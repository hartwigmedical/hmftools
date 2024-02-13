package com.hartwig.hmftools.wisp;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.wisp.TestUtils.MOCK_REF_GENOME;
import static com.hartwig.hmftools.wisp.TestUtils.REF_BASES_CHR_1;
import static com.hartwig.hmftools.wisp.TestUtils.REF_BASES_CHR_2;
import static com.hartwig.hmftools.wisp.TestUtils.TEST_CONFIG;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.gene.TranscriptCodingType;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;
import com.hartwig.hmftools.common.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.sv.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.wisp.probe.StructuralVariant;

import org.junit.Test;

public class StructuralVariantTest
{
    @Test
    public void testStructuralVariantSequences()
    {
        MOCK_REF_GENOME.RefGenomeMap.put(CHR_1, REF_BASES_CHR_1);
        MOCK_REF_GENOME.RefGenomeMap.put(CHR_2, REF_BASES_CHR_2);

        // a DEL
        StructuralVariant var = new StructuralVariant(
                createSv(0, "001", CHR_1, CHR_1, 20, 30, POS_ORIENT, NEG_ORIENT, StructuralVariantType.DEL,
                        2, 2, ""),
                Lists.newArrayList(), Lists.newArrayList());
        var.markAmpDelDriver(false);

        var.generateSequences(MOCK_REF_GENOME, TEST_CONFIG);
        String sequence = REF_BASES_CHR_1.substring(11, 21) + REF_BASES_CHR_1.substring(30, 40);
        assertEquals(sequence, var.sequence());
        assertEquals(REF_BASES_CHR_1.substring(10, 10 + TEST_CONFIG.ProbeLength), var.refSequences().get(0));
        assertEquals(REF_BASES_CHR_1.substring(20, 20 + TEST_CONFIG.ProbeLength), var.refSequences().get(1));

        // DUP with an insert
        String insertSeq = "AAAAA";
        var = new StructuralVariant(
                createSv(0, "001", CHR_1, CHR_1, 20, 40, NEG_ORIENT, POS_ORIENT, StructuralVariantType.DUP,
                        2, 2, insertSeq),
                Lists.newArrayList(), Lists.newArrayList());

        var.generateSequences(MOCK_REF_GENOME, TEST_CONFIG);
        sequence = REF_BASES_CHR_1.substring(33, 41) + insertSeq + REF_BASES_CHR_1.substring(20, 27);
        assertEquals(sequence, var.sequence());

        // BND with -1/-1
        insertSeq = "AAAAAA";
        var = new StructuralVariant(
                createSv(0, "001", CHR_1, CHR_2, 20, 40, NEG_ORIENT, NEG_ORIENT, StructuralVariantType.BND,
                        2, 2, insertSeq),
                Lists.newArrayList(), Lists.newArrayList());

        var.generateSequences(MOCK_REF_GENOME, TEST_CONFIG);
        sequence = Nucleotides.reverseComplementBases(REF_BASES_CHR_2.substring(40, 47)) + insertSeq + REF_BASES_CHR_1.substring(20, 27);
        assertEquals(sequence, var.sequence());

        // SGL with positive orientation
        insertSeq = "AAAAA";
        LinxBreakend reportableBreakend = createBreakend(true);
        var = new StructuralVariant(
                createSv(0, "001", CHR_1, "-1", 20, 0, POS_ORIENT, 0, StructuralVariantType.SGL,
                        2, 0, insertSeq),
                Lists.newArrayList(reportableBreakend), Lists.newArrayList());

        var.generateSequences(MOCK_REF_GENOME, TEST_CONFIG);
        sequence = REF_BASES_CHR_1.substring(6, 21) + insertSeq;
        assertEquals(sequence, var.sequence());
        assertEquals(REF_BASES_CHR_1.substring(10, 10 + TEST_CONFIG.ProbeLength), var.refSequences().get(0));

        // with a longer insert sequence - will only allocate half the probe length
        insertSeq = "AAAAAGGGGGCCCCCTTTTT";
        var = new StructuralVariant(
                createSv(0, "001", CHR_1, "-1", 20, 0, POS_ORIENT, 0, StructuralVariantType.SGL,
                        2, 0, insertSeq),
                Lists.newArrayList(reportableBreakend), Lists.newArrayList());

        var.generateSequences(MOCK_REF_GENOME, TEST_CONFIG);
        sequence = REF_BASES_CHR_1.substring(11, 21) + insertSeq.substring(0, 10);
        assertEquals(sequence, var.sequence());
        assertEquals(REF_BASES_CHR_1.substring(10, 10 + TEST_CONFIG.ProbeLength), var.refSequences().get(0));

        // negative orientation
        var = new StructuralVariant(
                createSv(0, "001", CHR_1, "-1", 20, 0, NEG_ORIENT, 0, StructuralVariantType.SGL,
                        2, 0, insertSeq),
                Lists.newArrayList(reportableBreakend), Lists.newArrayList());

        var.generateSequences(MOCK_REF_GENOME, TEST_CONFIG);
        sequence = insertSeq.substring(insertSeq.length() - 10) + REF_BASES_CHR_1.substring(20, 30);
        assertEquals(sequence, var.sequence());
        assertEquals(REF_BASES_CHR_1.substring(10, 10 + TEST_CONFIG.ProbeLength), var.refSequences().get(0));
    }

    public static StructuralVariantData createSv(
            final int varId, final String vcfId, final String chrStart, final String chrEnd,
            int posStart, int posEnd, int orientStart, int orientEnd, StructuralVariantType type,
            double cnChgStart, double cnChgEnd, final String insertSeq)
    {
        return ImmutableStructuralVariantData.builder()
                        .id(varId)
                        .vcfId(vcfId)
                        .startChromosome(chrStart)
                        .endChromosome(chrEnd)
                        .startPosition(posStart)
                        .endPosition(posEnd)
                        .startOrientation((byte) orientStart)
                        .endOrientation((byte) orientEnd)
                        .startHomologySequence("")
                        .endHomologySequence("")
                        .startAF(1.0)
                        .endAF(1.0)
                        .junctionCopyNumber(1)
                        .adjustedStartAF(1.0)
                        .adjustedEndAF(1.0)
                        .adjustedStartCopyNumber(2)
                        .adjustedEndCopyNumber(2)
                        .adjustedStartCopyNumberChange(cnChgStart)
                        .adjustedEndCopyNumberChange(cnChgEnd)
                        .insertSequence(insertSeq)
                        .type(type)
                        .filter(PASS)
                        .imprecise(false)
                        .qualityScore(0.0)
                        .event("")
                        .startTumorVariantFragmentCount(10)
                        .startTumorReferenceFragmentCount(10)
                        .startNormalVariantFragmentCount(10)
                        .startNormalReferenceFragmentCount(10)
                        .endTumorVariantFragmentCount(10)
                        .endTumorReferenceFragmentCount(10)
                        .endNormalVariantFragmentCount(10)
                        .endNormalReferenceFragmentCount(10)
                        .startIntervalOffsetStart(0)
                        .startIntervalOffsetEnd(0)
                        .endIntervalOffsetStart(0)
                        .endIntervalOffsetEnd(0)
                        .inexactHomologyOffsetStart(0)
                        .inexactHomologyOffsetEnd(0)
                        .startLinkedBy("")
                        .endLinkedBy("")
                        .startRefContext("")
                        .endRefContext("")
                        .recovered(false)
                        .recoveryMethod("")
                        .recoveryFilter("")
                        .insertSequenceAlignments("")
                        .insertSequenceRepeatClass("")
                        .insertSequenceRepeatType("")
                        .insertSequenceRepeatOrientation((byte) 0)
                        .insertSequenceRepeatCoverage(0.0)
                        .startAnchoringSupportDistance(0)
                        .endAnchoringSupportDistance(0)
                        .ponCount(0)
                        .build();
    }

    private static LinxBreakend createBreakend(boolean reportable)
    {
        return ImmutableLinxBreakend.builder()
                .id(0)
                .svId(0)
                .isStart(true)
                .gene("GENE")
                .transcriptId("")
                .canonical(true)
                .geneOrientation("")
                .disruptive(true)
                .reportedDisruption(reportable)
                .undisruptedCopyNumber(1.0)
                .regionType(TranscriptRegionType.INTRONIC)
                .codingType(TranscriptCodingType.CODING)
                .biotype("")
                .exonicBasePhase(0)
                .nextSpliceExonRank(1)
                .nextSpliceExonPhase(1)
                .nextSpliceDistance(1000)
                .totalExonCount(5)
                .type(StructuralVariantType.DEL)
                .chromosome("1")
                .orientation(POS_ORIENT)
                .strand(1)
                .chrBand("")
                .exonUp(1)
                .exonDown(2)
                .junctionCopyNumber(1)
                .build();

    }

}

package com.hartwig.hmftools.panelbuilder.samplevariants;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;

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

import org.junit.Test;

// TODO

public class StructuralVariantTest
{
    @Test
    public void testStructuralVariantSequences()
    {
        TestUtils.MOCK_REF_GENOME.RefGenomeMap.put(CHR_1, TestUtils.REF_BASES_CHR_1);
        TestUtils.MOCK_REF_GENOME.RefGenomeMap.put(CHR_2, TestUtils.REF_BASES_CHR_2);

        // a DEL
        StructuralVariant var = new StructuralVariant(
                createSv(0, "001", CHR_1, CHR_1, 20, 30, ORIENT_FWD, ORIENT_REV, StructuralVariantType.DEL,
                        2, 2, ""),
                Lists.newArrayList(), Lists.newArrayList());
        var.markAmpDelDriver(false);

        var.generateSequences(TestUtils.MOCK_REF_GENOME);
        String sequence = TestUtils.REF_BASES_CHR_1.substring(11, 21) + TestUtils.REF_BASES_CHR_1.substring(30, 40);
        assertEquals(sequence, var.sequence());
        assertEquals(TestUtils.REF_BASES_CHR_1.substring(10, 10 + PROBE_LENGTH), var.refSequences().get(0));
        assertEquals(TestUtils.REF_BASES_CHR_1.substring(20, 20 + PROBE_LENGTH), var.refSequences().get(1));

        // DUP with an insert
        String insertSeq = "AAAAA";
        var = new StructuralVariant(
                createSv(0, "001", CHR_1, CHR_1, 20, 40, ORIENT_REV, ORIENT_FWD, StructuralVariantType.DUP,
                        2, 2, insertSeq),
                Lists.newArrayList(), Lists.newArrayList());

        var.generateSequences(TestUtils.MOCK_REF_GENOME);
        sequence = TestUtils.REF_BASES_CHR_1.substring(33, 41) + insertSeq + TestUtils.REF_BASES_CHR_1.substring(20, 27);
        assertEquals(sequence, var.sequence());

        // BND with -1/-1
        insertSeq = "AAAAAA";
        var = new StructuralVariant(
                createSv(0, "001", CHR_1, CHR_2, 20, 40, ORIENT_REV, ORIENT_REV, StructuralVariantType.BND,
                        2, 2, insertSeq),
                Lists.newArrayList(), Lists.newArrayList());

        var.generateSequences(TestUtils.MOCK_REF_GENOME);
        sequence = Nucleotides.reverseComplementBases(TestUtils.REF_BASES_CHR_2.substring(40, 47)) + insertSeq
                + TestUtils.REF_BASES_CHR_1.substring(20, 27);
        assertEquals(sequence, var.sequence());

        // SGL with positive orientation
        insertSeq = "AAAAA";
        LinxBreakend reportableBreakend = createBreakend(true);
        var = new StructuralVariant(
                createSv(0, "001", CHR_1, "-1", 20, 0, ORIENT_FWD, 0, StructuralVariantType.SGL,
                        2, 0, insertSeq),
                Lists.newArrayList(reportableBreakend), Lists.newArrayList());

        var.generateSequences(TestUtils.MOCK_REF_GENOME);
        sequence = TestUtils.REF_BASES_CHR_1.substring(6, 21) + insertSeq;
        assertEquals(sequence, var.sequence());
        assertEquals(TestUtils.REF_BASES_CHR_1.substring(10, 10 + PROBE_LENGTH), var.refSequences().get(0));

        // with a longer insert sequence - will only allocate half the probe length
        insertSeq = "AAAAAGGGGGCCCCCTTTTT";
        var = new StructuralVariant(
                createSv(0, "001", CHR_1, "-1", 20, 0, ORIENT_FWD, 0, StructuralVariantType.SGL,
                        2, 0, insertSeq),
                Lists.newArrayList(reportableBreakend), Lists.newArrayList());

        var.generateSequences(TestUtils.MOCK_REF_GENOME);
        sequence = TestUtils.REF_BASES_CHR_1.substring(11, 21) + insertSeq.substring(0, 10);
        assertEquals(sequence, var.sequence());
        assertEquals(TestUtils.REF_BASES_CHR_1.substring(10, 10 + PROBE_LENGTH), var.refSequences().get(0));

        // negative orientation
        var = new StructuralVariant(
                createSv(0, "001", CHR_1, "-1", 20, 0, ORIENT_REV, 0, StructuralVariantType.SGL,
                        2, 0, insertSeq),
                Lists.newArrayList(reportableBreakend), Lists.newArrayList());

        var.generateSequences(TestUtils.MOCK_REF_GENOME);
        sequence = insertSeq.substring(insertSeq.length() - 10) + TestUtils.REF_BASES_CHR_1.substring(20, 30);
        assertEquals(sequence, var.sequence());
        assertEquals(TestUtils.REF_BASES_CHR_1.substring(10, 10 + PROBE_LENGTH), var.refSequences().get(0));
    }

    public static StructuralVariantData createSv(
            final int varId, final String vcfId, final String chrStart, final String chrEnd,
            int posStart, int posEnd, int orientStart, int orientEnd, StructuralVariantType type,
            double cnChgStart, double cnChgEnd, final String insertSeq)
    {
        return ImmutableStructuralVariantData.builder()
                .id(varId)
                .vcfIdStart(vcfId)
                .vcfIdEnd(vcfId)
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
                .exonUp(1)
                .exonDown(2)
                .build();

    }

}

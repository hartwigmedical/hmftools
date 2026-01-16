package com.hartwig.hmftools.esvee.caller;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_INFO;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SPLIT_FRAGS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.STRAND_BIAS;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.setIlluminaSequencing;
import static com.hartwig.hmftools.esvee.caller.CallerTestUtils.createSv;
import static com.hartwig.hmftools.esvee.caller.FiltersTest.FILTER_CONSTANTS;
import static com.hartwig.hmftools.esvee.caller.FiltersTest.FRAG_LENGTHS;
import static com.hartwig.hmftools.esvee.caller.SeqTechUtils.SBX_HEURISTIC_SHORT_LENGTH;
import static com.hartwig.hmftools.esvee.caller.SeqTechUtils.SBX_INV_INSERT_MOTIFS;
import static com.hartwig.hmftools.esvee.caller.SeqTechUtils.SBX_INV_INSERT_MOTIF_MAX_DIFFS;
import static com.hartwig.hmftools.esvee.caller.SeqTechUtils.sequenceMismatches;
import static com.hartwig.hmftools.esvee.common.SvConstants.SEQUENCING_TYPE;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.esvee.common.FilterType;
import com.hartwig.hmftools.esvee.prep.types.DiscordantStats;

import org.junit.After;
import org.junit.Ignore;
import org.junit.Test;

public class SeqTechTest
{
    private final VariantFilters mVariantFilters = new VariantFilters(FILTER_CONSTANTS, FRAG_LENGTHS, new DiscordantStats());

    public SeqTechTest()
    {
        SEQUENCING_TYPE = SequencingType.SBX;
    }

    @After
    public void resetSequencingType() { setIlluminaSequencing(); }

    @Test
    public void testSbxStrandBias()
    {
        Map<String, Object> commonAttributes = Maps.newHashMap();

        Map<String, Object> tumorAttributes = Maps.newHashMap();
        tumorAttributes.put(SPLIT_FRAGS, 4);
        tumorAttributes.put(STRAND_BIAS, 1);

        // too few SF for an INV
        Variant var = createSv(
                "01", CHR_1, CHR_1, 100, 150, ORIENT_FWD, ORIENT_FWD, "",
                commonAttributes, tumorAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertFalse(var.filters().contains(FilterType.SBX_STRAND_BIAS));

        // then sufficient
        tumorAttributes.put(SPLIT_FRAGS, 12);

        var = createSv(
                "01", CHR_1, CHR_1, 100, 150, ORIENT_FWD, ORIENT_FWD, "",
                commonAttributes, tumorAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertTrue(var.filters().contains(FilterType.SBX_STRAND_BIAS));

        // any SF for a BND as long as strand bias is 0 or 1
        tumorAttributes.put(SPLIT_FRAGS, 1);
        tumorAttributes.put(STRAND_BIAS, 0.95);

        var = createSv(
                "01", CHR_1, CHR_2, 100, 150, ORIENT_FWD, ORIENT_FWD, "",
                commonAttributes, tumorAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertFalse(var.filters().contains(FilterType.SBX_STRAND_BIAS));

        tumorAttributes.put(SPLIT_FRAGS, 1);
        tumorAttributes.put(STRAND_BIAS, 0);

        var = createSv(
                "01", CHR_1, CHR_2, 100, 150, ORIENT_FWD, ORIENT_FWD, "",
                commonAttributes, tumorAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertTrue(var.filters().contains(FilterType.SBX_STRAND_BIAS));
    }

    @Test
    public void testSbxArtefacts()
    {
        // first test indels

        // test 1: non-line and stranded
        Map<String, Object> commonAttributes = Maps.newHashMap();

        commonAttributes.put(ASM_LENGTH, 1200);

        Map<String, Object> tumorAttributes = Maps.newHashMap();
        tumorAttributes.put(SPLIT_FRAGS, 10);
        tumorAttributes.put(STRAND_BIAS, 1);

        int posStart = 100;
        int posEnd = 200;

        Variant var = createSv(
                "01", CHR_1, CHR_1, posStart, posEnd, ORIENT_FWD, ORIENT_REV, "",
                commonAttributes, tumorAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertTrue(var.filters().contains(FilterType.SBX_ARTEFACT));

        // test 2: not applied for longer indels
        var = createSv(
                "01", CHR_1, CHR_1, posStart, posEnd + SBX_HEURISTIC_SHORT_LENGTH, ORIENT_REV, ORIENT_FWD, "",
                commonAttributes, tumorAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertFalse(var.filters().contains(FilterType.SBX_ARTEFACT));

        // test 3: penalty for breakends away from original junctions
        tumorAttributes.put(STRAND_BIAS, 0.5);

        commonAttributes.put(ASM_INFO, "2:100:1_2:2000:-1"); // note diff chromosome

        var = createSv(
                "01", CHR_1, CHR_1, posStart, posEnd, ORIENT_REV, ORIENT_FWD, "",
                commonAttributes, tumorAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertTrue(var.filters().contains(FilterType.SBX_ARTEFACT));

        // very short DUPs
        commonAttributes.put(ASM_INFO, "1:100:-1_1:130:1"); // matching

        var = createSv(
                "01", CHR_1, CHR_1, posStart, posStart + 30, ORIENT_REV, ORIENT_FWD, "",
                commonAttributes, tumorAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertTrue(var.filters().contains(FilterType.SBX_ARTEFACT));

        // INVs
        tumorAttributes.put(SPLIT_FRAGS, 12);

        // test 1: short INV with long insert sequence
        String insSequence = "G".repeat(21);

        var = createSv(
                "01", CHR_1, CHR_1, posStart, posStart + 60, ORIENT_FWD, ORIENT_FWD, insSequence,
                commonAttributes, tumorAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertTrue(var.filters().contains(FilterType.SBX_ARTEFACT));

        // test 2: zero-length INV
        var = createSv(
                "01", CHR_1, CHR_1, posStart, posStart, ORIENT_FWD, ORIENT_FWD, insSequence,
                commonAttributes, tumorAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertTrue(var.filters().contains(FilterType.SBX_ARTEFACT));

        // BNDs
        insSequence = "A".repeat(20);

        commonAttributes.put(ASM_INFO, "1:100:1_2:200:1"); // matching

        var = createSv(
                "01", CHR_1, CHR_2, posStart, posEnd, ORIENT_FWD, ORIENT_FWD, insSequence,
                commonAttributes, tumorAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertFalse(var.filters().contains(FilterType.SBX_ARTEFACT));

        insSequence = "G".repeat(20);

        var = createSv(
                "01", CHR_1, CHR_2, posStart, posEnd, ORIENT_FWD, ORIENT_FWD, insSequence,
                commonAttributes, tumorAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertTrue(var.filters().contains(FilterType.SBX_ARTEFACT));

        // SGLs
        commonAttributes.put(ASM_INFO, "1:100:1"); // matching

        var = createSv(
                "01", CHR_1, null, posStart, -1, ORIENT_FWD, ORIENT_FWD, "",
                commonAttributes, tumorAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertFalse(var.filters().contains(FilterType.SBX_ARTEFACT));

        // non-local and shorter assembly
        commonAttributes.put(ASM_INFO, "2:100:1"); // not matching
        commonAttributes.put(ASM_LENGTH, 400);

        var = createSv(
                "01", CHR_1, null, posStart, -1, ORIENT_FWD, ORIENT_FWD, "",
                commonAttributes, tumorAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertTrue(var.filters().contains(FilterType.SBX_ARTEFACT));
    }

    @Test
    public void testSbxZeroLengthInversions()
    {
        // test 1: non-line and stranded
        Map<String, Object> commonAttributes = Maps.newHashMap();

        Map<String, Object> tumorAttributes = Maps.newHashMap();
        tumorAttributes.put(SPLIT_FRAGS, 10);
        tumorAttributes.put(STRAND_BIAS, 1);

        int posStart = 100;

        Variant var = createSv(
                "01", CHR_1, CHR_1, posStart, posStart, ORIENT_FWD, ORIENT_FWD, "",
                commonAttributes, tumorAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertTrue(var.filters().contains(FilterType.SBX_INV_ZERO_MOTIF));

        // test 2: insert sequence with known motif
        String insSequence = "A" + SBX_INV_INSERT_MOTIFS.get(0);

        var = createSv(
                "01", CHR_1, CHR_1, posStart, posStart, ORIENT_FWD, ORIENT_FWD, insSequence,
                commonAttributes, tumorAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertTrue(var.filters().contains(FilterType.SBX_INV_ZERO_MOTIF));

        // reverse-complimented
        insSequence = Nucleotides.reverseComplementBases(SBX_INV_INSERT_MOTIFS.get(0)).substring(2) + "A";

        var = createSv(
                "01", CHR_1, CHR_1, posStart, posStart, ORIENT_FWD, ORIENT_FWD, insSequence,
                commonAttributes, tumorAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertTrue(var.filters().contains(FilterType.SBX_INV_ZERO_MOTIF));

        // handle mismatches: insert, SNV and delete
        String motif =  SBX_INV_INSERT_MOTIFS.get(1);
        insSequence = motif.substring(0, 10) + "G" + motif.substring(10, 15) + "A" + motif.substring(16, 20) + motif.substring(21, motif.length() - 2);

        var = createSv(
                "01", CHR_1, CHR_1, posStart, posStart, ORIENT_FWD, ORIENT_FWD, insSequence,
                commonAttributes, tumorAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertTrue(var.filters().contains(FilterType.SBX_INV_ZERO_MOTIF));

        // too many mismatches
        motif =  SBX_INV_INSERT_MOTIFS.get(2); // ACCCGTAGTATTAATGTATGTAATACTACGGGT
        insSequence = "TTTTT" + motif.substring(5, motif.length() - 6) + "GAAAAA";

        var = createSv(
                "01", CHR_1, CHR_1, posStart, posStart, ORIENT_FWD, ORIENT_FWD, insSequence,
                commonAttributes, tumorAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertFalse(var.filters().contains(FilterType.SBX_INV_ZERO_MOTIF));
    }

    @Ignore
    @Test
    public void testSbxInversionSequenceMatches()
    {
        List<String> testSequences = List.of(
                "ACCATCAATCCGGATGTATGCCGGATTGATGGT",
                "AGGAGTAACATCCATGTATGGGATGTTACTCCT",
                "AGGCGAGTATTCCATGTATGGGAATACTCGCCT",
                "ACCCGTAGTATTAATGTATGTAATACTACGGGT",
                "ACTATTGAAGGCTATGTATGAGCCTTCAATAGT",
                "AGGTATTGCCTGGATGTATGCCAGGCAATACCT",
                "AGAATCGCAGTTAATGTATGTAACTGCGATTCT",
                "CCATCAATCCGGATGTATGCCGGATTGATGG",
                "AGGATTATGGCAAATGTATGTTGCCATAATCCT",
                "CCCGTAGTATTAATGTATGTAATACTACGGG",
                "GGCGAGTATTCCATGTATGGGAATACTCGCC",
                "ACCTGCGGTTGAAATGTATGTTCAACCGCAGGT",
                "AGGATTATGGCAAAATGTATGTTGCCATAATCCT",
                "CTATTGAAGGCTATGTATGAGCCTTCAATAG",
                "GGAGTAACATCCATGTATGGGATGTTACTCC",
                "GGTATTGCCTGGATGTATGCCAGGCAATACC",
                "AGGCGAGTATTCCATACATGGAATACTCGCCT",
                "AGGTTATCAGCTTATGTATGAAGCTGATAACCT",
                "ACCTGCGGTTGAAAATGTATGTTCAACCGCAGGT",
                "GAATCGCAGTTAATGTATGTAACTGCGATTC",
                "GCGAGTATTCCATGTATGGGAATACTCGC",
                "GGATTATGGCAAATGTATGTTGCCATAATCC",
                "AGGGACATTACCAATGTATGTGGTAATGTCCCT",
                "CCTGCGGTTGAAATGTATGTTCAACCGCAGG",
                "GGATTATGGCAAAATGTATGTTGCCATAATCC",
                "AGGCGAGTATCCCATACATGGAATACTCGCCT",
                "GGCGAGTATTCCATACATGGAATACTCGCC",
                "ACCACTTCGCCAAATGTATGTTGGCGAAGTGGT",
                "AGGCTTACTACGGATGTATGCCGTAGTAAGCCT",
                "AGGATCTCCTATTATGTATGAATAGGAGATCCT",
                "AGGCGAGTATTCCATGTATGGAAATACTCGCCT",
                "AGGCGAGTATTCCATGTATGGGAAACTCGCCT",
                "AGGCGAGTATTCCATGTAGGGAATACTCGCCT",
                "AGGAGTAACATCCATACATGGATGTTACTCCT",
                "ACCCGTAGTATTAAATGTATGTAATACTACGGGT",
                "GTATTGCCTGGATGTATGCCAGGCAATAC",
                "GAGTAACATCCATGTATGGGATGTTACTC",
                "AGGAGTAACATCCCATACATGGGATGTTACTCCT",
                "AGGCGAGTATTCCCATACAGGGAATACTCGCCT",
                "GGTTATCAGCTTATGTATGAAGCTGATAACC",
                "CATCAATCCGGATGTATGCCGGATTGATG",
                "ACTATTGAAGGCTATGTATGAACCTTCAATAGT",
                "TCTTCTGTATTTCTAGTGTTCAATAGAAATACAGAAGA",
                "AATCGCAGTTAATGTATGTAACTGCGATT",
                "ACCATCAATCCGGAATGTATGCCGGATTGATGGT",
                "ACCATCAATCCGCATACATCCGGATTGATGGT",
                "CCTGCGGTTGAAAATGTATGTTCAACCGCAGG",
                "AGAATCGCAGTTAAATGTATGTAACTGCGATTCT",
                "AGGATTATGGCAAATGTATGTGCCATAATCCT",
                "GATTATGGCAAATGTATGTTGCCATAATC",
                "ACCCGTAGTATTAACATACATTAATACTACGGGT",
                "CCACTTCGCCAAATGTATGTTGGCGAAGTGG",
                "ACAACAGCCGAAGATGTATGCTTCGGCTGTTGT",
                "AGGAGTAACATCCCATACACGGATGTTACTCCT",
                "ACATAGAAGGTAGATGTATGCTACCTTCTATGT",
                "AGGAGTAACATCCATGTAATGGGATGTTACTCCT",
                "GGAGTAACATCCATGTATGAGATGTTACTCC");

        for(String testSequence : testSequences)
        {
            boolean hasMatch = false;
            String lowestMotifSeqMatch = "";
            int lowestMismatches = -1;

            String testSeqReversed = Nucleotides.reverseComplementBases(testSequence);

            for(String motifSeq : SBX_INV_INSERT_MOTIFS)
            {
                int lengthDiff = testSequence.length() - motifSeq.length();

                for(int i = 0; i <= abs(lengthDiff); ++i)
                {
                    int s1Offset = lengthDiff > 0 ? i : 0;
                    int s2Offset = lengthDiff < 0 ? i : 0;
                    int diffs = sequenceMismatches(testSequence, motifSeq, s1Offset, s2Offset);
                    diffs = min(diffs, sequenceMismatches(testSeqReversed, motifSeq, s1Offset, s2Offset));

                    if(diffs <= SBX_INV_INSERT_MOTIF_MAX_DIFFS)
                    {
                        hasMatch = true;
                        break;
                    }

                    if(lowestMismatches < 0 || diffs < lowestMismatches)
                    {
                        lowestMismatches = diffs;
                        lowestMotifSeqMatch = motifSeq;
                    }
                }

                if(hasMatch)
                    break;
            }

            if(!hasMatch)
                assertTrue(false);
        }
    }
}

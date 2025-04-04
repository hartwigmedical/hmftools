package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.pave.ImpactTestUtils.generateTestBases;
import static com.hartwig.hmftools.pave.VariantData.NO_LOCAL_PHASE_SET;

import static junit.framework.TestCase.assertEquals;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.pave.impact.ImpactClassifier;
import com.hartwig.hmftools.pave.impact.VariantTransImpact;

import org.junit.Test;

public class HgvsCodingAddressesTest
{
    // all exons and introns of length 12
    String exon1 = "CTAGGACACGAG";
    String intron1 = "CGAGGGCCCAAA";
    String exon2 = "CGACACGAGTAA";
    String intron2 = "GGGCCCAAATTT";
    String exon3 = "ATGAAGGAACCT";
    String intron3 = "AAAGGGCCCTTT";
    String exon4 = "AAGATGGAACCT";
    String intron4 = "AACCGGTTACGT";
    String exon5 = "GTACACAAGCCT";
    String intron6 = "CACCGGTTACGT";
    String exon6 = "ATACACAAGCCT";

    private final String refBases = generateTestBases(20)
            + exon1 + intron1 + exon2 + intron2 + exon3 + intron3 + exon4 + intron4 + exon5 + intron6 + exon6
            + generateTestBases(20);
    MockRefGenome refGenome = new MockRefGenome();
    private final ImpactClassifier classifier;
    private final TranscriptData transcriptData;
    private final TranscriptData transcriptDataRS;

    public HgvsCodingAddressesTest()
    {
        refGenome.RefGenomeMap.put(CHR_1, refBases);
        classifier = new ImpactClassifier(refGenome);
        //            21           33           45           57           69           81           93           105          117          129          141
        //            |   exon1    |   intron1  |   exon2    |   intron2  |   exon3    |   intron4  |   exon4    |   intron5  |   exon5    |   intron6  |   exon6
        // random 20  CTAGGACACGAG CGAGGGCCCAAA CGACACGAGTAA GGGCCCAAATTT ATGAAGGAACCT AAAGGGCCCTTT AAGATGGAACCT AACCGGTTACGT GTACACAAGCCT CACCGGTTACGT ATACACAAGCCT random 20
        //            |                         |          |              |                         |                         |                         |
        //  +         -24                       -12        -1             1                         13                        *1                        *13
        //  -         *24        *13            *12        *1             24         13             12         1              -1                        -13
        transcriptData = createTransExons(GENE_ID_1, TRANS_ID_1, POS_STRAND, new int[] { 21, 45, 69, 93, 117, 141 }, 11, 69, 104, false, "");
        transcriptDataRS = createTransExons(GENE_ID_1, TRANS_ID_1, NEG_STRAND, new int[] { 21, 45, 69, 93, 117, 141 }, 11, 69, 104, false, "");
    }

    @Test
    public void secondDownstreamExonRS()
    {
        checkDupRS(30,  2, "c.*13_*14dupCT");
//        checkDelRS(30,  2, "c.*13_*14delCT");
        checkInsRS(31,  "C", "c.*13_*14insG");
        checkSnvRS(31,  "A", "T", "c.*14T>A");

        checkSnvRS(32, "G", "A", "c.*13C>T");
    }

    @Test
    public void secondDownstreamIntronExonBoundaryRS()
    {
//        checkDupRS(31,  2, "c.*13-1_*13dupGC");
//        checkDelRS(31,  2, "c.*13-1_*13delGC");
//        checkInsRS(32,  "AA", "c.*13_*13+1insTT");
    }

    @Test
    public void secondDownstreamIntronRS()
    {
        checkDupRS(34,  4, "c.*13-6_*13-3dupCCCT");
//        checkDelRS(34,  4, "c.*13-6_*13-3delCCCT");
        checkInsRS(35,  "GG", "c.*13-4_*13-3insCC");
        checkSnvRS(35, "A", "C", "c.*13-3T>G");

        checkSnvRS(44, "A", "C", "c.*12+1T>G");
    }

    @Test
    public void downstreamExonIntronBoundaryRS()
    {
        checkDupRS(43,  2, "c.*12-0_*12+1dupGT"); // c.*12_*12+1dupGT
        checkDelRS(43,  2, "c.*12_*12+1delGT");
        checkInsRS(44, "CG",  "c.*12-0_*12+1insCG"); // c.*12_*12+1insCG ??
        checkSnvRS(44, "A", "G", "c.*12+1T>C");
    }

    @Test
    public void firstDownstreamExonRS()
    {
        checkDupRS(45,  2, "c.*10_*11dupTC");
        checkDelRS(45,  2, "c.*10_*11delTC");
        checkInsRS(45,  "TA", "c.*11_*12insTA");
        checkSnvRS(45, "C", "G", "c.*12G>C");

        checkDupRS(54,  2, "c.*1_*2dupTT");
        //        checkDelRS(54,  2, "c.*1_*2delTT");
        checkInsRS(55,  "CC", "c.*1_*2insGG");
        checkSnvRS(55, "A", "C", "c.*2T>G");

        checkDupRS(54,  2, "c.*1_*2dupTT");
//        checkDelRS(54,  2, "c.*1_*2delTT");
        checkInsRS(55,  "CC", "c.*1_*2insGG");
        checkSnvRS(55, "A", "C", "c.*2T>G");

        checkSnvRS(56, "A", "C", "c.*1T>G");
    }

    @Test
    public void firstDownstreamIntronRS()
    {
//        checkDupRS(64,  2, "c.24+3_24+4dupAT");
//        checkDelRS(64,  2, "c.24+3_24+4delAT");
//        checkInsRS(60,  "AAA", "c.22_23insTTT");
//        checkSnvRS(60, "T", "C", "c.23A>G");

//        checkSnvRS(68, "T", "C", "c.24+1A>G");
//        checkSnvRS(67, "T", "C", "c.24+2A>G");
    }

    @Test
    public void secondCodingExonRS()
    {
        checkDupRS(70,  2, "c.21_22dupTC");
        checkDelRS(70,  2, "c.21_22delTC");
        checkInsRS(70,  "AAA", "c.22_23insTTT");
        checkSnvRS(70, "T", "C", "c.23A>G");

//        checkDupRS(68,  2, "c.23_24dupAT");
//        checkDelRS(68,  2, "c.23_24delAT");
        checkInsRS(69,  "AAA", "c.23_24insTTT");
        checkSnvRS(69, "A", "C", "c.24T>G");
    }

    @Test
    public void intronExonBoundary2RS()
    {
//        checkDupRS(79, 2, "c.13-1_13dupTA");
        checkDelRS(79, 2, "c.13-1_13delTA");
        checkInsRS(79, "CAG",  "c.13_14insCTG"); // in the exon but useful comparison with ins at 80 below
//        checkInsRS(80, "CAG",  "c.13-1_13insCTG");
        checkSnvRS(80, "T", "C", "c.13A>G");
    }

    @Test
    public void firstIntronRS()
    {
        checkDupRS(89, 3, "c.12+1_12+3dupAAA");
        checkDelRS(89, 3, "c.12+1_12+3delAAA");
        checkInsRS(89, "AA",  "c.12+3_12+4insTT");
        checkSnvRS(89, "C", "T", "c.12+4G>A");

        checkDupRS(81, 2, "c.13-3_13-2dupTT");
        checkDelRS(81, 3, "c.13-4_13-2delCTT");
        checkInsRS(81, "CAG",  "c.13-2_13-1insCTG");
        checkSnvRS(81, "A", "C", "c.13-1T>G");
    }

    @Test
    public void exonBoundaryRS()
    {
        checkSnvRS(93, "A", "C", "c.12T>G");

        checkDelRS(91, 2, "c.12_12+1delTA");
        checkDupRS(91, 2, "c.12-0_12+1dupTA"); // c.12_12+1dupTA ?

        checkDelRS(90, 3, "c.12_12+2delTAA");
        checkDupRS(90, 3, "c.12-0_12+2dupTAA");

        checkDelRS(90, 4, "c.11_12+2delTTAA");
//        checkDupRS(90, 4, "c.11_12+2dupTTAA");
    }

    @Test
    public void firstCodingExonRS()
    {
        checkSnvRS(104, "T", "C", "c.1A>G");

        checkDupRS(102,  2, "c.1_2dupAG");
        checkDelRS(102,  2, "c.1_2delAG");
        checkInsRS(103,  "C", "c.1_2insG");
        checkSnvRS(103, "C", "G", "c.2G>C");
    }

    @Test
    public void firstUpstreamIntronRS()
    {
//        checkDupRS(114,  2, "c.-1+1_-1+2dupAC");
        checkDelRS(114,  2, "c.-1+1_-1+2delAC");
        checkInsRS(115,  "AA", "c.-1+1_-1+2insTT");
        checkSnvRS(116, "T", "G", "c.-1+1A>C");

//        checkDupRS(109,  2, "c.-1+6_1-6dupAC");
        checkDelRS(109,  2, "c.-1+6_-1+7delAC"); // c.-1+6_1-6delAC
        checkInsRS(110,  "AA", "c.-1+6_-1+7insTT");
        checkSnvRS(111, "T", "G", "c.-1+6A>C");

//        checkDupRS(107,  2, "c.1-5_1-4dupCG");
        checkDelRS(107,  2, "c.1-5_1-4delCG"); // c.-1+6_1-6delAC
        checkInsRS(108,  "AA", "c.1-5_1-4insTT");
        checkSnvRS(108, "C", "G", "c.1-4G>C");

//        checkDupRS(105,  2, "c.1-3_1-2dupGT");
        checkDelRS(105,  2, "c.1-3_1-2delGT");
        checkInsRS(106,  "G", "c.1-3_1-2insC");
        checkSnvRS(103, "C", "G", "c.2G>C");

        checkSnvRS(105, "A", "G", "c.1-1T>C");
    }

    @Test
    public void firstUpstreamExonIntronBoundaryRS()
    {
//        checkDupRS(115,  2, "c.-1_-1+1dupCA");
        checkDelRS(115,  2, "c.-1_-1+1delCA");
        checkInsRS(116,  "CC", "c.-1-0_-1+1insGG"); // c.-1_-1+1insGG
        checkSnvRS(117, "G", "A", "c.-1C>T");
    }

    @Test
    public void firstUpstreamExonRS()
    {
        checkDupRS(120,  2, "c.-6_-5dupGT");
        checkDelRS(120,  2, "c.-6_-5delGT");
        checkInsRS(120,  "G", "c.-5_-4insC");
        checkSnvRS(120, "C", "A", "c.-4G>T");
    }

    @Test
    public void secondUpstreamExonIntronBoundaryRS()
    {
//        checkDupRS(127,  2, "c.-12-1_-12dupGA");
        checkDelRS(127,  2, "c.-12-1_-12delGA");
//        checkInsRS(128,  "A", "c.-12-1_-12insT");
        checkInsRS(127,  "A", "c.-12_-11insT");
        checkSnvRS(128, "T", "A", "c.-12A>T");
    }

    @Test
    public void secondUpstreamIntronRS()
    {
//        checkDupRS(137,  2, "c.-13+2_-13+3dupCG");
        checkDelRS(137,  2, "c.-13+2_-13+3delCG");
        checkInsRS(137,  "G", "c.-13+3_-13+4insC");
        checkSnvRS(137, "A", "G", "c.-13+4T>C");

//        checkDupRS(132,  2, "c.-12-6_-12-5dupCC");
//        checkDelRS(132,  2, "c.-12-6_-12-5delCC");
        checkInsRS(132,  "AA", "c.-12-5_-12-4insTT");
        checkSnvRS(132, "C", "A", "c.-12-4G>T");

        checkSnvRS(129, "C", "A", "c.-12-1G>T");
    }

    @Test
    public void thirdUpstreamBoundaryRS()
    {
//        checkDupRS(139,  2, "c.-13_-13+1dupTA");
        checkDelRS(139,  2, "c.-13_-13+1delTA"); // c.-13_-13+1delTA ?
        checkInsRS(140,  "G", "c.-13-0_-13+1insC"); // ?
    }

    @Test
    public void secondUpstreamExonRS()
    {
        checkDupRS(150,  2, "c.-24_-23dupAG");
//        checkDelRS(150,  2, "c.-24_-23delAG");
        checkInsRS(150,  "A", "c.-23_-22insT");
        checkSnvRS(150, "C", "G", "c.-22G>C");

        checkDupRS(141,  2, "c.-15_-14dupTA");
        checkDelRS(141,  2, "c.-15_-14delTA");
        checkInsRS(141,  "G", "c.-14_-13insC");
        checkSnvRS(141, "A", "G", "c.-13T>C");
    }

    @Test
    public void secondDownstreamIntronAndExon()
    {
        // boundary
//        checkDupAt(127, 2, "c.*12_*12+1dupTC");
        checkDel(127, 2, "c.*12_*12+1delTC");
        checkIns(128, "GG", "c.*12-0_*12+1insGG"); // c.*12_*12+1insGG

        // in intron6
        checkDup(128, 2, "c.*12+1_*12+2dupCA");
        checkDel(128, 2, "c.*12+1_*12+2delCA");
        checkIns(129, "GG",  "c.*12+1_*12+2insGG");
        checkSnv(129, "C", "G", "c.*12+1C>G");

        checkDup(133, 2, "c.*12+6_*12+7dupGT");
        checkDel(133, 2, "c.*12+6_*12+7delGT");
        checkIns(134, "AA",  "c.*12+6_*12+7insAA");
        checkSnv(134, "G", "A", "c.*12+6G>A");

//        checkDupAt(134, 2, "c.*13-6_*13-5dupTT");
//        checkDelAt(134, 2, "c.*13-6_*13-5delTT");
        checkIns(135, "AA",  "c.*13-6_*13-5insAA");
        checkSnv(135, "T", "A", "c.*13-6T>A");

        // exon6
        checkDup(140, 2, "c.*13_*14dupAT");
//        checkDelAt(140, 2, "c.*13_*14delAT");
        checkIns(141, "CC",  "c.*13_*14insCC");
        checkSnv(141, "A", "T", "c.*13A>T");

        checkDup(141, 2, "c.*14_*15dupTA");
        checkDel(141, 2, "c.*14_*15delTA");
        checkIns(142, "CC",  "c.*14_*15insCC");
        checkSnv(142, "T", "G", "c.*14T>G");

        checkDup(150, 2, "c.*23_*24dupCT");
        checkDel(150, 2, "c.*23_*24delCT");
        checkIns(151, "AA",  "c.*23_*24insAA");
        checkSnv(151, "C", "G", "c.*23C>G");

        checkSnv(152, "T", "G", "c.*24T>G");
    }

    @Test
    public void downstreamExon()
    {
        //        checkDupAt(116, 2, "c.*1_*2dupGT");
        //        checkDelAt(116, 2, "c.*1_*2delGT");
        checkIns(117, "CC", "c.*1_*2insCC");
        checkSnv(117, "G", "C", "c.*1G>C");

        checkDup(119, 2, "c.*4_*5dupCA");
        checkDel(119, 2, "c.*4_*5delCA");
        checkIns(120, "GG", "c.*4_*5insGG");
        checkSnv(120, "C", "A", "c.*4C>A");

        checkDup(126, 2, "c.*11_*12dupCT");
        checkDel(126, 2, "c.*11_*12delCT");
        checkIns(127, "GG", "c.*11_*12insGG");
        checkSnv(127, "C", "A", "c.*11C>A");

        checkDup(127, 1, "c.*12dupT");
        checkDel(127, 1, "c.*12delT");
        checkIns(128, "GG", "c.*12-0_*12+1insGG"); // c.*12_*12+1insGG ?
        checkSnv(128, "T", "A", "c.*12T>A");
    }

    @Test
    public void downnstreamIntronExonBoundary()
    {
        checkDup(115, 2, "c.*1-1_*1-0dupTG");
        //        checkDelAt(115, 2, "c.*1-1_*1delTG");
        //        checkInsAt(116, "CC", "c.*1-1_*1insCC");
    }

    @Test
    public void downstreamIntron()
    {
        checkDup(104, 2, "c.24+1_24+2dupAA");
        checkDel(104, 2, "c.24+1_24+2delAA");
        checkIns(105, "CC", "c.24+1_24+2insCC");
        checkSnv(105, "A", "C", "c.24+1A>C");

        checkDup(109, 2, "c.24+6_24+7dupGT"); //c.24+6_c.*1-6dupGT?
        checkDel(109, 2, "c.24+6_24+7delGT");
        checkIns(110, "CC", "c.24+6_24+7insCC");
        checkSnv(110, "G", "C", "c.24+6G>C");

        //        checkDupAt(109, 2, "c.*1-2_*1-1dupGT");
        checkDel(114, 2, "c.*1-2_*1-1delGT");
        checkIns(115, "CC", "c.*1-2_*1-1insCC");
        checkSnv(115, "G", "C", "c.*1-2G>C");
    }

    @Test
    public void secondExonDownstreamIntronBoundary()
    {
        //        checkDupAt(103, 2, "c.24_24+1dupTA");
        //        checkDelAt(103, 2, "c.24_24+1delTA");
        checkIns(104, "GG", "c.24-0_24+1insGG"); // c.24_24+1insGG?
    }

    @Test
    public void secondCodingExon()
    {
        checkDup(102, 2, "c.23_24dupCT");
        //        checkDelAt(102, 2, "c.23_24delCT");
        checkIns(103, "A", "c.23_24insA");
        checkSnv(103, "C", "G", "c.23C>G");

        checkDup(93, 3, "c.14_16dupAGA");
        checkDel(93, 3, "c.14_16delAGA");
        checkIns(94, "C", "c.14_15insC");
        checkSnv(94, "A", "C", "c.14A>C");

        //        checkDupAt(92, 3, "c.13_15dupAAG");
        checkDel(92, 3, "c.13_15delAAG");
        checkIns(93, "C", "c.13_14insC");
        checkSnv(93, "A", "C", "c.13A>C");
    }

    @Test
    public void firstIntronSecondExonBoundary()
    {
        checkDup(90, 3, "c.13-2_13-0dupTTA"); // c.13-2_13dupTTA ?
        checkDel(90, 3, "c.13-2_13delTTA");
    }

    @Test
    public void firstIntron()
    {
        checkDup(90, 2, "c.13-2_13-1dupTT");
        checkDel(90, 2, "c.13-2_13-1delTT");
        checkIns(91, "AC", "c.13-2_13-1insAC");
        checkSnv(92, "T", "A", "c.13-1T>A");

        //        checkDupAt(86, 3, "c.13-6_13-4dupCCC");
        //        checkDelAt(86, 3, "c.13-6_13-4delCCC");
        checkIns(87, "TT", "c.13-6_13-5insTT");
        checkSnv(87, "C", "T", "c.13-6C>T");

        checkDup(85, 3, "c.12+6_12+8dupGCC");
        checkDel(85, 3, "c.12+6_12+8delGCC");
        checkIns(86, "TT", "c.12+6_12+7insTT");
        checkSnv(86, "G", "T", "c.12+6G>T");

        checkDup(84, 3, "c.12+5_12+7dupGGC");
        checkDel(84, 3, "c.12+5_12+7delGGC");
        checkIns(85, "TT", "c.12+5_12+6insTT");
        checkSnv(85, "G", "T", "c.12+5G>T");

        checkDup(80, 2, "c.12+1_12+2dupAA");
        checkDel(80, 2, "c.12+1_12+2delAA");
        checkIns(81, "CC", "c.12+1_12+2insCC");
        checkSnv(81, "A", "C", "c.12+1A>C");
    }

    @Test
    public void firstExonFirstIntronBoundary()
    {
        //        checkDupAt(79, 2, "c.12_12+1dupTT");
        checkDel(79, 3, "c.12_12+2delTAA");
        checkIns(80, "CC", "c.12-0_12+1insCC"); // c.12_12+1insCC ??
    }

    @Test
    public void firstCodingExon()
    {
        checkDup(79, 1, "c.12dupT");
        checkDel(79, 1, "c.12delT");
        checkSnv(80, "T", "A", "c.12T>A");

        checkDup(78, 2, "c.11_12dupCT");
        checkDel(78, 2, "c.11_12delCT");
        checkIns(79, "AA", "c.11_12insAA");
        checkSnv(79, "C", "T", "c.11C>T");

        checkDup(76, 2, "c.9_10dupAC");
        checkDel(76, 2, "c.9_10delAC");
        checkIns(77, "TTT", "c.9_10insTTT");
        checkSnv(77, "A", "C", "c.9A>C");

        checkDup(71, 2, "c.4_5dupAA");
        checkDel(71, 2, "c.4_5delAA");
        checkIns(72, "TTT", "c.4_5insTTT");
        checkSnv(72, "A", "C", "c.4A>C");

        checkDup(69, 2, "c.2_3dupTG");
        checkDel(69, 2, "c.2_3delTG");
        checkIns(70, "C", "c.2_3insC");
        checkSnv(70, "T", "C", "c.2T>C");

        //        checkDupAt(68, 3, "c.1_3dupATG");
        checkDel(68, 3, "c.1_3delATG");
        checkIns(69, "CC", "c.1_2insCC");
        checkSnv(69, "A", "C", "c.1A>C");
        checkDel(68, 1, "c.1delA");
    }

    @Test
    public void firstUpstreamIntron()
    {
        //        checkDupAt(66, 2, "c.1-2_1-1dupTT");
        checkDel(66, 2, "c.1-2_1-1delTT");
        //        checkDupAt(64, 4, "c.1-4_1-1dupATTT");
        checkDel(64, 4, "c.1-4_1-1delATTT");
        //        checkDupAt(63, 3, "c.1-5_1-3dupAAT");
        checkDel(63, 3, "c.1-5_1-3delAAT");
        checkDel(63, 1, "c.1-5delA");
        checkDel(62, 1, "c.-1+7delA");
        checkDel(61, 1, "c.-1+6delC");
        checkDel(59, 1, "c.-1+4delC");
        checkDel(58, 1, "c.-1+3delG");
        checkDel(57, 1, "c.-1+2delG");
        checkDel(56, 1, "c.-1+1delG");
    }

    @Test
    public void firstUpstreamExon()
    {
        checkDel(55, 1, "c.-1delA");
        //        checkDupAt(62, 3, "c.-1+6_-1+8dupAAA"); // -1+6_1-4 ??
        //        checkDelAt(62, 3, "c.-1+6_-1+8delAAA");
        //        checkDupAt(59, 3, "c.-1+3_-1+5dupCCC");
        //        checkDelAt(59, 3, "c.-1+3_-1+5delCCC");
        //        checkDupAt(57, 3, "c.-1+1_-1+3dupGGC");
        //        checkDelAt(57, 3, "c.-1+1_-1+3delGGC");
        //        checkDupAt(56, 3, "c.-1+0_-1+2dupGGG");
        //        checkDelAt(56, 3, "c.-1_-1+2delGGG");
        //        checkDupAt(55, 2, "c.-2_-1dupAA");
        checkDel(55, 1, "c.-1delA");
        checkDel(55, 2, "c.-1_-1+1delAG");
        //        checkDupAt(54, 2, "c.-2_-1dupAA");
        checkDel(54, 2, "c.-2_-1delAA");
        //        checkDupAt(46, 2, "c.-11_-10dupAC");
        checkDel(46, 2, "c.-10_-9delAC");
        //        checkDupAt(45, 2, "c.-12_-11dupAC");
        checkDel(45, 2, "c.-11_-10delGA");
        checkDel(45, 1, "c.-11delG");
        //        checkDelAt(44, 1, "c.-12delC");
    }

    @Test
    public void secondUpstreamIntron()
    {
        checkDel(43, 1, "c.-12-1delA");
        checkDel(40, 4, "c.-12-4_-12-1delCAAA");
        checkDel(36, 2, "c.-13+5_-13+6delGG");
        checkDel(32, 2, "c.-13+1_-13+2delCG");
    }

    @Test
    public void secondUpstreamExon()
    {
        checkDel(31, 2, "c.-13_-13+1delGC");
        checkDel(27, 3, "c.-17_-15delACG");
        checkDel(23, 3, "c.-21_-19delGGA");
        checkDel(21, 1, "c.-23delT");
        //        checkDelAt(20, 1, "c.-24delC");
    }

    private void checkDup(int position, int length, String expected)
    {
        checkVariantImpact(expected, createDupVariant(refBases, position, length));
    }

    private void checkDupRS(int position, int length, String expected)
    {
        checkVariantImpact(expected, createDupVariant(refBases, position, length), false);
    }

    private void checkSnv(int position, String ref, String alt, String expected)
    {
        checkSnv(position, ref, alt, expected, true);
    }

    private void checkSnvRS(int position, String ref, String alt, String expected)
    {
        checkSnv(position, ref, alt, expected, false);
    }

    private void checkSnv(int position, String ref, String alt, String expected, boolean forwardStrand)
    {
        checkBaseAt(position, ref);
        checkVariantImpact(expected, createSnvVariant(position, ref, alt), forwardStrand);
    }

    private void checkBaseAt(int position, String ref)
    {
        assertEquals(refBases.substring(position - 1, position), ref);
    }

    private VariantData createDupVariant(String refBases, int position, int length)
    {
        int pos = position - 1;
        String ref = refBases.substring(pos, pos + 1);
        String alt = refBases.substring(pos, pos + length + 1);
        VariantData var = new VariantData(CHR_1, position, ref, alt);
        String altBases = alt.substring(1);
        var.setVariantDetails(NO_LOCAL_PHASE_SET, altBases, altBases, 1);
        return var;
    }

    private VariantData createSnvVariant(int position, String ref, String alt)
    {
        return new VariantData(CHR_1, position, ref, alt);
    }

    private void checkIns(int position, String bases, String expected)
    {
        checkVariantImpact(expected, createInsVariant(position, bases));
    }

    private void checkInsRS(int position, String bases, String expected)
    {
        checkVariantImpact(expected, createInsVariant(position, bases), false);
    }

    private void checkDel(int position, int length, String expected)
    {
        checkVariantImpact(expected, createDelVariant(refBases, position, length));
    }

    private void checkDelRS(int position, int length, String expected)
    {
        checkVariantImpact(expected, createDelVariant(refBases, position, length), false);
    }

    private void checkVariantImpact(final String expected, final VariantData var)
    {
        checkVariantImpact(expected,var, true);
    }

    private void checkVariantImpact(String expected, VariantData var, boolean forwardStrand)
    {
        TranscriptData td = forwardStrand ? transcriptData : transcriptDataRS;
        VariantTransImpact impact = classifier.classifyVariant(var, td);
        assertEquals(expected, impact.codingContext().Hgvs);
    }

    private VariantData createInsVariant(int pos, String bases)
    {
        String ref = refBases.substring(pos, pos + 1);
        String alt = ref + bases;
        return new VariantData(CHR_1, pos, ref, alt);
    }

    private VariantData createDelVariant(String refBases, int position, int length)
    {
        int pos = position - 1;
        String ref = refBases.substring(pos, pos + length + 1);
        String alt = refBases.substring(pos, pos + 1);
        return new VariantData(CHR_1, position, ref, alt);
    }
}

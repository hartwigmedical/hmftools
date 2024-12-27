package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.common.TestUtils.TEST_SAMPLE;
import static com.hartwig.hmftools.sage.common.TestUtils.buildCigarString;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadCounter;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;
import com.hartwig.hmftools.sage.common.VariantTier;
import com.hartwig.hmftools.sage.phase.AppendVariantPhaser;
import com.hartwig.hmftools.sage.phase.LpsReadCounts;
import com.hartwig.hmftools.sage.vcf.CandidateSerialisation;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class PhaseAppendTest
{
    private static final String REF_BASES =
        //             10        20        30        40        50        60        70        80        90        100
        //   012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901
            "CGCAATATTCGAGTGGAAGTGACTCGATTAATCCGGTGCGTTCGTCACCGCTGTCTATTAAGCCGGTGCGTTCGTCACCGCTGTCTATTTACCCGGTGCGTTCGTC";
        //                  A                   T                   C                   G

    private static final String ALT_BASES = REF_BASES.substring(0, 15) + "A"
            + REF_BASES.substring(16, 35) + "T" + REF_BASES.substring(36, 55) + "C" + REF_BASES.substring(56);

    private static RefSequence REF_SEQUENCE = new RefSequence(1, REF_BASES.getBytes());

    @Test
    public void testAppendPhasedEvidence()
    {
        AppendVariantPhaser variantPhaser = new AppendVariantPhaser();

        variantPhaser.initialise(new ChrBaseRegion(CHR_1, 1, 1000), TEST_SAMPLE);

        VariantReadContextBuilder builder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        int position = 15;
        SimpleVariant var1 = new SimpleVariant(CHR_1, position, "G", "A");

        String readBases = ALT_BASES.substring(1, 35);
        String readCigar = buildCigarString(readBases.length());
        SAMRecord read = buildSamRecord(1, readCigar, readBases);

        VariantReadContext readContext = builder.createContext(var1, read, 14, REF_SEQUENCE);
        Candidate candidate1 = new Candidate(VariantTier.PANEL, readContext, 1, 0);

        VariantContextBuilder variantContextBuilder = CandidateSerialisation.toContext(candidate1);
        VariantContext varContext1 = variantContextBuilder.make();

        position = 35;
        SimpleVariant var2 = new SimpleVariant(CHR_1, position, "G", "T");

        readBases = ALT_BASES.substring(16, 55);
        readCigar = buildCigarString(readBases.length());
        read = buildSamRecord(16, readCigar, readBases);

        readContext = builder.createContext(var2, read, 19, REF_SEQUENCE);
        Candidate candidate2 = new Candidate(VariantTier.PANEL, readContext, 1, 0);

        variantContextBuilder = CandidateSerialisation.toContext(candidate2);
        VariantContext varContext2 = variantContextBuilder.make();

        position = 55;
        SimpleVariant var3 = new SimpleVariant(CHR_1, position, "T", "C");

        readBases = ALT_BASES.substring(36, 80);
        readCigar = buildCigarString(readBases.length());
        read = buildSamRecord(36, readCigar, readBases);

        readContext = builder.createContext(var3, read, 19, REF_SEQUENCE);
        Candidate candidate3 = new Candidate(VariantTier.PANEL, readContext, 1, 0);

        variantContextBuilder = CandidateSerialisation.toContext(candidate3);
        VariantContext varContext3 = variantContextBuilder.make();

        position = 75;
        SimpleVariant var4 = new SimpleVariant(CHR_1, position, "C", "G");

        readBases = ALT_BASES.substring(56, 99);
        readCigar = buildCigarString(readBases.length());
        read = buildSamRecord(56, readCigar, readBases);

        readContext = builder.createContext(var4, read, 19, REF_SEQUENCE);
        Candidate candidate4 = new Candidate(VariantTier.PANEL, readContext, 1, 0);

        variantContextBuilder = CandidateSerialisation.toContext(candidate4);
        VariantContext varContext4 = variantContextBuilder.make();

        int lpsId1 = 1; // variants 1 & 2
        int lpsId2 = 2; // variants 1, 2 & 3
        int lpsId3 = 3; // variants 2, 3 & 4

        int[] lpsIds = new int[] {lpsId1, lpsId2};
        varContext1.getCommonInfo().putAttribute(LOCAL_PHASE_SET, lpsIds);

        lpsIds = new int[] {lpsId1, lpsId2, lpsId3};
        varContext2.getCommonInfo().putAttribute(LOCAL_PHASE_SET, lpsIds);

        lpsIds = new int[] {lpsId2, lpsId3};
        varContext3.getCommonInfo().putAttribute(LOCAL_PHASE_SET, lpsIds);

        lpsIds = new int[] {lpsId3};
        varContext4.getCommonInfo().putAttribute(LOCAL_PHASE_SET, lpsIds);

        List<Candidate> candidates = Lists.newArrayList(candidate1, candidate2, candidate3, candidate4);
        List<VariantContext> variantContexts = Lists.newArrayList(varContext1, varContext2, varContext3, varContext4);

        variantPhaser.registerLocalPhaseSets(candidates, variantContexts);

        ReadContextCounter rcCounter1 = createReadCounter(candidate1.readContext());
        ReadContextCounter rcCounter2 = createReadCounter(candidate2.readContext());
        ReadContextCounter rcCounter3 = createReadCounter(candidate3.readContext());
        ReadContextCounter rcCounter4 = createReadCounter(candidate4.readContext());

        // LPS #1
        // depth support only
        List<ReadContextCounter> altCounters = Lists.newArrayList();
        List<ReadContextCounter> refCounters = Lists.newArrayList(rcCounter1, rcCounter2);

        variantPhaser.registeredPhasedVariants(altCounters, refCounters);
        variantPhaser.registeredPhasedVariants(altCounters, refCounters);

        // alt support
        altCounters = Lists.newArrayList(rcCounter1, rcCounter2);
        refCounters.clear();

        variantPhaser.registeredPhasedVariants(altCounters, refCounters);

        // LPS #2

        // depth support
        refCounters = Lists.newArrayList(rcCounter1, rcCounter2, rcCounter3);
        altCounters.clear();

        variantPhaser.registeredPhasedVariants(altCounters, refCounters);

        // alt support
        altCounters = Lists.newArrayList(rcCounter1, rcCounter2, rcCounter3);
        refCounters.clear();

        variantPhaser.registeredPhasedVariants(altCounters, refCounters);
        variantPhaser.registeredPhasedVariants(altCounters, refCounters);

        // LPS #3

        // alt support
        altCounters = Lists.newArrayList(rcCounter2, rcCounter3, rcCounter4);
        refCounters.clear();

        variantPhaser.registeredPhasedVariants(altCounters, refCounters);
        variantPhaser.registeredPhasedVariants(altCounters, refCounters);

        // subset of LPS 3
        refCounters = Lists.newArrayList(rcCounter3, rcCounter4);
        altCounters.clear();

        variantPhaser.registeredPhasedVariants(altCounters, refCounters);

        altCounters = Lists.newArrayList(rcCounter3, rcCounter4);
        refCounters.clear();

        variantPhaser.registeredPhasedVariants(altCounters, refCounters);

        Map<Integer,LpsReadCounts> lpsReadCountMap = variantPhaser.sampleLpsReadCounts().get(TEST_SAMPLE);
        assertNotNull(lpsReadCountMap);

        LpsReadCounts lpsReadCounts = lpsReadCountMap.get(lpsId1);
        assertNotNull(lpsReadCounts);
        assertNull(lpsReadCounts.SubsetVariantCounts);
        assertEquals(3, lpsReadCounts.Depth);
        assertEquals(1, lpsReadCounts.AltSupport);

        lpsReadCounts = lpsReadCountMap.get(lpsId2);
        assertNotNull(lpsReadCounts);
        assertNull(lpsReadCounts.SubsetVariantCounts);
        assertEquals(3, lpsReadCounts.Depth);
        assertEquals(2, lpsReadCounts.AltSupport);

        lpsReadCounts = lpsReadCountMap.get(lpsId3);
        assertNotNull(lpsReadCounts);
        assertNotNull(lpsReadCounts.SubsetVariantCounts);

        assertEquals(2, lpsReadCounts.Depth);
        assertEquals(2, lpsReadCounts.AltSupport);

        LpsReadCounts subLpsReadCounts = lpsReadCounts.SubsetVariantCounts.get(0);
        assertEquals(2, subLpsReadCounts.Depth);
        assertEquals(1, subLpsReadCounts.AltSupport);

        // variantPhaser.populateLocalPhaseSetInfo(candidates, variantContexts);
    }


}

package com.hartwig.hmftools.sage;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.sage.common.TestUtils.findInitialCandidates;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.quality.UltimaLocalRealigner;
import com.hartwig.hmftools.sage.testutil.MutatedBases;
import com.hartwig.hmftools.sage.testutil.MutatedBasesBuilder;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

// TODO: Move to better place.
public class UltimaLocalRealignmentIntegrationTest
{
    private static final int READ_DEPTH = 40;
    private static final int READ_LENGTH = 151;
    private static final int RANDOM_BASES_PADDING_LENGTH = 2 * READ_LENGTH;
    private static final String CHROMOSOME = CHR_1;

    // TODO: Rename
    @Test
    public void testCoreExtensionToIncludeHomopolymers()
    {
        Random random = new Random(0);

        String initRefBases = "A".repeat(RANDOM_BASES_PADDING_LENGTH);
        String middleRefBases = "ACTCGGCCTCCCGGAGGTGCCGGA";
        String finalRefBases = "A".repeat(RANDOM_BASES_PADDING_LENGTH);

        String refBases = initRefBases + middleRefBases + finalRefBases;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHROMOSOME, refBases);
        refGenome.ChromosomeLengths.put(CHROMOSOME, refBases.length());
        ChrBaseRegion region = new ChrBaseRegion(CHROMOSOME, 1, refBases.length());

        MutatedBasesBuilder mutatedBasesBuilder = new MutatedBasesBuilder(refBases);
        mutatedBasesBuilder.delBases(initRefBases.length() + 10, 1);
        mutatedBasesBuilder.mutateBase(initRefBases.length() + 13, 'T');
        MutatedBases mutatedBases = mutatedBasesBuilder.build();

        int leftMiddleMutIndex = mutatedBases.leftMutIndexFromRefPos(initRefBases.length());
        int rightMiddleMutIndex = mutatedBases.rightMutIndexFromRefPos(refBases.length() - finalRefBases.length() + 1);
        int centreMiddleMutIndex = (leftMiddleMutIndex + rightMiddleMutIndex)/2;
        int mutStartIndex = centreMiddleMutIndex - READ_LENGTH/2;

        SAMRecord read = mutatedBases.getUnpairedRead("READ_001", CHROMOSOME, mutStartIndex, READ_LENGTH, true);

        List<SAMRecord> reads = Lists.newArrayList();
        while(reads.size() < READ_DEPTH)
        {
            reads.add(read);
        }

        SageConfig config = new SageConfig(false, SequencingType.ULTIMA, null);
        List<Candidate> initialCandidates = findInitialCandidates(config, region, refGenome, reads);

        assertEquals(2, initialCandidates.size());

        List<Candidate> indelCandidates = initialCandidates.stream().filter(x -> x.variant().isIndel()).collect(Collectors.toList());
        List<Candidate> snvCandidates = initialCandidates.stream().filter(x -> !x.variant().isIndel()).collect(Collectors.toList());

        assertEquals(1, indelCandidates.size());
        assertEquals(1, snvCandidates.size());

        Candidate snvCandidate = snvCandidates.get(0);
        VariantReadContext snvReadContext = snvCandidate.readContext();

        assertTrue(UltimaLocalRealigner.isVariantExplained(snvReadContext));
    }

    // TODO: Rename
    @Test
    public void testSnvCandidate()
    {
        Random random = new Random(0);

        String initRefBases = "A".repeat(RANDOM_BASES_PADDING_LENGTH);
        String middleRefBases = "ACTCGGCCTCCCGGAGGTGCCGGA";
        String finalRefBases = "A".repeat(RANDOM_BASES_PADDING_LENGTH);

        String refBases = initRefBases + middleRefBases + finalRefBases;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHROMOSOME, refBases);
        refGenome.ChromosomeLengths.put(CHROMOSOME, refBases.length());
        ChrBaseRegion region = new ChrBaseRegion(CHROMOSOME, 1, refBases.length());

        MutatedBasesBuilder mutatedBasesBuilder = new MutatedBasesBuilder(refBases);
        // TODO: change this back once we get read base expansion working.
        // mutatedBasesBuilder.delBases(initRefBases.length() + 10, 1);
        mutatedBasesBuilder.delBases(initRefBases.length() + 11, 1);
        mutatedBasesBuilder.mutateBase(initRefBases.length() + 13, 'T');
        MutatedBases mutatedBases = mutatedBasesBuilder.build();

        int leftMiddleMutIndex = mutatedBases.leftMutIndexFromRefPos(initRefBases.length());
        int rightMiddleMutIndex = mutatedBases.rightMutIndexFromRefPos(refBases.length() - finalRefBases.length() + 1);
        int centreMiddleMutIndex = (leftMiddleMutIndex + rightMiddleMutIndex)/2;
        int mutStartIndex = centreMiddleMutIndex - READ_LENGTH/2;

        SAMRecord read = mutatedBases.getUnpairedRead("READ_001", CHROMOSOME, mutStartIndex, READ_LENGTH, true);

        List<SAMRecord> reads = Lists.newArrayList();
        while(reads.size() < READ_DEPTH)
        {
            reads.add(read);
        }

        SageConfig config = new SageConfig(false, SequencingType.ULTIMA, null);
        List<Candidate> initialCandidates = findInitialCandidates(config, region, refGenome, reads);

        assertEquals(2, initialCandidates.size());

        List<Candidate> indelCandidates = initialCandidates.stream().filter(x -> x.variant().isIndel()).collect(Collectors.toList());
        List<Candidate> snvCandidates = initialCandidates.stream().filter(x -> !x.variant().isIndel()).collect(Collectors.toList());

        assertEquals(1, indelCandidates.size());
        assertEquals(1, snvCandidates.size());

        Candidate snvCandidate = snvCandidates.get(0);
        VariantReadContext snvReadContext = snvCandidate.readContext();

        assertTrue(UltimaLocalRealigner.isVariantExplained(snvReadContext));
    }

    // TODO: Rename
    @Test
    public void testSingleBaseHomopolymerContractionInMiddle()
    {
        Random random = new Random(0);

        String initRefBases = "A".repeat(RANDOM_BASES_PADDING_LENGTH);
        String middleRefBases = "ACTCGGCCTCCCGGAGGTGCCGGA";
        String finalRefBases = "A".repeat(RANDOM_BASES_PADDING_LENGTH);

        String refBases = initRefBases + middleRefBases + finalRefBases;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHROMOSOME, refBases);
        refGenome.ChromosomeLengths.put(CHROMOSOME, refBases.length());
        ChrBaseRegion region = new ChrBaseRegion(CHROMOSOME, 1, refBases.length());

        MutatedBasesBuilder mutatedBasesBuilder = new MutatedBasesBuilder(refBases);
        // TODO: change this back once we get read base expansion working.
        //            mutatedBasesBuilder.delBases(initRefBases.length() + 10, 1);
        mutatedBasesBuilder.delBases(initRefBases.length() + 11, 1);
        mutatedBasesBuilder.mutateBase(initRefBases.length() + 13, 'T');
        MutatedBases mutatedBases = mutatedBasesBuilder.build();

        int leftMiddleMutIndex = mutatedBases.leftMutIndexFromRefPos(initRefBases.length());
        int rightMiddleMutIndex = mutatedBases.rightMutIndexFromRefPos(refBases.length() - finalRefBases.length() + 1);
        int centreMiddleMutIndex = (leftMiddleMutIndex + rightMiddleMutIndex)/2;
        int mutStartIndex = centreMiddleMutIndex - READ_LENGTH/2;

        SAMRecord read = mutatedBases.getUnpairedRead("READ_001", CHROMOSOME, mutStartIndex, READ_LENGTH, true);

        List<SAMRecord> reads = Lists.newArrayList();
        while(reads.size() < READ_DEPTH)
        {
            reads.add(read);
        }

        SageConfig config = new SageConfig(false, SequencingType.ULTIMA, null);
        List<Candidate> initialCandidates = findInitialCandidates(config, region, refGenome, reads);

        assertEquals(2, initialCandidates.size());

        List<Candidate> indelCandidates = initialCandidates.stream().filter(x -> x.variant().isIndel()).collect(Collectors.toList());
        List<Candidate> snvCandidates = initialCandidates.stream().filter(x -> !x.variant().isIndel()).collect(Collectors.toList());

        assertEquals(1, indelCandidates.size());
        assertEquals(1, snvCandidates.size());

        Candidate indelCandidate = indelCandidates.get(0);
        VariantReadContext indelReadContext = indelCandidate.readContext();

        assertTrue(UltimaLocalRealigner.isVariantExplained(indelReadContext));
    }

    // TODO: Rename
    @Test
    public void testSingleBaseHomopolymerContractionAtBoundary()
    {
        Random random = new Random(0);

        String initRefBases = "A".repeat(RANDOM_BASES_PADDING_LENGTH);
        String middleRefBases = "ACTCGGCCTCCCGGAGGTGCCGGA";
        String finalRefBases = "A".repeat(RANDOM_BASES_PADDING_LENGTH);

        String refBases = initRefBases + middleRefBases + finalRefBases;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHROMOSOME, refBases);
        refGenome.ChromosomeLengths.put(CHROMOSOME, refBases.length());
        ChrBaseRegion region = new ChrBaseRegion(CHROMOSOME, 1, refBases.length());

        MutatedBasesBuilder mutatedBasesBuilder = new MutatedBasesBuilder(refBases);
        mutatedBasesBuilder.delBases(initRefBases.length() + 10, 1);
        mutatedBasesBuilder.mutateBase(initRefBases.length() + 13, 'T');
        MutatedBases mutatedBases = mutatedBasesBuilder.build();

        int leftMiddleMutIndex = mutatedBases.leftMutIndexFromRefPos(initRefBases.length());
        int rightMiddleMutIndex = mutatedBases.rightMutIndexFromRefPos(refBases.length() - finalRefBases.length() + 1);
        int centreMiddleMutIndex = (leftMiddleMutIndex + rightMiddleMutIndex)/2;
        int mutStartIndex = centreMiddleMutIndex - READ_LENGTH/2;

        SAMRecord read = mutatedBases.getUnpairedRead("READ_001", CHROMOSOME, mutStartIndex, READ_LENGTH, true);

        List<SAMRecord> reads = Lists.newArrayList();
        while(reads.size() < READ_DEPTH)
        {
            reads.add(read);
        }

        SageConfig config = new SageConfig(false, SequencingType.ULTIMA, null);
        List<Candidate> initialCandidates = findInitialCandidates(config, region, refGenome, reads);

        assertEquals(2, initialCandidates.size());

        List<Candidate> indelCandidates = initialCandidates.stream().filter(x -> x.variant().isIndel()).collect(Collectors.toList());
        List<Candidate> snvCandidates = initialCandidates.stream().filter(x -> !x.variant().isIndel()).collect(Collectors.toList());

        assertEquals(1, indelCandidates.size());
        assertEquals(1, snvCandidates.size());

        Candidate indelCandidate = indelCandidates.get(0);
        VariantReadContext indelReadContext = indelCandidate.readContext();

        assertTrue(UltimaLocalRealigner.isVariantExplained(indelReadContext));
    }

    // TODO: Rename
    @Test
    public void testMultiBaseHomopolymerContractionInMiddle()
    {
        Random random = new Random(0);

        String initRefBases = "A".repeat(RANDOM_BASES_PADDING_LENGTH);
        String middleRefBases = "ACTCGGCCTCCCCCGGAGGTGCCGGA";
        String finalRefBases = "A".repeat(RANDOM_BASES_PADDING_LENGTH);

        String refBases = initRefBases + middleRefBases + finalRefBases;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHROMOSOME, refBases);
        refGenome.ChromosomeLengths.put(CHROMOSOME, refBases.length());
        ChrBaseRegion region = new ChrBaseRegion(CHROMOSOME, 1, refBases.length());

        MutatedBasesBuilder mutatedBasesBuilder = new MutatedBasesBuilder(refBases);
        // TODO: change this back once we get read base expansion working.
        //            mutatedBasesBuilder.delBases(initRefBases.length() + 10, 1);
        mutatedBasesBuilder.delBases(initRefBases.length() + 11, 3);
        mutatedBasesBuilder.mutateBase(initRefBases.length() + 15, 'T');
        MutatedBases mutatedBases = mutatedBasesBuilder.build();

        int leftMiddleMutIndex = mutatedBases.leftMutIndexFromRefPos(initRefBases.length());
        int rightMiddleMutIndex = mutatedBases.rightMutIndexFromRefPos(refBases.length() - finalRefBases.length() + 1);
        int centreMiddleMutIndex = (leftMiddleMutIndex + rightMiddleMutIndex)/2;
        int mutStartIndex = centreMiddleMutIndex - READ_LENGTH/2;

        SAMRecord read = mutatedBases.getUnpairedRead("READ_001", CHROMOSOME, mutStartIndex, READ_LENGTH, true);

        List<SAMRecord> reads = Lists.newArrayList();
        while(reads.size() < READ_DEPTH)
        {
            // TODO: Different read name for each read. Other tests as well.
            reads.add(read);
        }

        SageConfig config = new SageConfig(false, SequencingType.ULTIMA, null);
        List<Candidate> initialCandidates = findInitialCandidates(config, region, refGenome, reads);

        assertEquals(2, initialCandidates.size());

        List<Candidate> indelCandidates = initialCandidates.stream().filter(x -> x.variant().isIndel()).collect(Collectors.toList());
        List<Candidate> snvCandidates = initialCandidates.stream().filter(x -> !x.variant().isIndel()).collect(Collectors.toList());

        assertEquals(1, indelCandidates.size());
        assertEquals(1, snvCandidates.size());

        Candidate indelCandidate = indelCandidates.get(0);
        VariantReadContext indelReadContext = indelCandidate.readContext();

        assertTrue(UltimaLocalRealigner.isVariantExplained(indelReadContext));
    }

    // TODO: Rename
    @Test
    public void testMultiBaseHomopolymerContractionAtBoundary()
    {
        Random random = new Random(0);

        String initRefBases = "A".repeat(RANDOM_BASES_PADDING_LENGTH);
        String middleRefBases = "ACTCGGCCTCCCCCGGAGGTGCCGGA";
        String finalRefBases = "A".repeat(RANDOM_BASES_PADDING_LENGTH);

        String refBases = initRefBases + middleRefBases + finalRefBases;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHROMOSOME, refBases);
        refGenome.ChromosomeLengths.put(CHROMOSOME, refBases.length());
        ChrBaseRegion region = new ChrBaseRegion(CHROMOSOME, 1, refBases.length());

        MutatedBasesBuilder mutatedBasesBuilder = new MutatedBasesBuilder(refBases);
        mutatedBasesBuilder.delBases(initRefBases.length() + 10, 3);
        mutatedBasesBuilder.mutateBase(initRefBases.length() + 15, 'T');
        MutatedBases mutatedBases = mutatedBasesBuilder.build();

        int leftMiddleMutIndex = mutatedBases.leftMutIndexFromRefPos(initRefBases.length());
        int rightMiddleMutIndex = mutatedBases.rightMutIndexFromRefPos(refBases.length() - finalRefBases.length() + 1);
        int centreMiddleMutIndex = (leftMiddleMutIndex + rightMiddleMutIndex)/2;
        int mutStartIndex = centreMiddleMutIndex - READ_LENGTH/2;

        SAMRecord read = mutatedBases.getUnpairedRead("READ_001", CHROMOSOME, mutStartIndex, READ_LENGTH, true);

        List<SAMRecord> reads = Lists.newArrayList();
        while(reads.size() < READ_DEPTH)
        {
            // TODO: Different read name for each read. Other tests as well.
            reads.add(read);
        }

        SageConfig config = new SageConfig(false, SequencingType.ULTIMA, null);
        List<Candidate> initialCandidates = findInitialCandidates(config, region, refGenome, reads);

        assertEquals(2, initialCandidates.size());

        List<Candidate> indelCandidates = initialCandidates.stream().filter(x -> x.variant().isIndel()).collect(Collectors.toList());
        List<Candidate> snvCandidates = initialCandidates.stream().filter(x -> !x.variant().isIndel()).collect(Collectors.toList());

        assertEquals(1, indelCandidates.size());
        assertEquals(1, snvCandidates.size());

        Candidate indelCandidate = indelCandidates.get(0);
        VariantReadContext indelReadContext = indelCandidate.readContext();

        assertTrue(UltimaLocalRealigner.isVariantExplained(indelReadContext));
    }

    // TODO: Rename
    @Test
    public void testSingleHomopolymerDeletion()
    {
        Random random = new Random(0);

        String initRefBases = "A".repeat(RANDOM_BASES_PADDING_LENGTH);
        String middleRefBases = "ACTCGGCCTCCCGGAGGTGCCGGA";
        String finalRefBases = "A".repeat(RANDOM_BASES_PADDING_LENGTH);

        String refBases = initRefBases + middleRefBases + finalRefBases;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHROMOSOME, refBases);
        refGenome.ChromosomeLengths.put(CHROMOSOME, refBases.length());
        ChrBaseRegion region = new ChrBaseRegion(CHROMOSOME, 1, refBases.length());

        MutatedBasesBuilder mutatedBasesBuilder = new MutatedBasesBuilder(refBases);
        mutatedBasesBuilder.delBases(initRefBases.length() + 10, 3);
        mutatedBasesBuilder.mutateBase(initRefBases.length() + 13, 'T');
        MutatedBases mutatedBases = mutatedBasesBuilder.build();

        int leftMiddleMutIndex = mutatedBases.leftMutIndexFromRefPos(initRefBases.length());
        int rightMiddleMutIndex = mutatedBases.rightMutIndexFromRefPos(refBases.length() - finalRefBases.length() + 1);
        int centreMiddleMutIndex = (leftMiddleMutIndex + rightMiddleMutIndex)/2;
        int mutStartIndex = centreMiddleMutIndex - READ_LENGTH/2;

        SAMRecord read = mutatedBases.getUnpairedRead("READ_001", CHROMOSOME, mutStartIndex, READ_LENGTH, true);

        List<SAMRecord> reads = Lists.newArrayList();
        while(reads.size() < READ_DEPTH)
        {
            reads.add(read);
        }

        SageConfig config = new SageConfig(false, SequencingType.ULTIMA, null);
        List<Candidate> initialCandidates = findInitialCandidates(config, region, refGenome, reads);

        assertEquals(2, initialCandidates.size());

        List<Candidate> indelCandidates = initialCandidates.stream().filter(x -> x.variant().isIndel()).collect(Collectors.toList());
        List<Candidate> snvCandidates = initialCandidates.stream().filter(x -> !x.variant().isIndel()).collect(Collectors.toList());

        assertEquals(1, indelCandidates.size());
        assertEquals(1, snvCandidates.size());

        Candidate indelCandidate = indelCandidates.get(0);
        VariantReadContext indelReadContext = indelCandidate.readContext();

        assertTrue(UltimaLocalRealigner.isVariantExplained(indelReadContext));
    }

    // TODO: Rename
    @Test
    public void testMultipleHomopolymerDeletion()
    {
        Random random = new Random(0);

        String initRefBases = "A".repeat(RANDOM_BASES_PADDING_LENGTH);
        String middleRefBases = "ACTCGGCCTCCCTAAGGAGGTGCCGGA";
        String finalRefBases = "A".repeat(RANDOM_BASES_PADDING_LENGTH);

        String refBases = initRefBases + middleRefBases + finalRefBases;
        MockRefGenome refGenome = new MockRefGenome(true);
        refGenome.RefGenomeMap.put(CHROMOSOME, refBases);
        refGenome.ChromosomeLengths.put(CHROMOSOME, refBases.length());
        ChrBaseRegion region = new ChrBaseRegion(CHROMOSOME, 1, refBases.length());

        MutatedBasesBuilder mutatedBasesBuilder = new MutatedBasesBuilder(refBases);
        mutatedBasesBuilder.delBases(initRefBases.length() + 10, 6);
        mutatedBasesBuilder.mutateBase(initRefBases.length() + 16, 'T');
        MutatedBases mutatedBases = mutatedBasesBuilder.build();

        int leftMiddleMutIndex = mutatedBases.leftMutIndexFromRefPos(initRefBases.length());
        int rightMiddleMutIndex = mutatedBases.rightMutIndexFromRefPos(refBases.length() - finalRefBases.length() + 1);
        int centreMiddleMutIndex = (leftMiddleMutIndex + rightMiddleMutIndex)/2;
        int mutStartIndex = centreMiddleMutIndex - READ_LENGTH/2;

        SAMRecord read = mutatedBases.getUnpairedRead("READ_001", CHROMOSOME, mutStartIndex, READ_LENGTH, true);

        List<SAMRecord> reads = Lists.newArrayList();
        while(reads.size() < READ_DEPTH)
        {
            reads.add(read);
        }

        SageConfig config = new SageConfig(false, SequencingType.ULTIMA, null);
        List<Candidate> initialCandidates = findInitialCandidates(config, region, refGenome, reads);

        assertEquals(2, initialCandidates.size());

        List<Candidate> indelCandidates = initialCandidates.stream().filter(x -> x.variant().isIndel()).collect(Collectors.toList());
        List<Candidate> snvCandidates = initialCandidates.stream().filter(x -> !x.variant().isIndel()).collect(Collectors.toList());

        assertEquals(1, indelCandidates.size());
        assertEquals(1, snvCandidates.size());

        Candidate indelCandidate = indelCandidates.get(0);
        VariantReadContext indelReadContext = indelCandidate.readContext();

        assertTrue(UltimaLocalRealigner.isVariantExplained(indelReadContext));
    }
}

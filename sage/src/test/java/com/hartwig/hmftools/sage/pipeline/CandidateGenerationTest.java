//package com.hartwig.hmftools.sage.pipeline;
//
//import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
//import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
//import static com.hartwig.hmftools.sage.common.TestUtils.findInitialCandidates;
//
//import static org.junit.Assert.assertEquals;
//import static org.junit.Assert.assertFalse;
//
//import java.util.List;
//import java.util.Random;
//import java.util.stream.Collectors;
//
//import com.hartwig.hmftools.common.region.BaseRegion;
//import com.hartwig.hmftools.common.region.ChrBaseRegion;
//import com.hartwig.hmftools.common.test.MockRefGenome;
//import com.hartwig.hmftools.sage.candidate.Candidate;
//import com.hartwig.hmftools.sage.testutil.MutatedBases;
//import com.hartwig.hmftools.sage.testutil.MutatedBasesBuilder;
//
//import org.junit.Test;
//
//import htsjdk.samtools.SAMRecord;
//
//public class CandidateGenerationTest
//{
//    private static final int FRAGMENT_LENGTH = 600;
//    private static final int RANDOM_BASES_PADDING_LENGTH = 2 * FRAGMENT_LENGTH;
//    private static final int READ_DEPTH = 40;
//    private static final int READ_LENGTH = 151;
//    private static final String CHROMOSOME = CHR_1;
//
//    @Test
//    public void testNonOverlappingCoresAreNotExtended()
//    {
//        Random random = new Random(0);
//        int insertLength = 20;
//        int middleRefBasesLength = 100;
//
//        String initRefBases = generateRandomBases(random, RANDOM_BASES_PADDING_LENGTH);
//        String middleRefBases = generateRandomBases(random, middleRefBasesLength);
//        String finalRefBases = generateRandomBases(random, RANDOM_BASES_PADDING_LENGTH);
//
//        String refBases = initRefBases + middleRefBases + finalRefBases;
//        MockRefGenome refGenome = new MockRefGenome(true);
//        refGenome.RefGenomeMap.put(CHROMOSOME, refBases);
//        refGenome.ChromosomeLengths.put(CHROMOSOME, refBases.length());
//        ChrBaseRegion region = new ChrBaseRegion(CHROMOSOME, 1, refBases.length());
//
//        MutatedBasesBuilder mutatedBasesBuilder = new MutatedBasesBuilder(refBases);
//        mutatedBasesBuilder.insertBases(initRefBases.length() + 1, "A".repeat(insertLength));
//        mutatedBasesBuilder.mutateBase(initRefBases.length() + 21, 'T');
//        MutatedBases mutatedBases = mutatedBasesBuilder.build();
//
//        BaseRegion readRegion = new BaseRegion(
//                mutatedBases.leftMutIndexFromRefPos(initRefBases.length()),
//                mutatedBases.rightMutIndexFromRefPos(refBases.length() - finalRefBases.length() + 1));
//        List<SAMRecord> reads = mutatedBases.generateRandomReads(random, CHROMOSOME, READ_DEPTH, READ_LENGTH, FRAGMENT_LENGTH, readRegion);
//
//        List<Candidate> initialCandidates = findInitialCandidates(region, refGenome, reads);
//        List<Candidate> indelCandidates =
//                initialCandidates.stream().filter(candidate -> candidate.variant().isIndel()).collect(Collectors.toList());
//        List<Candidate> snvCandidates =
//                initialCandidates.stream().filter(candidate -> candidate.variant().isSNV()).collect(Collectors.toList());
//
//        assertEquals(2, initialCandidates.size());
//        assertEquals(1, indelCandidates.size());
//        assertEquals(1, snvCandidates.size());
//
//        assertFalse(indelCandidates.get(0).coreBaseRegion().overlaps(snvCandidates.get(0).coreBaseRegion()));
//    }
//}

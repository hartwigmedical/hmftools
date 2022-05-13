package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxFunction.FUSIONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.TRANSCRIPT_COUNTS;
import static com.hartwig.hmftools.isofox.ReadCountsTest.REF_BASE_STR_1;
import static com.hartwig.hmftools.isofox.TestUtils.CHR_1;
import static com.hartwig.hmftools.isofox.TestUtils.CHR_2;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_3;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_5;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.isofox.TestUtils.NEG_STRAND;
import static com.hartwig.hmftools.isofox.TestUtils.POS_STRAND;
import static com.hartwig.hmftools.isofox.TestUtils.TRANS_1;
import static com.hartwig.hmftools.isofox.TestUtils.addRacReadGroup;
import static com.hartwig.hmftools.isofox.TestUtils.addTestGenes;
import static com.hartwig.hmftools.isofox.TestUtils.addTestTranscripts;
import static com.hartwig.hmftools.isofox.TestUtils.createCigar;
import static com.hartwig.hmftools.isofox.TestUtils.createGeneCollection;
import static com.hartwig.hmftools.isofox.TestUtils.createGeneReadData;
import static com.hartwig.hmftools.isofox.TestUtils.createMappedRead;
import static com.hartwig.hmftools.isofox.TestUtils.createReadPair;
import static com.hartwig.hmftools.isofox.TestUtils.createReadRecord;
import static com.hartwig.hmftools.isofox.TestUtils.createSupplementaryReadPair;
import static com.hartwig.hmftools.isofox.TestUtils.populateRefGenome;
import static com.hartwig.hmftools.isofox.fusion.FusionConstants.HIGH_LOG_COUNT;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.MATCHED_JUNCTION;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGNED;
import static com.hartwig.hmftools.isofox.fusion.FusionTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.isofox.fusion.FusionTestUtils.createGroup;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.isofox.BamFragmentAllocator;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.BaseDepth;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.results.ResultsWriter;

import org.junit.Test;

public class FusionFiltersTest
{
    @Test
    public void testHardFilterCache()
    {
        HardFilteredCache hardFilteredCache = new HardFilteredCache();

        Map<String,Set<String>> hardFilteredReadIds = Maps.newHashMap();

        String CHR_X = "X";
        String chrPair1 = HardFilteredCache.formChromosomePairString(CHR_2, CHR_1);
        String chrPair2 = HardFilteredCache.formChromosomePairString(CHR_X, CHR_1);
        String chrPair3 = HardFilteredCache.formChromosomePairString(CHR_2, CHR_X);

        // read IDs indicating which chromosomes they link
        String readId1 = "1_2_01";
        String readId2 = "1_2_02";
        String readId3 = "1_X_01";
        String readId4 = "1_X_02";
        String readId5 = "2_X_01";
        String readId6 = "2_X_02";

        Set<String> readIds = Sets.newHashSet();
        readIds.add(readId1);
        readIds.add(readId2);
        hardFilteredReadIds.put(chrPair1, readIds);

        readIds = Sets.newHashSet();
        readIds.add(readId3);
        readIds.add(readId4);
        hardFilteredReadIds.put(chrPair2, readIds);

        hardFilteredCache.addHardFilteredReads(hardFilteredReadIds);

        assertEquals(4, hardFilteredCache.cacheCount());
        assertEquals(2, hardFilteredCache.chrPairCount());

        // now test 2's incomplete groups
        Map<String,Map<String,FusionReadGroup>> chrIncompleteGroups = Maps.newHashMap();
        Map<String,FusionReadGroup> otherChrGroups = Maps.newHashMap();
        otherChrGroups.put(readId1, createGroup(readId1));
        otherChrGroups.put(readId2, createGroup(readId2));
        chrIncompleteGroups.put(CHR_1, otherChrGroups);

        otherChrGroups = Maps.newHashMap();
        otherChrGroups.put(readId5, createGroup(readId5));
        otherChrGroups.put(readId6, createGroup(readId6));
        chrIncompleteGroups.put(CHR_X, otherChrGroups);

        // add 2's hard-filtered reads
        hardFilteredReadIds = Maps.newHashMap();

        readIds = Sets.newHashSet();
        readIds.add(readId1);
        readIds.add(readId2);

        hardFilteredReadIds.put(chrPair1, readIds);

        readIds = Sets.newHashSet();
        readIds.add(readId5);
        readIds.add(readId6);

        hardFilteredReadIds.put(chrPair3, readIds);

        assertEquals(4, chrIncompleteGroups.values().stream().mapToInt(x -> x.size()).sum());
        assertEquals(4, hardFilteredReadIds.values().stream().mapToInt(x -> x.size()).sum());

        // 2 of the new hard-filtered groups are now removed
        hardFilteredCache.removeHardFilteredReads(CHR_2, chrIncompleteGroups, hardFilteredReadIds);

        assertEquals(2, chrIncompleteGroups.values().stream().mapToInt(x -> x.size()).sum());
        assertEquals(2, hardFilteredReadIds.values().stream().mapToInt(x -> x.size()).sum());
        assertEquals(2, hardFilteredCache.cacheCount());

        hardFilteredCache.purgeChromosomeEntries(CHR_2, CHR_1);
        assertEquals(1, hardFilteredCache.chrPairCount());

        hardFilteredCache.addHardFilteredReads(hardFilteredReadIds);

        assertEquals(4, hardFilteredCache.cacheCount());

        chrIncompleteGroups.clear();
        otherChrGroups = Maps.newHashMap();
        otherChrGroups.put(readId3, createGroup(readId3));
        otherChrGroups.put(readId4, createGroup(readId4));
        chrIncompleteGroups.put(CHR_1, otherChrGroups);

        otherChrGroups = Maps.newHashMap();
        otherChrGroups.put(readId5, createGroup(readId5));
        otherChrGroups.put(readId6, createGroup(readId6));
        chrIncompleteGroups.put(CHR_2, otherChrGroups);

        // add X's hard-filtered reads
        hardFilteredReadIds = Maps.newHashMap();
        readIds = Sets.newHashSet();
        readIds.add(readId3);
        readIds.add(readId4);

        hardFilteredReadIds.put(chrPair2, readIds);

        readIds = Sets.newHashSet();
        readIds.add(readId5);
        readIds.add(readId6);

        hardFilteredReadIds.put(chrPair3, readIds);

        assertEquals(4, chrIncompleteGroups.values().stream().mapToInt(x -> x.size()).sum());

        // all of the new hard-filtered groups are now removed
        hardFilteredCache.removeHardFilteredReads(CHR_X, chrIncompleteGroups, hardFilteredReadIds);
        hardFilteredCache.purgeChromosomeEntries(CHR_X, CHR_1);
        hardFilteredCache.purgeChromosomeEntries(CHR_X, CHR_2);

        assertEquals(0, chrIncompleteGroups.values().stream().mapToInt(x -> x.size()).sum());
        assertEquals(0, hardFilteredReadIds.values().stream().mapToInt(x -> x.size()).sum());
        assertEquals(0, hardFilteredCache.cacheCount());

        hardFilteredCache.addHardFilteredReads(hardFilteredReadIds);
        assertEquals(0, hardFilteredCache.cacheCount());
        assertEquals(0, hardFilteredCache.chrPairCount());
    }

    private static FusionReadGroup createGroup(final String readId)
    {
        return new FusionReadGroup(readId, Lists.newArrayList());
    }

    @Test
    public void testFusionHardFilters()
    {
        final EnsemblDataCache geneTransCache = createGeneDataCache();

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        IsofoxConfig config = new IsofoxConfig();
        config.Functions.clear();
        config.Functions.add(FUSIONS);
        config.Fusions.MinHardFilterFrags = 2;

        populateRefGenome(config.RefGenome);

        // create supplementary read groups across chromosomes 1 & 2 where for fusion A it is hard-filtered on chr 1 but not on 2,
        // being left as a partial group, and then vice versa for fusion A on chr 2

        int gcId = 0;

        final GeneCollection gc3 =
                createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_3)));
        final GeneCollection gc5 =
                createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_5)));

        BamFragmentAllocator bamReader1 = new BamFragmentAllocator(config, new ResultsWriter(config));

        FusionTaskManager fusionTaskManager = new FusionTaskManager(config, geneTransCache);

        FusionFinder finderChr1 = fusionTaskManager.createFusionFinder(gc3.chromosome());

        int readId1 = 1;
        int readId2 = 2;
        int readId3 = 3;

        ReadRecord read1 = createMappedRead(readId1, gc5, 10210, 10249, createCigar(0, 40, 20));

        ReadRecord[] readPair1 = createSupplementaryReadPair(readId1, gc5, gc3, 10200, 10219, 20281, 20300,
                createCigar(20, 20, 0), createCigar(0, 20, 20), true);

        readPair1[0].setStrand(true, false);

        bamReader1.processReadRecords(gc5, Lists.newArrayList(read1, readPair1[0]));

        // supporting the first junction and enough to avoid being hard-filtered
        ReadRecord read2 = createMappedRead(readId2, gc5, 10210, 10249, createCigar(0, 40, 20));

        ReadRecord[] readPair2 = createSupplementaryReadPair(readId2, gc5, gc3, 10200, 10219, 20281, 20300,
                createCigar(20, 20, 0), createCigar(0, 20, 20), true);

        readPair2[0].setStrand(true, false);

        bamReader1.processReadRecords(gc5, Lists.newArrayList(read2, readPair2[0]));

        // single read for a new junction, hard-filtered - cannot be at a known splice site
        ReadRecord read3 = createMappedRead(readId3, gc5, 10410, 10449, createCigar(0, 40, 20));

        ReadRecord[] readPair3 = createSupplementaryReadPair(readId3, gc5, gc3, 10401, 10420, 20480, 20499,
                createCigar(20, 20, 0), createCigar(0, 20, 20), true);

        readPair3[0].setStrand(true, false);

        bamReader1.processReadRecords(gc5, Lists.newArrayList(read3, readPair3[0]));

        bamReader1.getChimericReadTracker().postProcessChimericReads(bamReader1.getBaseDepth(), bamReader1.getFragmentTracker());

        List<FusionReadGroup> completeGroups = finderChr1.processNewChimericReadGroups(
                gc5, bamReader1.getBaseDepth(), bamReader1.getChimericReadTracker().getReadMap());

        fusionTaskManager.addRacFragments(gc5.chromosome(), gc5.id(), bamReader1.getChimericReadTracker().extractJunctionRacFragments());

        finderChr1.processLocalReadGroups(completeGroups);

        assertEquals(2, finderChr1.getChimericPartialReadGroups().size());
        assertEquals(0, completeGroups.size());

        // finish chr 1
        assertEquals(1, bamReader1.getChimericReadTracker().getHardFilteredReadIds().values().stream().mapToInt(x -> x.size()).sum());

        Map<String,Map<String,FusionReadGroup>> chrIncompleteReadsGroups = finderChr1.extractIncompleteReadGroups(
                gc5.chromosome(), bamReader1.getChimericReadTracker().getHardFilteredReadIds());

        List<FusionReadGroup> interChromosomalGroups = fusionTaskManager.addIncompleteReadGroup(
                gc5.chromosome(), chrIncompleteReadsGroups, bamReader1.getChimericReadTracker().getHardFilteredReadIds());

        assertEquals(1, fusionTaskManager.hardFilteredCache().cacheCount());
        assertEquals(1, fusionTaskManager.hardFilteredCache().chrPairCount());

        finderChr1.processInterChromosomalReadGroups(interChromosomalGroups);

        assertTrue(finderChr1.getFusionCandidates().isEmpty());

        // now chromosome 2
        BamFragmentAllocator bamReader2 = new BamFragmentAllocator(config, new ResultsWriter(config));
        bamReader2.processReadRecords(gc3, Lists.newArrayList(readPair1[1]));
        bamReader2.processReadRecords(gc3, Lists.newArrayList(readPair2[1]));
        bamReader2.processReadRecords(gc3, Lists.newArrayList(readPair3[1]));

        bamReader2.getChimericReadTracker().postProcessChimericReads(bamReader2.getBaseDepth(), bamReader2.getFragmentTracker());

        FusionFinder finderChr2 = fusionTaskManager.createFusionFinder(gc3.chromosome());

        completeGroups = finderChr2.processNewChimericReadGroups(
                gc3, bamReader2.getBaseDepth(), bamReader2.getChimericReadTracker().getReadMap());

        fusionTaskManager.addRacFragments(gc3.chromosome(), gc3.id(), bamReader2.getChimericReadTracker().extractJunctionRacFragments());

        finderChr2.processLocalReadGroups(completeGroups);

        assertEquals(2, finderChr2.getChimericPartialReadGroups().size());
        assertEquals(0, completeGroups.size());

        // finish chr 2
        assertEquals(1, bamReader2.getChimericReadTracker().getHardFilteredReadIds().values().stream().mapToInt(x -> x.size()).sum());

        chrIncompleteReadsGroups = finderChr2.extractIncompleteReadGroups(
                gc3.chromosome(), bamReader2.getChimericReadTracker().getHardFilteredReadIds());

        // reconcile the 2 chromosomes:
        // the two hard-filtered refs to readId 3 are cancelled out, so there no remaining hard-filtered reads
        // reads 1 & 2 are merged from partial to complete groups
        interChromosomalGroups = fusionTaskManager.addIncompleteReadGroup(
                gc3.chromosome(), chrIncompleteReadsGroups, bamReader2.getChimericReadTracker().getHardFilteredReadIds());

        assertEquals(0, fusionTaskManager.incompleteReadGroups().values().stream().mapToInt(x -> x.size()).sum());
        assertEquals(0, fusionTaskManager.hardFilteredCache().cacheCount());
        assertEquals(0, fusionTaskManager.hardFilteredCache().chrPairCount());

        finderChr2.processInterChromosomalReadGroups(interChromosomalGroups);

        assertEquals(1, finderChr2.getFusionCandidates().values().stream().mapToInt(x -> x.size()).sum());
    }

}

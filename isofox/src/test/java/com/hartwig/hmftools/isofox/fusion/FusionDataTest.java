package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.isofox.TestUtils.CHR_1;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_1;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_2;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_3;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_5;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_NAME_2;
import static com.hartwig.hmftools.isofox.TestUtils.NEG_STRAND;
import static com.hartwig.hmftools.isofox.TestUtils.POS_STRAND;
import static com.hartwig.hmftools.isofox.TestUtils.addRacReadGroup;
import static com.hartwig.hmftools.isofox.TestUtils.addTestGenes;
import static com.hartwig.hmftools.isofox.TestUtils.createFusionFinder;
import static com.hartwig.hmftools.isofox.TestUtils.createGeneCollection;
import static com.hartwig.hmftools.isofox.TestUtils.addTestTranscripts;
import static com.hartwig.hmftools.isofox.TestUtils.createCigar;
import static com.hartwig.hmftools.isofox.TestUtils.createMappedRead;
import static com.hartwig.hmftools.isofox.TestUtils.createReadPair;
import static com.hartwig.hmftools.isofox.TestUtils.createSupplementaryReadPair;
import static com.hartwig.hmftools.isofox.TestUtils.overrideRefGenome;
import static com.hartwig.hmftools.isofox.TestUtils.populateRefGenome;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT_JUNCTION;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.MATCHED_JUNCTION;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGNED;
import static com.hartwig.hmftools.isofox.fusion.FusionTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.isofox.fusion.FusionTestUtils.createGroup;

import static htsjdk.samtools.SAMFlag.FIRST_OF_PAIR;
import static htsjdk.samtools.SAMFlag.SECOND_OF_PAIR;
import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.BaseDepth;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.ReadRecord;

import org.junit.Test;

public class FusionDataTest
{
    @Test
    public void testInvalidFusions()
    {
        final EnsemblDataCache geneTransCache = createGeneDataCache();

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        IsofoxConfig config = new IsofoxConfig();

        FusionFinder finder = createFusionFinder(config, geneTransCache);

        int gcId = 0;

        final GeneCollection gc1 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));

        // a DEL within the same gene collection is invalid
        ReadRecord read1 = createMappedRead(1, gc1, 1081, 1100, createCigar(0, 20, 20));
        ReadRecord read2 = createMappedRead(1, gc1, 1200, 1219, createCigar(20, 20, 0));

        final List<FusionReadGroup> chimericReadGroups = Lists.newArrayList();
        chimericReadGroups.add(createGroup(read1, read2));
        finder.processLocalReadGroups(chimericReadGroups);

        assertTrue(finder.getFusionCandidates().isEmpty());
    }

    @Test
    public void testFusionFragmentAssignment()
    {
        // Configurator.setRootLevel(Level.DEBUG);

        final EnsemblDataCache geneTransCache = createGeneDataCache();

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        IsofoxConfig config = new IsofoxConfig();
        FusionFinder finder = createFusionFinder(config, geneTransCache);

        int gcId = 0;

        final GeneCollection gc1 =
                createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));
        final GeneCollection gc2 =
                createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_2)));

        // 2 spliced fragments
        int readId = 0;
        ReadRecord read1 = createMappedRead(readId, gc1, 1050, 1089, createCigar(0, 40, 0));

        ReadRecord[] readPair = createSupplementaryReadPair(readId, gc1, gc2, 1081, 1100, 10200, 10219,
                createCigar(0, 20, 20), createCigar(20, 20, 0), true);

        readPair[1].setStrand(true, false);

        final Map<String, FusionReadGroup> readGroups1 = Maps.newHashMap();
        final Map<String, FusionReadGroup> readGroups2 = Maps.newHashMap();

        readGroups1.put(read1.Id, createGroup(read1, readPair[0]));
        readGroups2.put(read1.Id, createGroup(readPair[1]));

        readPair = createSupplementaryReadPair(++readId, gc1, gc2, 1081, 1100, 10200, 10219,
                createCigar(0, 20, 20), createCigar(20, 20, 0), true);

        ReadRecord read3 = createMappedRead(readId, gc2, 10210, 10249, createCigar(0, 40, 0));

        readPair[1].setStrand(true, false);
        read3.setStrand(true, false);

        readGroups1.put(read3.Id, createGroup(readPair[0]));
        readGroups2.put(read3.Id, createGroup(readPair[1], read3));

        // 1 unspliced fragment - will remain its own fusion
        read1 = createMappedRead(++readId, gc1, 1110, 1149, createCigar(0, 40, 0));

        readPair = createSupplementaryReadPair(readId, gc1, gc2, 1131, 1150, 10150, 10169,
                createCigar(0, 20, 20), createCigar(20, 20, 0), true);

        readPair[1].setStrand(true, false);

        readGroups1.put(read1.Id, createGroup(read1, readPair[0]));
        readGroups2.put(read1.Id, createGroup(readPair[1]));

        // 1 discordant fragment supporting the spliced fusion
        read1 = createMappedRead(++readId, gc1, 1055, 1084, createCigar(0, 40, 0));
        ReadRecord read2 = createMappedRead(readId, gc2, 10220, 10259, createCigar(0, 40, 0));
        read2.setStrand(true, false);

        readGroups1.put(read1.Id, createGroup(read1));
        readGroups2.put(read1.Id, createGroup(read2));

        // and 1 discordant fragment supporting the unspliced fusion
        read1 = createMappedRead(++readId, gc1, 1110, 1149, createCigar(0, 40, 0));
        read2 = createMappedRead(readId, gc2, 10160, 10199, createCigar(0, 40, 0));
        read2.setStrand(true, false);

        readGroups1.put(read1.Id, createGroup(read1));
        readGroups2.put(read1.Id, createGroup(read2));

        BaseDepth baseDepth = new BaseDepth();
        List<FusionReadGroup> completeGroups = finder.processNewChimericReadGroups(gc1, baseDepth, readGroups1);
        assertTrue(completeGroups.isEmpty());
        assertEquals(5, finder.getChimericPartialReadGroups().size());
        finder.processLocalReadGroups(completeGroups);

        // and the second GC's reads
        completeGroups = finder.processNewChimericReadGroups(gc2, baseDepth, readGroups2);
        assertTrue(finder.getChimericPartialReadGroups().isEmpty());
        finder.processLocalReadGroups(completeGroups);

        assertEquals(1, finder.getFusionCandidates().size());
        List<FusionReadData> fusions = finder.getFusionCandidates().values().iterator().next();
        assertEquals(2, fusions.size());

        FusionReadData fusion = fusions.stream().filter(x -> x.isKnownSpliced()).findFirst().orElse(null);
        assertTrue(fusion != null);
        assertEquals(2, fusion.getFragments(MATCHED_JUNCTION).size());
        assertEquals(1, fusion.getFragments(DISCORDANT).size());
        assertEquals(20, fusion.getMaxSplitLengths()[SE_START]);
        assertEquals(20, fusion.getMaxSplitLengths()[SE_END]);

        fusion = fusions.stream().filter(x -> x.isUnspliced()).findFirst().orElse(null);
        assertTrue(fusion != null);
        assertEquals(1, fusion.getFragments(MATCHED_JUNCTION).stream().filter(x -> x.isUnspliced()).count());
        assertEquals(2, fusion.getFragments(DISCORDANT).size());

        // check again with the discordant read having to fall within the correct transcript & exon
        finder.clearState();
        readGroups1.clear();
        readGroups2.clear();

        read1 = createMappedRead(++readId, gc1, 1050, 1089, createCigar(0, 40, 0));

        readPair = createSupplementaryReadPair(readId, gc1, gc2, 1081, 1100, 10200, 10219,
                createCigar(0, 20, 20), createCigar(20, 20, 0), true);

        readPair[1].setStrand(true, false);

        readGroups1.put(read1.Id, createGroup(read1, readPair[0]));
        readGroups2.put(read1.Id, createGroup(readPair[1]));

        // 1 discordant read in a valid location
        read1 = createMappedRead(++readId, gc1, 1020, 1059, createCigar(0, 40, 0));
        read2 = createMappedRead(readId, gc2, 10220, 10259, createCigar(0, 40, 0));
        read2.setStrand(true, false);

        readGroups1.put(read1.Id, createGroup(read1));
        readGroups2.put(read1.Id, createGroup(read2));

        // and another too many exons away
        read1 = createMappedRead(++readId, gc1, 1020, 1059, createCigar(0, 40, 0));
        read2 = createMappedRead(readId, gc2, 10820, 10859, createCigar(0, 40, 0));
        read2.setStrand(true, false);

        readGroups1.put(read1.Id, createGroup(read1));
        readGroups2.put(read1.Id, createGroup(read2));

        // and another too far away
        config.MaxFragmentLength = 200;
        read1 = createMappedRead(++readId, gc1, 1020, 1059, createCigar(0, 40, 0));
        read2 = createMappedRead(readId, gc2, 10720, 10759, createCigar(0, 40, 0));
        read2.setStrand(true, false);

        readGroups1.put(read1.Id, createGroup(read1));
        readGroups2.put(read1.Id, createGroup(read2));

        completeGroups = finder.processNewChimericReadGroups(gc1, baseDepth, readGroups1);
        assertEquals(4, finder.getChimericPartialReadGroups().size());
        finder.processLocalReadGroups(completeGroups);

        // and the second GC's reads
        completeGroups = finder.processNewChimericReadGroups(gc2, baseDepth, readGroups2);
        assertTrue(finder.getChimericPartialReadGroups().isEmpty());
        finder.processLocalReadGroups(completeGroups);

        fusions = finder.getFusionCandidates().values().iterator().next();
        assertEquals(1, fusions.size());
        fusion = fusions.get(0);
        assertEquals(1, fusion.getFragments(MATCHED_JUNCTION).size());
        assertEquals(1, fusion.getFragments(DISCORDANT).size());
        assertEquals(2, finder.getUnfusedFragments().get(fusion.locationId()).size());
    }

    @Test
    public void testInterChromosomalFusion()
    {
        // test with 2 -ve strand genes across chromosomes

        final EnsemblDataCache geneTransCache = createGeneDataCache();

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        IsofoxConfig config = new IsofoxConfig();
        populateRefGenome(config.RefGenome);

        int gcId = 0;

        final GeneCollection gc3 =
                createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_3)));
        final GeneCollection gc5 =
                createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_5)));

        FusionTaskManager fusionTaskManager = new FusionTaskManager(config, geneTransCache);
        FusionFinder finder = fusionTaskManager.createFusionFinder(gc3.chromosome());
        RacFragmentCache racFragmentCache = fusionTaskManager.racFragmentCache();

        // 2 spliced fragments
        int readId = 0;
        final Map<String,FusionReadGroup> readGroups1 = Maps.newHashMap();
        final Map<String,FusionReadGroup> readGroups2 = Maps.newHashMap();

        ReadRecord read1 = createMappedRead(readId, gc5, 10210, 10249, createCigar(0, 40, 20));

        ReadRecord[] readPair = createSupplementaryReadPair(readId, gc5, gc3, 10200, 10219, 20281, 20300,
                createCigar(20, 20, 0), createCigar(0, 20, 20), true);

        readPair[0].setStrand(true, false);

        readGroups1.put(read1.Id, createGroup(read1, readPair[0]));
        readGroups2.put(read1.Id, createGroup(readPair[1]));

        // RAC fragment for GC3
        String junctionBases = config.RefGenome.getBaseString(gc3.chromosome(), 20264, 20300)
                + config.RefGenome.getBaseString(gc5.chromosome(), 10200, 10202);

        read1 = createMappedRead(++readId, gc3, 20264, 20300, createCigar(0, 37, 3), junctionBases);
        ReadRecord read2 = createMappedRead(readId, gc3, 20250, 20289, createCigar(0, 40, 0));
        read2.setStrand(true, false);

        addRacReadGroup(racFragmentCache, new ChimericReadGroup(read1, read2), POS_ORIENT, 20300);

        // RAC fragment for GC5
        junctionBases = config.RefGenome.getBaseString(gc5.chromosome(), 20298, 20300)
                + config.RefGenome.getBaseString(gc3.chromosome(), 10200, 10236);

        read1 = createMappedRead(++readId, gc5, 10200, 10236, createCigar(3, 37, 0), junctionBases);
        read2 = createMappedRead(readId, gc5, 10210, 10259, createCigar(0, 40, 0));
        read2.setStrand(true, false);

        addRacReadGroup(racFragmentCache, new ChimericReadGroup(read1, read2), NEG_ORIENT, 10200);

        // 1 intronic discordant read
        ReadRecord[] discordantReads = createReadPair(++readId, gc3, gc5, 20150, 20189, 10320, 10359,
                createCigar(0, 40, 0),  createCigar(0, 40, 0), POS_STRAND, NEG_STRAND);

        readGroups1.put(discordantReads[0].Id, createGroup(discordantReads[1]));
        readGroups2.put(discordantReads[0].Id, createGroup(discordantReads[0]));

        // 1 exonic discordant read
        discordantReads = createReadPair(++readId, gc3, gc5, 20250, 20289, 10210, 10249,
                createCigar(0, 40, 0),  createCigar(0, 40, 0), POS_STRAND, NEG_STRAND);

        readGroups1.put(discordantReads[0].Id, createGroup(discordantReads[1]));
        readGroups2.put(discordantReads[0].Id, createGroup(discordantReads[0]));

        config.MaxFragmentLength = 500;

        // test our inter-chromosomal fusion calling

        BaseDepth baseDepth = new BaseDepth();
        List<FusionReadGroup> completeGroups = finder.processNewChimericReadGroups(gc3, baseDepth, readGroups1);
        assertEquals(3, finder.getChimericPartialReadGroups().size());
        finder.processLocalReadGroups(completeGroups);

        List<FusionReadGroup> interChromosomalGroups = fusionTaskManager.addIncompleteReadGroup(
                gc3.chromosome(), finder.extractIncompleteReadGroups(gc3.chromosome()));

        assertTrue(interChromosomalGroups.isEmpty());

        // and the second GC's reads (using the same FusionFinder is ok)
        FusionFinder finder2 = fusionTaskManager.createFusionFinder(gc5.chromosome());

        completeGroups = finder2.processNewChimericReadGroups(gc5, baseDepth, readGroups2);
        assertEquals(3, finder2.getChimericPartialReadGroups().size());
        finder2.processLocalReadGroups(completeGroups);

        interChromosomalGroups = fusionTaskManager.addIncompleteReadGroup(
                gc5.chromosome(), finder2.extractIncompleteReadGroups(gc5.chromosome()));

        finder2.processInterChromosomalReadGroups(interChromosomalGroups);

        assertEquals(3, interChromosomalGroups.size());

        List<FusionReadData> fusions = finder2.getFusionCandidates().values().iterator().next();
        assertEquals(1, fusions.size());
        FusionReadData fusion = fusions.get(0);
        assertEquals(1, fusion.getFragments(MATCHED_JUNCTION).size());
        assertEquals(2, fusion.getFragments(DISCORDANT).size());
        assertEquals(2, fusion.getFragments(REALIGNED).size());
        assertEquals(20, fusion.getMaxSplitLengths()[SE_START]);
        assertEquals(20, fusion.getMaxSplitLengths()[SE_END]);
    }

    @Test
    public void testSoftClippedFragmentRealignment()
    {
        final EnsemblDataCache geneTransCache = createGeneDataCache();

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        IsofoxConfig config = new IsofoxConfig();
        populateRefGenome(config.RefGenome);

        FusionFinder finder = createFusionFinder(config, geneTransCache);
        RacFragmentCache racFragmentCache = finder.racFragmentCache();

        int gcId = 0;

        final GeneCollection gc1 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));
        final GeneCollection gc2 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_2)));

        // a spliced fragment to establish the fusion
        int readId = 0;
        ReadRecord read1 = createMappedRead(readId, gc1, 1050, 1089, createCigar(0, 40, 0));

        ReadRecord[] readPair = createSupplementaryReadPair(readId, gc1, gc2, 1081, 1100, 10200, 10219,
                createCigar(0, 20, 20), createCigar(20, 20, 0), true);

        readPair[1].setStrand(true, false);

        final Map<String, FusionReadGroup> readGroups1 = Maps.newHashMap();
        final Map<String, FusionReadGroup> readGroups2 = Maps.newHashMap();

        readGroups1.put(read1.Id, createGroup(read1, readPair[0]));
        readGroups2.put(read1.Id, createGroup(readPair[1]));

        // a soft-clipped read matching the other side of the fusion junction
        String junctionBases = config.RefGenome.getBaseString(gc1.chromosome(), 1071, 1100)
                + config.RefGenome.getBaseString(gc1.chromosome(), 10200, 10209);
        ReadRecord read4 = createMappedRead(++readId, gc1, 1071, 1100, createCigar(0, 30, 10), junctionBases);
        ReadRecord read5 = createMappedRead(readId, gc1, 1051, 1090, createCigar(0, 40, 0));
        read5.setStrand(true, false);

        addRacReadGroup(racFragmentCache, new ChimericReadGroup(read4, read4), POS_ORIENT, 1100);

        // a soft-clipped read matching 2 bases into the ref due to homology with the other side of the fusion junction
        junctionBases = config.RefGenome.getBaseString(gc1.chromosome(), 1091, 1100)
                + config.RefGenome.getBaseString(gc1.chromosome(), 10200, 10229);

        ReadRecord read6 = createMappedRead(++readId, gc2, 10198, 10229, createCigar(8, 32, 0), junctionBases);
        ReadRecord read7 = createMappedRead(readId, gc2, 10210, 10249, createCigar(0, 40, 0));
        read7.setStrand(true, false);

        addRacReadGroup(racFragmentCache, new ChimericReadGroup(read6, read7), NEG_ORIENT, 10200);

        BaseDepth baseDepth = new BaseDepth();
        List<FusionReadGroup> completeGroups = finder.processNewChimericReadGroups(gc1, baseDepth, readGroups1);
        finder.processLocalReadGroups(completeGroups);

        completeGroups = finder.processNewChimericReadGroups(gc2, baseDepth, readGroups2);
        finder.processLocalReadGroups(completeGroups);

        assertEquals(1, finder.getFusionCandidates().size());
        List<FusionReadData> fusions = finder.getFusionCandidates().values().iterator().next();
        assertEquals(1, fusions.size());

        FusionReadData fusion = fusions.stream().filter(x -> x.isKnownSpliced()).findFirst().orElse(null);
        assertTrue(fusion != null);
        assertEquals(1, fusion.getFragments(MATCHED_JUNCTION).size());
        assertEquals(2, fusion.getFragments(REALIGNED).size());
    }

    @Test
    public void testLocalDelFusion()
    {
        // split reads at known junctions but too short to have supplementary data
        final EnsemblDataCache geneTransCache = createGeneDataCache();

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        IsofoxConfig config = new IsofoxConfig();
        populateRefGenome(config.RefGenome);

        FusionFinder finder = createFusionFinder(config, geneTransCache);
        RacFragmentCache racFragmentCache = finder.racFragmentCache();

        config.Fusions.KnownFusions.addData(new KnownFusionData(KNOWN_PAIR, GENE_NAME_1, GENE_NAME_2, "", ""));

        int gcId = 0;

        final GeneCollection gc1 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));
        final GeneCollection gc2 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_2)));

        final Map<String, FusionReadGroup> readGroups1 = Maps.newHashMap();
        final Map<String, FusionReadGroup> readGroups2 = Maps.newHashMap();

        // a simple DEL, supported by 2 split fragments - the discordant read forming the known-pair fusion
        int readId = 0;
        ReadRecord read1 = createMappedRead(readId, gc1, 1081, 1100, createCigar(0, 20, 20));
        ReadRecord read2 = createMappedRead(readId, gc2, 10200, 10219, createCigar(20, 20, 0));
        read2.setStrand(true, false);

        readGroups1.put(read1.Id, createGroup(read1));
        readGroups2.put(read1.Id, createGroup(read2));

        // a second one
        read1 = createMappedRead(++readId, gc1, 1081, 1100, createCigar(0, 20, 20));
        read2 = createMappedRead(readId, gc2, 10200, 10219, createCigar(20, 20, 0));
        read2.setStrand(true, false);

        readGroups1.put(read1.Id, createGroup(read1));
        readGroups2.put(read1.Id, createGroup(read2));

        // realigned and discordant reads which support it

        // 2 soft-clipped fragments matching the other side of the fusion junction
        String junctionBases = config.RefGenome.getBaseString(gc1.chromosome(), 1071, 1100)
                + config.RefGenome.getBaseString(gc2.chromosome(), 10200, 10209);

        ReadRecord read3 = createMappedRead(++readId, gc1, 1071, 1100, createCigar(0, 30, 10), junctionBases);
        ReadRecord read4 = createMappedRead(readId, gc1, 1051, 1090, createCigar(0, 40, 0));
        read4.setStrand(true, false);

        // readGroups1.put(read3.Id, new ReadGroup(read3, read4));

        addRacReadGroup(racFragmentCache, new ChimericReadGroup(read3, read4), POS_ORIENT, 1100);

        junctionBases = config.RefGenome.getBaseString(gc1.chromosome(), 1091, 1100)
                + config.RefGenome.getBaseString(gc2.chromosome(), 10200, 10229);

        ReadRecord read5 = createMappedRead(++readId, gc2, 10200, 10229, createCigar(10, 30, 0), junctionBases);
        ReadRecord read6 = createMappedRead(readId, gc2, 10210, 10249, createCigar(0, 40, 0));
        read6.setStrand(true, false);

        readGroups1.put(read5.Id, createGroup(read5, read6));

        addRacReadGroup(racFragmentCache, new ChimericReadGroup(read5, read6), NEG_ORIENT, 10200);

        // and a discordant fragment
        read3 = createMappedRead(++readId, gc1, 1050, 1089, createCigar(0, 40, 0));
        read4 = createMappedRead(readId, gc2, 10210, 10249, createCigar(0, 40, 0));
        read4.setStrand(true, false);

        readGroups1.put(read3.Id, createGroup(read3));
        readGroups2.put(read4.Id, createGroup(read4));

        BaseDepth baseDepth = new BaseDepth();
        List<FusionReadGroup> completeGroups = finder.processNewChimericReadGroups(gc1, baseDepth, readGroups1);
        finder.processLocalReadGroups(completeGroups);

        completeGroups = finder.processNewChimericReadGroups(gc2, baseDepth, readGroups2);
        finder.processLocalReadGroups(completeGroups);

        assertEquals(1, finder.getFusionCandidates().size());
        List<FusionReadData> fusions = finder.getFusionCandidates().values().iterator().next();
        assertEquals(1, fusions.size());

        FusionReadData fusion = fusions.stream().filter(x -> x.isKnownSpliced()).findFirst().orElse(null);
        assertTrue(fusion != null);
        assertEquals(DISCORDANT_JUNCTION, fusion.getInitialFragment().type());
        assertEquals(2, fusion.getFragmentTypeCount(DISCORDANT_JUNCTION));
        assertEquals(1, fusion.getFragmentTypeCount(DISCORDANT));
        assertEquals(2, fusion.getFragments(REALIGNED).size());
        assertEquals(20, fusion.getMaxSplitLengths()[SE_START]);
        assertEquals(20, fusion.getMaxSplitLengths()[SE_END]);
    }

    @Test
    public void testLocalSplitReadDelFusions()
    {
        // a split read (DEL) spanning 2 genes
        final EnsemblDataCache geneTransCache = createGeneDataCache();

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        IsofoxConfig config = new IsofoxConfig();
        populateRefGenome(config.RefGenome);

        FusionFinder finder = createFusionFinder(config, geneTransCache);
        RacFragmentCache racFragmentCache = finder.racFragmentCache();

        int gcId = 0;

        final GeneCollection gc1 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));
        final GeneCollection gc2 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_2)));

        // handle reads from GC 1 and then 2 as the BamFragmentReader would

        int readId = 0;
        ReadRecord read1 = createMappedRead(readId, gc1, 1081, 10219, createCigar(0, 20, 9099, 20, 0));
        read1.setGeneCollection(SE_START, gc1.id(), true);

        ReadRecord read2 = createMappedRead(readId, gc2, 10210, 10249, createCigar(0, 40, 0));
        read2.setStrand(true, false);

        // realigned and discordant reads which support it

        // 2 soft-clipped fragments matching the other side of the fusion junction
        String junctionBases = config.RefGenome.getBaseString(gc1.chromosome(), 1071, 1100)
                + config.RefGenome.getBaseString(gc2.chromosome(), 10200, 10209);

        // first on GC1
        ReadRecord read3 = createMappedRead(++readId, gc1, 1071, 1100, createCigar(0, 30, 10), junctionBases);
        ReadRecord read4 = createMappedRead(readId, gc1, 1051, 1090, createCigar(0, 40, 0));
        read4.setStrand(true, false);

        addRacReadGroup(racFragmentCache, new ChimericReadGroup(read3, read4), POS_ORIENT, 1100);

        // then on GC2
        junctionBases = config.RefGenome.getBaseString(gc1.chromosome(), 1091, 1100)
                + config.RefGenome.getBaseString(gc2.chromosome(), 10200, 10229);
        ReadRecord read5 = createMappedRead(++readId, gc2, 10200, 10229, createCigar(10, 30, 0), junctionBases);
        ReadRecord read6 = createMappedRead(readId, gc2, 10210, 10249, createCigar(0, 40, 0));
        read6.setStrand(true, false);

        addRacReadGroup(racFragmentCache, new ChimericReadGroup(read5, read6), NEG_ORIENT, 10200);

        // and a discordant fragment
        ReadRecord read7 = createMappedRead(++readId, gc1, 1050, 1089, createCigar(0, 40, 0));
        ReadRecord read8 = createMappedRead(readId, gc2, 10210, 10249, createCigar(0, 40, 0));
        read8.setStrand(true, false);

        // another split read with both reads in the first GC
        ReadRecord read9 = createMappedRead(++readId, gc1, 1081, 10219, createCigar(0, 20, 9099, 20, 0));
        ReadRecord read10 = createMappedRead(readId, gc1, 1051, 1090, createCigar(0, 40, 0));
        read10.setStrand(true, false);

        final Map<String, FusionReadGroup> chimericReadGroups = Maps.newHashMap();

        chimericReadGroups.put(read1.Id, createGroup(read1));
        chimericReadGroups.put(read7.Id, createGroup(read7));
        chimericReadGroups.put(read9.Id, createGroup(read9, read10));

        BaseDepth baseDepth = new BaseDepth();

        List<FusionReadGroup> completeGroups = finder.processNewChimericReadGroups(gc1, baseDepth, chimericReadGroups);
        assertEquals(2, finder.getSpanningReadGroups().size());
        assertTrue(completeGroups.isEmpty());
        finder.processLocalReadGroups(completeGroups); // will result in an unassigned RA fragment from reads 3 & 4

        // GC 2 read handling
        chimericReadGroups.clear();
        chimericReadGroups.put(read2.Id, createGroup(read2));
        chimericReadGroups.put(read8.Id, createGroup(read8));

        completeGroups = finder.processNewChimericReadGroups(gc2, baseDepth, chimericReadGroups);
        assertEquals(0, finder.getSpanningReadGroups().size());
        assertEquals(3, completeGroups.size());
        finder.processLocalReadGroups(completeGroups);

        assertEquals(1, finder.getFusionCandidates().size());
        List<FusionReadData> fusions = finder.getFusionCandidates().values().iterator().next();
        assertEquals(1, fusions.size());

        FusionReadData fusion = fusions.stream().filter(x -> x.isKnownSpliced()).findFirst().orElse(null);
        assertTrue(fusion != null);
        assertEquals(2, fusion.getFragments(MATCHED_JUNCTION).size());
        assertEquals(2, fusion.getFragments(REALIGNED).size());
        assertEquals(1, fusion.getFragments(DISCORDANT).size());
        assertEquals(20, fusion.getMaxSplitLengths()[SE_START]);
        assertEquals(20, fusion.getMaxSplitLengths()[SE_END]);
    }

    @Test
    public void testNonGenicFusions()
    {
        // fragments spanning from pre-gene to the following gene
        final EnsemblDataCache geneTransCache = createGeneDataCache();

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        IsofoxConfig config = new IsofoxConfig();
        populateRefGenome(config.RefGenome);

        FusionFinder finder = createFusionFinder(config, geneTransCache);
        RacFragmentCache racFragmentCache = finder.racFragmentCache();

        int gcId = 0;

        final GeneCollection gc1 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));
        final GeneCollection gc2 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_2)));

        final Map<String, FusionReadGroup> readGroups1 = Maps.newHashMap();
        final Map<String, FusionReadGroup> readGroups2 = Maps.newHashMap();

        // a simple DEL, supported by a split-read
        int readId = 0;
        ReadRecord read1 = createMappedRead(readId, gc1, 581, 10219, createCigar(0, 20, 9599, 20, 0));
        ReadRecord read2 = createMappedRead(readId, gc2, 10205, 10224, createCigar(0, 20, 0));
        read2.setStrand(true, false);

        readGroups1.put(read1.Id, createGroup(read1));
        readGroups2.put(read2.Id, createGroup(read2));

        // realigned and discordant reads which support it

        // 2 soft-clipped fragments matching the other side of the fusion junction
        String junctionBases = config.RefGenome.getBaseString(gc1.chromosome(), 571, 600)
                + config.RefGenome.getBaseString(gc2.chromosome(), 10200, 10209);

        ReadRecord read3 = createMappedRead(++readId, gc1, 571, 600, createCigar(0, 30, 10), junctionBases);
        ReadRecord read4 = createMappedRead(readId, gc1, 551, 590, createCigar(0, 40, 0));
        read4.setStrand(true, false);

        // readGroups1.put(read3.Id, new ReadGroup(read3, read4));
        addRacReadGroup(racFragmentCache, new ChimericReadGroup(read3, read4), POS_ORIENT, 600);

        junctionBases = config.RefGenome.getBaseString(gc1.chromosome(), 591, 600)
                + config.RefGenome.getBaseString(gc2.chromosome(), 10200, 10229);

        ReadRecord read5 = createMappedRead(++readId, gc2, 10200, 10229, createCigar(10, 30, 0), junctionBases);
        ReadRecord read6 = createMappedRead(readId, gc2, 10210, 10249, createCigar(0, 40, 0));
        read6.setStrand(true, false);

        readGroups1.put(read5.Id, createGroup(read5 ,read6));
        addRacReadGroup(racFragmentCache, new ChimericReadGroup(read5, read6), NEG_ORIENT, 10200);

        // and a discordant fragment
        read3 = createMappedRead(++readId, gc1, 550, 589, createCigar(0, 40, 0));
        read4 = createMappedRead(readId, gc2, 10210, 10249, createCigar(0, 40, 0));
        read4.setStrand(true, false);

        readGroups1.put(read3.Id, createGroup(read3));
        readGroups2.put(read4.Id, createGroup(read4));

        BaseDepth baseDepth = new BaseDepth();
        List<FusionReadGroup> completeGroups = finder.processNewChimericReadGroups(gc1, baseDepth, readGroups1);
        finder.processLocalReadGroups(completeGroups); // will result in an unassigned RA fragment from reads 3 & 4

        completeGroups = finder.processNewChimericReadGroups(gc2, baseDepth, readGroups2);
        finder.processLocalReadGroups(completeGroups); // will result in an unassigned RA fragment from reads 3 & 4

        assertEquals(1, finder.getFusionCandidates().size());
        List<FusionReadData> fusions = finder.getFusionCandidates().values().iterator().next();
        assertEquals(1, fusions.size());

        FusionReadData fusion = fusions.stream().findFirst().orElse(null);
        assertTrue(fusion != null);
        assertEquals(MATCHED_JUNCTION, fusion.getInitialFragment().type());
        assertEquals(1, fusion.getFragmentTypeCount(MATCHED_JUNCTION));
        assertEquals(2, fusion.getFragmentTypeCount(REALIGNED));
    }

    @Test
    public void testHomologyMerging()
    {
        // Configurator.setRootLevel(Level.DEBUG);

        final EnsemblDataCache geneTransCache = createGeneDataCache();

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        IsofoxConfig config = new IsofoxConfig();
        populateRefGenome(config.RefGenome);

        FusionFinder finder = createFusionFinder(config, geneTransCache);

        int gcId = 0;
        final GeneCollection gc1 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));
        final GeneCollection gc2 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_2)));

        final Map<String, FusionReadGroup> readGroups1 = Maps.newHashMap();
        final Map<String, FusionReadGroup> readGroups2 = Maps.newHashMap();

        // a spliced fragment to establish the fusion
        int readId = 0;
        ReadRecord read1 = createMappedRead(readId, gc1, 1050, 1089, createCigar(0, 40, 0));

        // create 3 bases of homology around the junction
        overrideRefGenome(config.RefGenome, CHR_1, 1101, "AGT");
        overrideRefGenome(config.RefGenome, CHR_1, 10200, "AGT");

        String junctionBases = config.RefGenome.getBaseString(gc1.chromosome(), 1071, 1100)
                + config.RefGenome.getBaseString(gc1.chromosome(), 10200, 10209);

        ReadRecord read2 = createMappedRead(readId, gc1, 1071, 1100, createCigar(0, 30, 10), junctionBases);

        junctionBases = config.RefGenome.getBaseString(gc1.chromosome(), 1091, 1100)
                + config.RefGenome.getBaseString(gc1.chromosome(), 10200, 10229);

        ReadRecord read3 = createMappedRead(readId, gc2, 10200, 102929, createCigar(10, 30, 0), junctionBases);

        read2.setFlag(FIRST_OF_PAIR, true);
        read3.setFlag(SECOND_OF_PAIR, false);
        read3.setStrand(true, false);
        read2.setSuppAlignment(String.format("%s;%d;%s", read3.Chromosome, read3.PosStart, read3.Cigar.toString()));
        read3.setSuppAlignment(String.format("%s;%d;%s", read2.Chromosome, read2.PosStart, read2.Cigar.toString()));

        readGroups1.put(read1.Id, createGroup(read1, read2));
        readGroups2.put(read1.Id, createGroup(read3));

        // a second fragment with 2 bases difference due to homology
        read1 = createMappedRead(++readId, gc1, 1050, 1089, createCigar(0, 40, 0));

        junctionBases = config.RefGenome.getBaseString(gc1.chromosome(), 1073, 1102)
                + config.RefGenome.getBaseString(gc1.chromosome(), 10202, 10211);

        read2 = createMappedRead(readId, gc1, 1073, 1102, createCigar(0, 30, 10), junctionBases);

        junctionBases = config.RefGenome.getBaseString(gc1.chromosome(), 1093, 1102)
                + config.RefGenome.getBaseString(gc1.chromosome(), 10202, 10231);

        read3 = createMappedRead(readId, gc2, 10202, 102931, createCigar(10, 30, 0), junctionBases);

        read2.setFlag(FIRST_OF_PAIR, true);
        read3.setFlag(SECOND_OF_PAIR, false);
        read3.setStrand(true, false);
        read2.setSuppAlignment(String.format("%s;%d;%s", read3.Chromosome, read3.PosStart, read3.Cigar.toString()));
        read3.setSuppAlignment(String.format("%s;%d;%s", read2.Chromosome, read2.PosStart, read2.Cigar.toString()));

        readGroups1.put(read1.Id, createGroup(read1, read2));
        readGroups2.put(read1.Id, createGroup(read3));

        BaseDepth baseDepth = new BaseDepth();
        List<FusionReadGroup> completeGroups = finder.processNewChimericReadGroups(gc1, baseDepth, readGroups1);
        finder.processLocalReadGroups(completeGroups); // will result in an unassigned RA fragment from reads 3 & 4

        completeGroups = finder.processNewChimericReadGroups(gc2, baseDepth, readGroups2);
        finder.processLocalReadGroups(completeGroups); // will result in an unassigned RA fragment from reads 3 & 4

        assertEquals(1, finder.getFusionCandidates().size());
        List<FusionReadData> fusions = finder.getFusionCandidates().values().iterator().next();
        assertEquals(1, fusions.size());

        FusionReadData fusion = fusions.get(0);
        assertTrue(fusion != null);
        assertEquals(2, fusion.getFragments(MATCHED_JUNCTION).size());
    }

    @Test
    public void testCloseMatchFiltering()
    {
        final EnsemblDataCache geneTransCache = createGeneDataCache();

        IsofoxConfig config = new IsofoxConfig();

        FusionFinder finder = createFusionFinder(config, geneTransCache);

        int gcId = 0;

        final GeneCollection gc1 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));
        final GeneCollection gc2 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_2)));

        final Map<String, FusionReadGroup> readGroups1 = Maps.newHashMap();
        final Map<String, FusionReadGroup> readGroups2 = Maps.newHashMap();

        // 3 different but close spliced fragments, one of them with more support than the others
        int readId = 0;
        ReadRecord read1 = createMappedRead(readId, gc1, 1050, 1089, createCigar(0, 40, 0));

        ReadRecord[] readPair = createSupplementaryReadPair(readId, gc1, gc2, 1081, 1100, 10200, 10219,
                createCigar(0, 20, 20), createCigar(20, 20, 0), true);

        readPair[1].setStrand(true, false);

        readGroups1.put(read1.Id, createGroup(read1, readPair[0]));
        readGroups2.put(read1.Id, createGroup(readPair[1]));

        // more support - needs to be 5x or more
        for(int i = 0; i < 4; ++i)
        {
            read1 = createMappedRead(++readId, gc1, 1050, 1089, createCigar(0, 40, 0));

            readPair = createSupplementaryReadPair(readId, gc1, gc2, 1081, 1100, 10200, 10219,
                    createCigar(0, 20, 20), createCigar(20, 20, 0), true);

            readPair[1].setStrand(true, false);

            readGroups1.put(read1.Id, createGroup(read1, readPair[0]));
            readGroups2.put(read1.Id, createGroup(readPair[1]));
        }

        // another close by
        read1 = createMappedRead(++readId, gc1, 1050, 1089, createCigar(0, 40, 0));

        readPair = createSupplementaryReadPair(readId, gc1, gc2, 1080, 1099, 10199, 10218,
                createCigar(0, 20, 20), createCigar(20, 20, 0), true);

        readPair[1].setStrand(true, false);

        readGroups1.put(read1.Id, createGroup(read1, readPair[0]));
        readGroups2.put(read1.Id, createGroup(readPair[1]));

        // and another
        read1 = createMappedRead(++readId, gc1, 1050, 1089, createCigar(0, 40, 0));

        readPair = createSupplementaryReadPair(readId, gc1, gc2, 1082, 1101, 10201, 10220,
                createCigar(0, 20, 20), createCigar(20, 20, 0), true);

        readPair[1].setStrand(true, false);

        readGroups1.put(read1.Id, createGroup(read1, readPair[0]));
        readGroups2.put(read1.Id, createGroup(readPair[1]));

        BaseDepth baseDepth = new BaseDepth();
        List<FusionReadGroup> completeGroups = finder.processNewChimericReadGroups(gc1, baseDepth, readGroups1);
        finder.processLocalReadGroups(completeGroups); // will result in an unassigned RA fragment from reads 3 & 4

        completeGroups = finder.processNewChimericReadGroups(gc2, baseDepth, readGroups2);
        finder.processLocalReadGroups(completeGroups); // will result in an unassigned RA fragment from reads 3 & 4

        assertEquals(1, finder.getFusionCandidates().size());
        List<FusionReadData> fusions = finder.getFusionCandidates().values().iterator().next();
        assertEquals(1, fusions.size());

        FusionReadData fusion = fusions.get(0);
        assertTrue(fusion != null);
        assertEquals(5, fusion.getFragments(MATCHED_JUNCTION).size());
        assertEquals(0, fusion.getFragments(REALIGNED).size());
        assertEquals(1100, fusion.junctionPositions()[SE_START]);
        assertEquals(10200, fusion.junctionPositions()[SE_END]);
    }
}

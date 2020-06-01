package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createGeneDataCache;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.isofox.TestUtils.CHR_1;
import static com.hartwig.hmftools.isofox.TestUtils.CHR_2;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_1;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_2;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_3;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_4;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_5;
import static com.hartwig.hmftools.isofox.TestUtils.NEG_STRAND;
import static com.hartwig.hmftools.isofox.TestUtils.POS_STRAND;
import static com.hartwig.hmftools.isofox.TestUtils.TRANS_1;
import static com.hartwig.hmftools.isofox.TestUtils.TRANS_2;
import static com.hartwig.hmftools.isofox.TestUtils.TRANS_3;
import static com.hartwig.hmftools.isofox.TestUtils.addTestGenes;
import static com.hartwig.hmftools.isofox.TestUtils.createGeneCollection;
import static com.hartwig.hmftools.isofox.TestUtils.addTestTranscripts;
import static com.hartwig.hmftools.isofox.TestUtils.createCigar;
import static com.hartwig.hmftools.isofox.TestUtils.createMappedRead;
import static com.hartwig.hmftools.isofox.TestUtils.createReadPair;
import static com.hartwig.hmftools.isofox.TestUtils.createSupplementaryReadPair;
import static com.hartwig.hmftools.isofox.TestUtils.generateRandomBases;
import static com.hartwig.hmftools.isofox.TestUtils.populateRefGenome;
import static com.hartwig.hmftools.isofox.common.TransExonRef.hasTranscriptExonMatch;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.MATCHED_JUNCTION;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.DISCORDANT;
import static com.hartwig.hmftools.isofox.fusion.FusionFragmentType.REALIGNED;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.TransExonRef;
import com.hartwig.hmftools.isofox.fusion.FusionFinder;
import com.hartwig.hmftools.isofox.fusion.FusionFragment;
import com.hartwig.hmftools.isofox.fusion.FusionReadData;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.junit.Assert;
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

        FusionFinder finder = new FusionFinder(config, geneTransCache);

        int gcId = 0;

        final GeneCollection gc1 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));

        Map<Integer,List<EnsemblGeneData>> gcMap = Maps.newHashMap();

        gcMap.put(gc1.id(), gc1.genes().stream().map(x -> x.GeneData).collect(Collectors.toList()));
        finder.addChromosomeGeneCollections(CHR_1, gcMap);

        // a DEL within the same gene collection is invalid
        ReadRecord read1 = createMappedRead(1, gc1, 1081, 1100, createCigar(0, 20, 20));
        ReadRecord read2 = createMappedRead(1, gc1, 1200, 1219, createCigar(20, 20, 0));

        final Map<String,List<ReadRecord>> chimericReadMap = Maps.newHashMap();
        chimericReadMap.put(read1.Id, Lists.newArrayList(read1, read2));
        finder.addChimericReads(chimericReadMap);

        finder.findFusions();

        assertTrue(finder.getFusionCandidates().isEmpty());
    }

    @Test
    public void testTransExonComparisons()
    {
        List<TransExonRef> transExons1 = Lists.newArrayList();
        List<TransExonRef> transExons2 = Lists.newArrayList();

        transExons1.add(new TransExonRef(GENE_ID_1, TRANS_1, "TRANS1", 1));
        transExons1.add(new TransExonRef(GENE_ID_2, TRANS_2, "TRANS2", 4));

        transExons2.add(new TransExonRef(GENE_ID_3, TRANS_3, "TRANS3", 3));
        transExons2.add(new TransExonRef(GENE_ID_1, TRANS_1, "TRANS1", 1));

        assertTrue(hasTranscriptExonMatch(transExons1, transExons2));

        transExons2.clear();

        assertFalse(hasTranscriptExonMatch(transExons1, transExons2));
        transExons2.add(new TransExonRef(GENE_ID_2, TRANS_2, "TRANS2", 2));

        assertFalse(hasTranscriptExonMatch(transExons1, transExons2));
        assertFalse(hasTranscriptExonMatch(transExons1, transExons2, -1));
        assertTrue(hasTranscriptExonMatch(transExons1, transExons2, -2));
        assertFalse(hasTranscriptExonMatch(transExons1, transExons2, 2));

        transExons2.clear();
        transExons2.add(new TransExonRef(GENE_ID_2, TRANS_2, "TRANS2", 6));

        assertFalse(hasTranscriptExonMatch(transExons1, transExons2));
        assertFalse(hasTranscriptExonMatch(transExons1, transExons2, 1));
        assertTrue(hasTranscriptExonMatch(transExons1, transExons2, 2));
        assertFalse(hasTranscriptExonMatch(transExons1, transExons2, -2));
    }

    @Test
    public void testFusionFragmentAssignment()
    {
        Configurator.setRootLevel(Level.DEBUG);

        final EnsemblDataCache geneTransCache = createGeneDataCache();

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        IsofoxConfig config = new IsofoxConfig();
        FusionFinder finder = new FusionFinder(config, geneTransCache);

        int gcId = 0;

        final GeneCollection gc1 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));
        final GeneCollection gc2 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_2)));
        final GeneCollection gc3 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_3)));
        final GeneCollection gc4 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_4)));
        final GeneCollection gc5 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_5)));

        Map<Integer,List<EnsemblGeneData>> gcMap = Maps.newHashMap();

        gcMap.put(gc1.id(), gc1.genes().stream().map(x -> x.GeneData).collect(Collectors.toList()));
        gcMap.put(gc2.id(), gc2.genes().stream().map(x -> x.GeneData).collect(Collectors.toList()));
        gcMap.put(gc3.id(), gc3.genes().stream().map(x -> x.GeneData).collect(Collectors.toList()));
        finder.addChromosomeGeneCollections(CHR_1, gcMap);

        gcMap = Maps.newHashMap();
        gcMap.put(gc4.id(), gc4.genes().stream().map(x -> x.GeneData).collect(Collectors.toList()));
        gcMap.put(gc5.id(), gc5.genes().stream().map(x -> x.GeneData).collect(Collectors.toList()));
        finder.addChromosomeGeneCollections(CHR_2, gcMap);

        final Map<String,List<ReadRecord>> chimericReadMap = Maps.newHashMap();

        // 2 spliced fragments
        int readId = 0;
        ReadRecord read1 = createMappedRead(readId, gc1, 1050, 1089, createCigar(0, 40, 0));

        ReadRecord[] readPair = createSupplementaryReadPair(readId, gc1, gc2, 1081, 1100, 10200, 10219,
                createCigar(0, 20, 20), createCigar(20, 20, 0), true);

        readPair[1].setStrand(true, false);

        chimericReadMap.put(read1.Id, Lists.newArrayList(read1, readPair[0], readPair[1]));

        readPair = createSupplementaryReadPair(++readId, gc1, gc2, 1081, 1100, 10200, 10219,
                createCigar(0, 20, 20), createCigar(20, 20, 0), true);

        ReadRecord read3 = createMappedRead(readId, gc2, 10210, 10249, createCigar(0, 40, 0));

        readPair[1].setStrand(true, false);
        read3.setStrand(true, false);

        chimericReadMap.put(read3.Id, Lists.newArrayList(readPair[0], readPair[1], read3));

        // 1 unspliced fragment - will remain its own fusion
        read1 = createMappedRead(++readId, gc1, 1110, 1149, createCigar(0, 40, 0));

        readPair = createSupplementaryReadPair(readId, gc1, gc2, 1131, 1150, 10150, 10169,
                createCigar(0, 20, 20), createCigar(20, 20, 0), true);

        readPair[1].setStrand(true, false);

        chimericReadMap.put(read1.Id, Lists.newArrayList(read1, readPair[0], readPair[1]));

        // 1 discordant fragment supporting the spliced fusion
        read1 = createMappedRead(++readId, gc1, 1055, 1084, createCigar(0, 40, 0));
        ReadRecord read2 = createMappedRead(readId, gc2, 10220, 10259, createCigar(0, 40, 0));
        read2.setStrand(true, false);
        chimericReadMap.put(read1.Id, Lists.newArrayList(read1, read2));

        // and 1 discordant fragment supporting the unspliced fusion
        read1 = createMappedRead(++readId, gc1, 1110, 1149, createCigar(0, 40, 0));
        read2 = createMappedRead(readId, gc2, 10160, 10199, createCigar(0, 40, 0));
        read2.setStrand(true, false);
        chimericReadMap.put(read1.Id, Lists.newArrayList(read1, read2));

        finder.addChimericReads(chimericReadMap);

        finder.findFusions();

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
        chimericReadMap.clear();

        read1 = createMappedRead(++readId, gc1, 1050, 1089, createCigar(0, 40, 0));

        readPair = createSupplementaryReadPair(readId, gc1, gc2, 1081, 1100, 10200, 10219,
                createCigar(0, 20, 20), createCigar(20, 20, 0), true);

        readPair[1].setStrand(true, false);
        chimericReadMap.put(read1.Id, Lists.newArrayList(read1, readPair[0], readPair[1]));

        // 1 discordant read in a valid location
        read1 = createMappedRead(++readId, gc1, 1020, 1059, createCigar(0, 40, 0));
        read2 = createMappedRead(readId, gc2, 10220, 10259, createCigar(0, 40, 0));
        read2.setStrand(true, false);
        chimericReadMap.put(read1.Id, Lists.newArrayList(read1, read2));

        // and another too many exons away
        read1 = createMappedRead(++readId, gc1, 1020, 1059, createCigar(0, 40, 0));
        read2 = createMappedRead(readId, gc2, 10820, 10859, createCigar(0, 40, 0));
        read2.setStrand(true, false);
        chimericReadMap.put(read1.Id, Lists.newArrayList(read1, read2));

        // and another too far away
        config.MaxFragmentLength = 200;
        read1 = createMappedRead(++readId, gc1, 1020, 1059, createCigar(0, 40, 0));
        read2 = createMappedRead(readId, gc2, 10720, 10759, createCigar(0, 40, 0));
        read2.setStrand(true, false);
        chimericReadMap.put(read1.Id, Lists.newArrayList(read1, read2));

        finder.addChimericReads(chimericReadMap);

        finder.findFusions();

        fusions = finder.getFusionCandidates().values().iterator().next();
        assertEquals(1, fusions.size());
        fusion = fusions.get(0);
        assertEquals(1, fusion.getFragments(MATCHED_JUNCTION).size());
        assertEquals(1, fusion.getFragments(DISCORDANT).size());
        assertEquals(2, finder.getUnfusedFragments().get(fusion.locationId()).size());

        // test again but with 2 -ve strand genes
        finder.clearState();
        chimericReadMap.clear();

        read2 = createMappedRead(readId, gc5, 10210, 10249, createCigar(0, 40, 20));

        readPair = createSupplementaryReadPair(readId, gc5, gc3, 10200, 10219, 20281, 20300,
                createCigar(20, 20, 0), createCigar(0, 20, 20), true);

        readPair[0].setStrand(true, false);

        chimericReadMap.put(read1.Id, Lists.newArrayList(readPair[0], read2, readPair[1]));

        // 1 intronic discordant read
        ReadRecord[] discordantReads = createReadPair(++readId, gc3, gc5, 20150, 20189, 10320, 10359,
                createCigar(0, 40, 0),  createCigar(0, 40, 0), POS_STRAND, NEG_STRAND);

        chimericReadMap.put(discordantReads[0].Id, Lists.newArrayList(discordantReads[0], discordantReads[1]));

        // 1 exonic discordant read
        discordantReads = createReadPair(++readId, gc3, gc5, 20250, 20289, 10210, 10249,
                createCigar(0, 40, 0),  createCigar(0, 40, 0), POS_STRAND, NEG_STRAND);

        chimericReadMap.put(discordantReads[0].Id, Lists.newArrayList(discordantReads[0], discordantReads[1]));

        finder.addChimericReads(chimericReadMap);

        config.MaxFragmentLength = 500;
        finder.findFusions();

        fusions = finder.getFusionCandidates().values().iterator().next();
        assertEquals(1, fusions.size());
        fusion = fusions.get(0);
        assertEquals(1, fusion.getFragments(MATCHED_JUNCTION).size());
        assertEquals(2, fusion.getFragments(DISCORDANT).size());
        assertEquals(20, fusion.getMaxSplitLengths()[SE_START]);
        assertEquals(20, fusion.getMaxSplitLengths()[SE_END]);
    }

    @Test
    public void testSoftClippedFragmentRealignment()
    {
        Configurator.setRootLevel(Level.DEBUG);

        final EnsemblDataCache geneTransCache = createGeneDataCache();

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        IsofoxConfig config = new IsofoxConfig();
        populateRefGenome(config.RefGenome);

        FusionFinder finder = new FusionFinder(config, geneTransCache);

        int gcId = 0;

        final GeneCollection gc1 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));
        final GeneCollection gc2 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_2)));

        Map<Integer,List<EnsemblGeneData>> gcMap = Maps.newHashMap();

        gcMap.put(gc1.id(), gc1.genes().stream().map(x -> x.GeneData).collect(Collectors.toList()));
        gcMap.put(gc2.id(), gc2.genes().stream().map(x -> x.GeneData).collect(Collectors.toList()));
        finder.addChromosomeGeneCollections(CHR_1, gcMap);

        final Map<String,List<ReadRecord>> chimericReadMap = Maps.newHashMap();

        // a spliced fragment to establish the fusion
        int readId = 0;
        ReadRecord read1 = createMappedRead(readId, gc1, 1050, 1089, createCigar(0, 40, 0));

        ReadRecord[] readPair = createSupplementaryReadPair(readId, gc1, gc2, 1081, 1100, 10200, 10219,
                createCigar(0, 20, 20), createCigar(20, 20, 0), true);

        readPair[1].setStrand(true, false);

        chimericReadMap.put(read1.Id, Lists.newArrayList(read1, readPair[0], readPair[1]));

        // a soft-clipped read matching the other side of the fusion junction
        String junctionBases = config.RefGenome.getBaseString(gc1.chromosome(), 1071, 1100)
                + config.RefGenome.getBaseString(gc1.chromosome(), 10200, 10209);
        ReadRecord read4 = createMappedRead(++readId, gc1, 1071, 1100, createCigar(0, 30, 10), junctionBases);
        ReadRecord read5 = createMappedRead(readId, gc1, 1051, 1090, createCigar(0, 40, 0));
        read5.setStrand(true, false);

        chimericReadMap.put(read4.Id, Lists.newArrayList(read4, read5));

        // a soft-clipped read matching 2 bases into the ref due to homology with the other side of the fusion junction
        junctionBases = config.RefGenome.getBaseString(gc1.chromosome(), 1091, 1100)
                + config.RefGenome.getBaseString(gc1.chromosome(), 10200, 10229);

        ReadRecord read6 = createMappedRead(++readId, gc2, 10198, 10229, createCigar(8, 32, 0), junctionBases);
        ReadRecord read7 = createMappedRead(readId, gc2, 10210, 10249, createCigar(0, 40, 0));
        read7.setStrand(true, false);

        chimericReadMap.put(read6.Id, Lists.newArrayList(read6, read7));

        finder.addChimericReads(chimericReadMap);

        finder.findFusions();

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

        FusionFinder finder = new FusionFinder(config, geneTransCache);

        int gcId = 0;

        final GeneCollection gc1 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));
        final GeneCollection gc2 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_2)));

        final Map<String,List<ReadRecord>> chimericReadMap = Maps.newHashMap();

        // a simple DEL, supported by 2 split fragments
        int readId = 0;
        ReadRecord read1 = createMappedRead(readId, gc1, 1081, 1100, createCigar(0, 20, 20));
        ReadRecord read2 = createMappedRead(readId, gc2, 10200, 10219, createCigar(20, 20, 0));
        read2.setStrand(true, false);

        chimericReadMap.put(read1.Id, Lists.newArrayList(read1, read2));

        // a second one
        read1 = createMappedRead(++readId, gc1, 1081, 1100, createCigar(0, 20, 20));
        read2 = createMappedRead(readId, gc2, 10200, 10219, createCigar(20, 20, 0));
        read2.setStrand(true, false);

        chimericReadMap.put(read1.Id, Lists.newArrayList(read1, read2));

        // realigned and discordant reads which support it

        // 2 soft-clipped fragments matching the other side of the fusion junction
        String junctionBases = config.RefGenome.getBaseString(gc1.chromosome(), 1071, 1100)
                + config.RefGenome.getBaseString(gc2.chromosome(), 10200, 10209);

        ReadRecord read3 = createMappedRead(++readId, gc1, 1071, 1100, createCigar(0, 30, 10), junctionBases);
        ReadRecord read4 = createMappedRead(readId, gc1, 1051, 1090, createCigar(0, 40, 0));
        read4.setStrand(true, false);
        chimericReadMap.put(read3.Id, Lists.newArrayList(read3, read4));

        junctionBases = config.RefGenome.getBaseString(gc1.chromosome(), 1091, 1100)
                + config.RefGenome.getBaseString(gc2.chromosome(), 10200, 10229);

        ReadRecord read5 = createMappedRead(++readId, gc2, 10200, 10229, createCigar(10, 30, 0), junctionBases);
        ReadRecord read6 = createMappedRead(readId, gc2, 10210, 10249, createCigar(0, 40, 0));
        read6.setStrand(true, false);
        chimericReadMap.put(read5.Id, Lists.newArrayList(read5, read6));

        // and a discordant fragment
        read3 = createMappedRead(++readId, gc1, 1050, 1089, createCigar(0, 40, 0));
        read4 = createMappedRead(readId, gc2, 10210, 10249, createCigar(0, 40, 0));
        read4.setStrand(true, false);
        chimericReadMap.put(read3.Id, Lists.newArrayList(read3, read4));

        finder.addChimericReads(chimericReadMap);

        finder.findFusions();

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
    public void testLocalSplitReadDelFusions()
    {
        // a split read (DEL) spanning 2 genes
        final EnsemblDataCache geneTransCache = createGeneDataCache();

        addTestGenes(geneTransCache);
        addTestTranscripts(geneTransCache);

        IsofoxConfig config = new IsofoxConfig();
        populateRefGenome(config.RefGenome);

        FusionFinder finder = new FusionFinder(config, geneTransCache);

        int gcId = 0;

        final GeneCollection gc1 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));
        final GeneCollection gc2 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_2)));

        final Map<String,List<ReadRecord>> chimericReadMap = Maps.newHashMap();

        int readId = 0;
        ReadRecord read1 = createMappedRead(readId, gc1, 1081, 10219, createCigar(0, 20, 9099, 20, 0));
        read1.setGeneCollection(SE_START, gc1.id(), true);
        read1.setGeneCollection(SE_END, gc2.id(), true);

        ReadRecord read2 = createMappedRead(readId, gc2, 10210, 10249, createCigar(0, 40, 0));
        read2.setStrand(true, false);

        chimericReadMap.put(read1.Id, Lists.newArrayList(read1, read2));

        // realigned and discordant reads which support it

        // 2 soft-clipped fragments matching the other side of the fusion junction
        String junctionBases = config.RefGenome.getBaseString(gc1.chromosome(), 1071, 1100)
                + config.RefGenome.getBaseString(gc2.chromosome(), 10200, 10209);

        ReadRecord read3 = createMappedRead(++readId, gc1, 1071, 1100, createCigar(0, 30, 10), junctionBases);
        ReadRecord read4 = createMappedRead(readId, gc1, 1051, 1090, createCigar(0, 40, 0));
        read4.setStrand(true, false);
        chimericReadMap.put(read3.Id, Lists.newArrayList(read3, read4));

        junctionBases = config.RefGenome.getBaseString(gc1.chromosome(), 1091, 1100)
                + config.RefGenome.getBaseString(gc2.chromosome(), 10200, 10229);
        ReadRecord read5 = createMappedRead(++readId, gc2, 10200, 10229, createCigar(10, 30, 0), junctionBases);
        ReadRecord read6 = createMappedRead(readId, gc2, 10210, 10249, createCigar(0, 40, 0));
        read6.setStrand(true, false);
        chimericReadMap.put(read5.Id, Lists.newArrayList(read5, read6));

        // and a discordant fragment
        read3 = createMappedRead(++readId, gc1, 1050, 1089, createCigar(0, 40, 0));
        read4 = createMappedRead(readId, gc2, 10210, 10249, createCigar(0, 40, 0));
        read4.setStrand(true, false);
        chimericReadMap.put(read3.Id, Lists.newArrayList(read3, read4));

        finder.addChimericReads(chimericReadMap);

        // populate the chr-gene-collection map
        Map<Integer,List<EnsemblGeneData>> gcGeneMap = Maps.newHashMap();
        gcGeneMap.put(gc1.id(), gc1.genes().stream().map(x -> x.GeneData).collect(Collectors.toList()));
        gcGeneMap.put(gc2.id(), gc2.genes().stream().map(x -> x.GeneData).collect(Collectors.toList()));
        finder.addChromosomeGeneCollections(CHR_1, gcGeneMap);

        finder.findFusions();

        assertEquals(1, finder.getFusionCandidates().size());
        List<FusionReadData> fusions = finder.getFusionCandidates().values().iterator().next();
        assertEquals(1, fusions.size());

        FusionReadData fusion = fusions.stream().filter(x -> x.isKnownSpliced()).findFirst().orElse(null);
        assertTrue(fusion != null);
        assertEquals(1, fusion.getFragments(MATCHED_JUNCTION).size());
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
        FusionFinder finder = new FusionFinder(config, geneTransCache);

        int gcId = 0;

        final GeneCollection gc1 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_1)));
        final GeneCollection gc2 = createGeneCollection(geneTransCache, gcId++, Lists.newArrayList(geneTransCache.getGeneDataById(GENE_ID_2)));

        final Map<String,List<ReadRecord>> chimericReadMap = Maps.newHashMap();

        // a simple DEL, supported by 1 split fragment
        int readId = 0;
        ReadRecord read1 = createMappedRead(readId, gc1, 581, 600, createCigar(0, 20, 20));
        ReadRecord read2 = createMappedRead(readId, gc2, 10200, 10219, createCigar(20, 20, 0));
        read2.setStrand(true, false);

        chimericReadMap.put(read1.Id, Lists.newArrayList(read1, read2));

        // realigned and discordant reads which support it

        // 2 soft-clipped fragments matching the other side of the fusion junction
        String junctionBases = config.RefGenome.getBaseString(gc1.chromosome(), 571, 600)
                + config.RefGenome.getBaseString(gc2.chromosome(), 10200, 10209);

        ReadRecord read3 = createMappedRead(++readId, gc1, 571, 600, createCigar(0, 30, 10), junctionBases);
        ReadRecord read4 = createMappedRead(readId, gc1, 551, 590, createCigar(0, 40, 0));
        read4.setStrand(true, false);
        chimericReadMap.put(read3.Id, Lists.newArrayList(read3, read4));

        junctionBases = config.RefGenome.getBaseString(gc1.chromosome(), 591, 600)
                + config.RefGenome.getBaseString(gc2.chromosome(), 10200, 10229);

        ReadRecord read5 = createMappedRead(++readId, gc2, 10200, 10229, createCigar(10, 30, 0), junctionBases);
        ReadRecord read6 = createMappedRead(readId, gc2, 10210, 10249, createCigar(0, 40, 0));
        read6.setStrand(true, false);
        chimericReadMap.put(read5.Id, Lists.newArrayList(read5, read6));

        // and a discordant fragment
        read3 = createMappedRead(++readId, gc1, 550, 589, createCigar(0, 40, 0));
        read4 = createMappedRead(readId, gc2, 10210, 10249, createCigar(0, 40, 0));
        read4.setStrand(true, false);
        chimericReadMap.put(read3.Id, Lists.newArrayList(read3, read4));

        finder.addChimericReads(chimericReadMap);

        finder.findFusions();

        assertEquals(1, finder.getFusionCandidates().size());
        List<FusionReadData> fusions = finder.getFusionCandidates().values().iterator().next();
        assertEquals(1, fusions.size());

        FusionReadData fusion = fusions.stream().findFirst().orElse(null);
        assertTrue(fusion != null);
        assertEquals(1, fusion.getFragments(MATCHED_JUNCTION).size());
        assertEquals(2, fusion.getFragments(REALIGNED).size());
        // assertEquals(1, fusion.getFragments(DISCORDANT).size()); // not yet handled
    }

}

package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.test.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.test.GeneTestUtils.generateExonStarts;
import static com.hartwig.hmftools.common.gene.TranscriptProteinData.BIOTYPE_PROTEIN_CODING;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConstants.SINGLE_MAP_QUALITY;
import static com.hartwig.hmftools.isofox.common.ReadRecord.findOverlappingRegions;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.N;
import static htsjdk.samtools.SAMFlag.FIRST_OF_PAIR;
import static htsjdk.samtools.SAMFlag.SECOND_OF_PAIR;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.fusion.ChimericReadGroup;
import com.hartwig.hmftools.isofox.fusion.FusionFinder;
import com.hartwig.hmftools.isofox.fusion.FusionWriter;
import com.hartwig.hmftools.isofox.fusion.JunctionRacFragments;
import com.hartwig.hmftools.isofox.fusion.PassingFusions;
import com.hartwig.hmftools.isofox.fusion.RacFragmentCache;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFlag;

public class TestUtils
{
    public static final String CHR_1 = "1";

    public static final String GENE_NAME_1 = "GENE1"; // +ve strand
    public static final String GENE_ID_1 = "ENSG0001";
    public static final int GENE_START_1 = 1000;
    public static final int[] GENE_1_BOUNDARIES = {0, 9999};

    public static final String GENE_NAME_2 = "GENE2"; // +ve strand
    public static final String GENE_ID_2 = "ENSG0002";
    public static final int GENE_START_2 = 10000;
    public static final int[] GENE_2_BOUNDARIES = {1501, 19999};

    public static final String GENE_NAME_3 = "GENE3"; // -ve strand
    public static final String GENE_ID_3 = "ENSG0003";
    public static final int GENE_START_3 = 20000;
    public static final int[] GENE_3_BOUNDARIES = {10901, 29999};

    public static final String CHR_2 = "2";

    public static final String GENE_NAME_4 = "GENE4"; // +ve strand
    public static final String GENE_ID_4 = "ENSG0004";
    public static final int GENE_START_4 = 1000;
    public static final int[] GENE_4_BOUNDARIES = {0, 9999};

    public static final String GENE_NAME_5 = "GENE5"; // -ve strand
    public static final String GENE_ID_5 = "ENSG0005";
    public static final int GENE_START_5 = 10000;
    public static final int[] GENE_5_BOUNDARIES = {1501, 19999};

    public static final String GENE_NAME_6 = "GENE6"; // -ve strand
    public static final String GENE_ID_6 = "ENSG0006";
    public static final int GENE_START_6 = 10400; // will overlap with previous gene and share some exons
    public static final int[] GENE_6_BOUNDARIES = {10901, 29999};

    public static final byte POS_STRAND = 1;
    public static final byte NEG_STRAND = -1;

    public static final int EXON_LENGTH = 100;
    public static final int INTRON_LENGTH = 200;

    public static IsofoxConfig createIsofoxConfig() { return new IsofoxConfig(new MockRefGenome()); }

    public static int getGeneCollection(final String geneId)
    {
        if(geneId == GENE_ID_1)
            return 0;
        if(geneId == GENE_ID_2)
            return 1;
        if(geneId == GENE_ID_3)
            return 2;
        if(geneId == GENE_ID_4)
            return 3;
        if(geneId == GENE_ID_5)
            return 4;
        if(geneId == GENE_ID_6)
            return 5;

        return -1;
    }

    public static void addTestGenes(EnsemblDataCache geneTransCache)
    {
        List<GeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(GENE_ID_1, GENE_NAME_1, CHR_1, POS_STRAND, GENE_START_1, GENE_START_1 + 500));
        geneList.add(createEnsemblGeneData(GENE_ID_2, GENE_NAME_2, CHR_1, POS_STRAND, GENE_START_2, GENE_START_2 + 900));
        geneList.add(createEnsemblGeneData(GENE_ID_3, GENE_NAME_3, CHR_1, NEG_STRAND, GENE_START_3, GENE_START_3 + 500));
        addGeneData(geneTransCache, CHR_1, geneList);

        geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(GENE_ID_4, GENE_NAME_4, CHR_2, POS_STRAND, GENE_START_4, GENE_START_4 + 500));
        geneList.add(createEnsemblGeneData(GENE_ID_5, GENE_NAME_5, CHR_2, NEG_STRAND, GENE_START_5, GENE_START_5 + 900));
        geneList.add(createEnsemblGeneData(GENE_ID_6, GENE_NAME_6, CHR_2, NEG_STRAND, GENE_START_6, GENE_START_6 + 1100));
        addGeneData(geneTransCache, CHR_2, geneList);
    }

    public static final int TRANS_1 = 1;
    public static final int TRANS_2 = 2;
    public static final int TRANS_3 = 3;
    public static final int TRANS_4 = 4;
    public static final int TRANS_5 = 5;
    public static final int TRANS_6 = 6;

    public static void addTestTranscripts(EnsemblDataCache geneTransCache)
    {
        Integer codingStart = null;
        Integer codingEnd = null;
        boolean canonical = true;

        List<TranscriptData> transDataList = Lists.newArrayList();

        // exons will be gene start + 0-100, 200-300, 400-500
        TranscriptData transData = createTransExons(
                GENE_ID_1, TRANS_1, POS_STRAND, generateExonStarts(GENE_START_1, 3, EXON_LENGTH, 100),
                EXON_LENGTH, codingStart, codingEnd, canonical, BIOTYPE_PROTEIN_CODING);

        transDataList.add(transData);

        addTransExonData(geneTransCache, GENE_ID_1, transDataList);

        transDataList = Lists.newArrayList();

        // exons will be gene start + 0-100, 200-300, 400-500, 600-700, 800-900
        transData = createTransExons(
                GENE_ID_2, TRANS_2, POS_STRAND, generateExonStarts(GENE_START_2, 5, EXON_LENGTH, 100),
                EXON_LENGTH, codingStart, codingEnd, canonical, BIOTYPE_PROTEIN_CODING);

        transDataList.add(transData);

        addTransExonData(geneTransCache, GENE_ID_2, transDataList);

        transDataList = Lists.newArrayList();

        transData = createTransExons(
                GENE_ID_3, TRANS_3, NEG_STRAND, generateExonStarts(GENE_START_3, 3, EXON_LENGTH, 100),
                EXON_LENGTH, codingStart, codingEnd, canonical, BIOTYPE_PROTEIN_CODING);

        transDataList.add(transData);

        addTransExonData(geneTransCache, GENE_ID_3, transDataList);

        transDataList = Lists.newArrayList();

        transData = createTransExons(
                GENE_ID_4, TRANS_4, POS_STRAND, generateExonStarts(GENE_START_4, 3, EXON_LENGTH, 100),
                EXON_LENGTH, codingStart, codingEnd, canonical, BIOTYPE_PROTEIN_CODING);

        transDataList.add(transData);

        addTransExonData(geneTransCache, GENE_ID_4, transDataList);

        transDataList = Lists.newArrayList();

        transData = createTransExons(
                GENE_ID_5, TRANS_5, NEG_STRAND, generateExonStarts(GENE_START_5, 5, EXON_LENGTH, 100),
                EXON_LENGTH, codingStart, codingEnd, canonical, BIOTYPE_PROTEIN_CODING);

        transDataList.add(transData);

        addTransExonData(geneTransCache, GENE_ID_5, transDataList);

        transDataList = Lists.newArrayList();

        transData = createTransExons(
                GENE_ID_6, TRANS_6, NEG_STRAND, generateExonStarts(GENE_START_6, 6, EXON_LENGTH, 100),
                EXON_LENGTH, codingStart, codingEnd, canonical, BIOTYPE_PROTEIN_CODING);

        transDataList.add(transData);

        addTransExonData(geneTransCache, GENE_ID_6, transDataList);
    }

    public static void populateRefGenome(final RefGenomeInterface refGenome)
    {
        final MockRefGenome mockRefGenome = (MockRefGenome)refGenome;
        final String refBases = generateRandomBases(GENE_START_3 + 10000);
        mockRefGenome.RefGenomeMap.put(CHR_1, refBases);
        mockRefGenome.RefGenomeMap.put(CHR_2, refBases);
    }

    public static void overrideRefGenome(final RefGenomeInterface refGenome, final String chromosome, int posStart, final String bases)
    {
        final MockRefGenome mockRefGenome = (MockRefGenome)refGenome;
        final String refBases = mockRefGenome.RefGenomeMap.get(chromosome);

        if(refBases == null || posStart + bases.length() > refBases.length())
            return;

        int posEnd = posStart + bases.length();
        final String newRefBases = refBases.substring(0, posStart) + bases + refBases.substring(posEnd);
        mockRefGenome.RefGenomeMap.put(chromosome, newRefBases);
    }

    public static Cigar createCigar(int preSc, int match, int postSc)
    {
        if(preSc == 0 && match == 0 && postSc == 0)
            return null;

        Cigar cigar = new Cigar();

        if(preSc > 0)
            cigar.add(new CigarElement(preSc, CigarOperator.SOFT_CLIP));

        if(match > 0)
            cigar.add(new CigarElement(match, CigarOperator.MATCH_OR_MISMATCH));

        if(postSc > 0)
            cigar.add(new CigarElement(postSc, CigarOperator.SOFT_CLIP));

        return cigar;
    }

    public static Cigar createCigar(int preSc, int preMatch, int nSplit, int postMatch, int postSc)
    {
        if(preSc == 0 && preMatch == 0 && nSplit == 0 && postMatch == 0 && postSc == 0)
            return null;

        Cigar cigar = new Cigar();

        if(preSc > 0)
            cigar.add(new CigarElement(preSc, CigarOperator.SOFT_CLIP));

        if(preMatch > 0)
            cigar.add(new CigarElement(preMatch, CigarOperator.MATCH_OR_MISMATCH));

        if(nSplit > 0)
            cigar.add(new CigarElement(nSplit, CigarOperator.SKIPPED_REGION));

        if(postMatch > 0)
            cigar.add(new CigarElement(postMatch, CigarOperator.MATCH_OR_MISMATCH));

        if(postSc > 0)
            cigar.add(new CigarElement(postSc, CigarOperator.SOFT_CLIP));

        return cigar;
    }

    public static ReadRecord createReadRecord(
            final int id, final String chromosome, int posStart, int posEnd, final String readBases, final Cigar cigar)
    {
        return createReadRecord(id, chromosome, posStart, posEnd, readBases, cigar, 0, chromosome, posStart);
    }

    public static ReadRecord createReadRecord(
            final int id, final String chromosome, int posStart, int posEnd, final String readBases, final Cigar cigar,
            int flags, final String mateChr, int mateStartPos)
    {
        Cigar readCigar = cigar != null ? cigar : createCigar(0, (int) (posEnd - posStart + 1), 0);

        ReadRecord read = new ReadRecord(String.valueOf(id), chromosome, posStart, posEnd, readBases, readCigar,
                0, flags, mateChr, mateStartPos);

        read.setFlag(SAMFlag.PROPER_PAIR, true);
        read.setFlag(SAMFlag.READ_PAIRED, true);
        read.setStrand(false, true);
        read.setMapQuality(SINGLE_MAP_QUALITY);
        return read;
    }

    public static ReadRecord[] createSupplementaryReadPair(final int id, final GeneCollection gc1, final GeneCollection gc2,
            int posStart1, int posEnd1, int posStart2, int posEnd2, final Cigar cigar1, final Cigar cigar2, boolean firstInPair)
    {
        int readBaseLength = cigar1.getCigarElements().stream()
                .filter(x -> x.getOperator() != N && x.getOperator() != D)
                .mapToInt(x -> x.getLength()).sum();

        String readBases = generateRandomBases(readBaseLength);

        ReadRecord read1 = createMappedRead(id, gc1, posStart1, posEnd1, cigar1, readBases);
        ReadRecord read2 = createMappedRead(id, gc2, posStart2, posEnd2, cigar2, readBases);
        read1.setFlag(FIRST_OF_PAIR, firstInPair);
        read2.setFlag(SECOND_OF_PAIR, !firstInPair);

        // note: strand is not currently set correctly
        SupplementaryReadData suppData1 = new SupplementaryReadData(
                read2.Chromosome, read2.PosStart, '+', read2.Cigar.toString(), 255);

        read1.setSuppAlignment(suppData1.asCsv());
        // read1.setSuppAlignment(String.format("%s;%d;%s", read2.Chromosome, read2.PosStart, read2.Cigar.toString()));

        SupplementaryReadData suppData2 = new SupplementaryReadData(
                read1.Chromosome, read1.PosStart, '+', read1.Cigar.toString(), 255);

        read2.setSuppAlignment(suppData2.asCsv());
        // read2.setSuppAlignment(String.format("%s;%d;%s", read1.Chromosome, read1.PosStart, read1.Cigar.toString()));

        return new ReadRecord[] { read1, read2 };
    }

    public static ReadRecord[] createReadPair(final int id, final GeneCollection gc1, final GeneCollection gc2,
            int posStart1, int posEnd1, int posStart2, int posEnd2, final Cigar cigar1, final Cigar cigar2, byte orient1, byte orient2)
    {
        int readBaseLength = cigar1.getCigarElements().stream()
                .filter(x -> x.getOperator() != N && x.getOperator() != D)
                .mapToInt(x -> x.getLength()).sum();

        String readBases = generateRandomBases(readBaseLength);

        ReadRecord read1 = createMappedRead(id, gc1, posStart1, posEnd1, cigar1, readBases);
        ReadRecord read2 = createMappedRead(id, gc2, posStart2, posEnd2, cigar2, readBases);
        read1.setFlag(FIRST_OF_PAIR, true);
        read2.setFlag(SECOND_OF_PAIR, true);

        read1.setStrand(orient1 == -1, orient2 == -1);
        read2.setStrand(orient2 == -1, orient1 == -1);

        return new ReadRecord[] { read1, read2 };
    }

    public static ReadRecord createMappedRead(final int id, final GeneCollection geneCollection, int posStart, int posEnd, final Cigar cigar)
    {
        int readBaseLength = cigar.getCigarElements().stream()
                .filter(x -> x.getOperator() != N && x.getOperator() != D)
                .mapToInt(x -> x.getLength()).sum();

        String readBases = generateRandomBases(readBaseLength);

        return createMappedRead(id, geneCollection, posStart, posEnd, cigar, readBases);
    }

    public static ReadRecord createMappedRead(
            final int id, final GeneCollection geneCollection, int posStart, int posEnd, final Cigar cigar, final String readBases)
    {
        ReadRecord read = createReadRecord(id, geneCollection.chromosome(), posStart, posEnd, readBases, cigar);

        read.processOverlappingRegions(findOverlappingRegions(geneCollection.getExonRegions(), read));

        if(read.getMappedRegions().isEmpty())
            read.addIntronicTranscriptRefs(geneCollection.getTranscripts());

        geneCollection.setReadGeneCollections(read, geneCollection.getNonGenicPositions());

        return read;
    }

    public static RegionReadData createRegion(final String geneId, int trans, int exonRank, final String chromosome, int posStart, int posEnd)
    {
        RegionReadData region = new RegionReadData(chromosome, posStart, posEnd);
        region.addExonRef(geneId, trans, String.valueOf(trans), exonRank);
        return region;
    }

    public static GeneReadData createGeneReadData(final String geneId, final String chromosome, byte strand, int posStart, int posEnd)
    {
        GeneData geneData = new GeneData(geneId, geneId, chromosome, strand, posStart, posEnd, "");
        return new GeneReadData(geneData);
    }

    public static GeneCollection createGeneCollection(final EnsemblDataCache geneTransCache, int id, final List<GeneData> genes)
    {
        List<GeneReadData> geneReadData = Lists.newArrayList();
        genes.forEach(x -> geneReadData.add(new GeneReadData(x)));
        geneReadData.forEach(x -> x.setTranscripts(geneTransCache.getTranscripts(x.GeneData.GeneId)));

        GeneCollection gc = new GeneCollection(id, geneReadData);

        if(gc.geneNames().equals(GENE_NAME_1))
        {
            gc.setNonGenicPosition(SE_START, GENE_1_BOUNDARIES[SE_START]);
            gc.setNonGenicPosition(SE_END, GENE_1_BOUNDARIES[SE_END]);
        }
        else if(gc.geneNames().equals(GENE_NAME_2))
        {
            gc.setNonGenicPosition(SE_START, GENE_2_BOUNDARIES[SE_START]);
            gc.setNonGenicPosition(SE_END, GENE_2_BOUNDARIES[SE_END]);
        }
        else if(gc.geneNames().equals(GENE_NAME_3))
        {
            gc.setNonGenicPosition(SE_START, GENE_3_BOUNDARIES[SE_START]);
            gc.setNonGenicPosition(SE_END, GENE_3_BOUNDARIES[SE_END]);
        }
        else if(gc.geneNames().equals(GENE_NAME_4))
        {
            gc.setNonGenicPosition(SE_START, GENE_4_BOUNDARIES[SE_START]);
            gc.setNonGenicPosition(SE_END, GENE_4_BOUNDARIES[SE_END]);
        }
        else if(gc.geneNames().equals(GENE_NAME_5))
        {
            gc.setNonGenicPosition(SE_START, GENE_5_BOUNDARIES[SE_START]);
            gc.setNonGenicPosition(SE_END, GENE_5_BOUNDARIES[SE_END]);
        }
        else if(gc.geneNames().equals(GENE_NAME_6))
        {
            gc.setNonGenicPosition(SE_START, GENE_6_BOUNDARIES[SE_START]);
            gc.setNonGenicPosition(SE_END, GENE_6_BOUNDARIES[SE_END]);
        }

        return gc;
    }

    public static FusionFinder createFusionFinder(final IsofoxConfig config, final EnsemblDataCache geneTransCache)
    {
        return createFusionFinder(config, geneTransCache, new RacFragmentCache());
    }

    public static FusionFinder createFusionFinder(
            final IsofoxConfig config, final EnsemblDataCache geneTransCache, final RacFragmentCache racFragmentCache)
    {
        config.Filters.buildGeneRegions(geneTransCache);

        return new FusionFinder(
                "FF", config, geneTransCache, racFragmentCache,
                new PassingFusions(config.Fusions.KnownFusions, null), new FusionWriter(config));
    }

    public static void addRacJunction(
            final RacFragmentCache racFragmentCache, final String chromosome, int gcId, byte juncOrient, int juncPosition)
    {
        JunctionRacFragments juncRacFragments = racFragmentCache.getRacFragments(chromosome, gcId);

        if(juncRacFragments == null)
        {
            juncRacFragments = new JunctionRacFragments();
            juncRacFragments.addJunction(juncPosition, juncOrient);
            racFragmentCache.addRacFragments(chromosome, gcId, juncRacFragments);
        }
    }

    public static void addRacReadGroup(
            final RacFragmentCache racFragmentCache, final ChimericReadGroup readGroup, byte juncOrient, int juncPosition)
    {
        String chromosome = readGroup.reads().get(0).Chromosome;
        int gcId = readGroup.reads().get(0).getGeneCollectons()[SE_START];

        JunctionRacFragments juncRacFragments = racFragmentCache.getRacFragments(chromosome, gcId);

        if(juncRacFragments == null)
        {
            juncRacFragments = new JunctionRacFragments();
            juncRacFragments.addJunction(juncPosition, juncOrient);
            juncRacFragments.checkAddCandidateGroup(readGroup);
            racFragmentCache.addRacFragments(chromosome, gcId, juncRacFragments);
        }
        else
        {
            juncRacFragments.checkAddCandidateGroup(readGroup);
        }
    }

}

package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.addGeneData;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.addTransExonData;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.common.ensemblcache.GeneTestUtils.generateExonStarts;
import static com.hartwig.hmftools.common.ensemblcache.TranscriptProteinData.BIOTYPE_PROTEIN_CODING;
import static com.hartwig.hmftools.isofox.ReadCountsTest.REF_BASE_STR_1;
import static com.hartwig.hmftools.isofox.common.ReadRecord.findOverlappingRegions;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.N;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionReadData;

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

    public static final String GENE_NAME_2 = "GENE2"; // +ve strand
    public static final String GENE_ID_2 = "ENSG0002";
    public static final int GENE_START_2 = 10000;

    public static final String GENE_NAME_3 = "GENE3"; // -ve strand
    public static final String GENE_ID_3 = "ENSG0003";
    public static final int GENE_START_3 = 20000;

    public static final String CHR_2 = "2";

    public static final String GENE_NAME_4 = "GENE4"; // +ve strand
    public static final String GENE_ID_4 = "ENSG0004";
    public static final int GENE_START_4 = 1000;

    public static final String GENE_NAME_5 = "GENE5"; // -ve strand
    public static final String GENE_ID_5 = "ENSG0005";
    public static final int GENE_START_5 = 10000;

    public static final byte POS_STRAND = 1;
    public static final byte NEG_STRAND = -1;

    public static final int EXON_LENGTH = 100;
    public static final int INTRON_LENGTH = 200;

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

        return -1;
    }

    public static void addTestGenes(EnsemblDataCache geneTransCache)
    {
        List<EnsemblGeneData> geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(GENE_ID_1, GENE_NAME_1, CHR_1, POS_STRAND, GENE_START_1, GENE_START_1 + 1000));
        geneList.add(createEnsemblGeneData(GENE_ID_2, GENE_NAME_2, CHR_1, POS_STRAND, GENE_START_2, GENE_START_2 + 1000));
        geneList.add(createEnsemblGeneData(GENE_ID_3, GENE_NAME_3, CHR_1, NEG_STRAND, GENE_START_3, GENE_START_3 + 1000));
        addGeneData(geneTransCache, CHR_1, geneList);

        geneList = Lists.newArrayList();
        geneList.add(createEnsemblGeneData(GENE_ID_4, GENE_NAME_4, CHR_2, POS_STRAND, GENE_START_4, GENE_START_4 + 1000));
        geneList.add(createEnsemblGeneData(GENE_ID_5, GENE_NAME_5, CHR_2, NEG_STRAND, GENE_START_5, GENE_START_5 + 1000));
        addGeneData(geneTransCache, CHR_2, geneList);
    }

    public static final int TRANS_1 = 1;
    public static final int TRANS_2 = 2;
    public static final int TRANS_3 = 3;
    public static final int TRANS_4 = 4;
    public static final int TRANS_5 = 5;

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
                GENE_ID_5, TRANS_5, NEG_STRAND, generateExonStarts(GENE_START_5, 3, EXON_LENGTH, 100),
                EXON_LENGTH, codingStart, codingEnd, canonical, BIOTYPE_PROTEIN_CODING);

        transDataList.add(transData);

        addTransExonData(geneTransCache, GENE_ID_5, transDataList);
    }

    public static Cigar createCigar(int preSc, int match, int postSc)
    {
        if (preSc == 0 && match == 0 && postSc == 0)
            return null;

        Cigar cigar = new Cigar();

        if (preSc > 0)
            cigar.add(new CigarElement(preSc, CigarOperator.SOFT_CLIP));

        if (match > 0)
            cigar.add(new CigarElement(match, CigarOperator.MATCH_OR_MISMATCH));

        if (postSc > 0)
            cigar.add(new CigarElement(postSc, CigarOperator.SOFT_CLIP));

        return cigar;

    }

    public static Cigar createCigar(int preSc, int preMatch, int nonRef, int postMatch, int postSc)
    {
        if (preSc == 0 && preMatch == 0 && nonRef == 0 && postMatch == 0 && postSc == 0)
            return null;

        Cigar cigar = new Cigar();

        if (preSc > 0)
            cigar.add(new CigarElement(preSc, CigarOperator.SOFT_CLIP));

        if (preMatch > 0)
            cigar.add(new CigarElement(preMatch, CigarOperator.MATCH_OR_MISMATCH));

        if (nonRef > 0)
            cigar.add(new CigarElement(nonRef, CigarOperator.SKIPPED_REGION));

        if (postMatch > 0)
            cigar.add(new CigarElement(postMatch, CigarOperator.MATCH_OR_MISMATCH));

        if (postSc > 0)
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
        read.setStrand(false, true);
        return read;
    }

    public static ReadRecord createMappedRead(final int id, final GeneCollection geneCollection, int posStart, int posEnd, final Cigar cigar)
    {
        int readBaseLength = cigar.getCigarElements().stream()
                .filter(x -> x.getOperator() != N || x.getOperator() != D)
                .mapToInt(x -> x.getLength()).sum();

        String readBases = generateRandomBases(readBaseLength);

        return createMappedRead(id, geneCollection, posStart, posEnd, cigar, readBases);
    }

    public static ReadRecord createMappedRead(final int id, final GeneCollection geneCollection, int posStart, int posEnd,
            final Cigar cigar, final String readBases)
    {
        ReadRecord read = createReadRecord(id, geneCollection.chromosome(), posStart, posEnd, readBases, cigar);

        read.processOverlappingRegions(findOverlappingRegions(geneCollection.getExonRegions(), read));

        if(read.getMappedRegions().isEmpty())
            read.addIntronicTranscriptRefs(geneCollection.getTranscripts());

        read.captureGeneInfo(geneCollection.id());

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
        EnsemblGeneData geneData = new EnsemblGeneData(geneId, geneId, chromosome, strand, posStart, posEnd, "");
        return new GeneReadData(geneData);
    }

    public static GeneCollection createGeneCollection(final EnsemblDataCache geneTransCache, int id, final List<EnsemblGeneData> genes)
    {
        List<GeneReadData> geneReadData = Lists.newArrayList();
        genes.forEach(x -> geneReadData.add(new GeneReadData(x)));
        geneReadData.forEach(x -> x.setTranscripts(geneTransCache.getTranscripts(x.GeneData.GeneId)));

        GeneCollection gc = new GeneCollection(id, geneReadData);
        return gc;
    }

    public static String generateRandomBases(int length)
    {
        char[] str = new char[length];
        String bases = "ACGT";


        int baseIndex = 0;
        for(int i = 0; i < length; ++i)
        {
            str[i] = bases.charAt(baseIndex);

            if(baseIndex == 3)
                baseIndex = 0;
            else
                ++baseIndex;
        }

        return String.valueOf(str);
    }

}

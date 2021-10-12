package com.hartwig.hmftools.common.test;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_0;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.codingBasesToPhase;
import static com.hartwig.hmftools.common.gene.TranscriptUtils.tickPhaseForward;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.BreakendGeneData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public final class GeneTestUtils
{
    public static final String CHR_1 = "1";
    public static final String CHR_2 = "2";
    public static final String CHR_3 = "3";

    public static final String GENE_ID_1 = "ENSG001";
    public static final String GENE_ID_2 = "ENSG002";
    public static final String GENE_ID_3 = "ENSG003";
    public static final String GENE_NAME_1 = "GENE_1";
    public static final String GENE_NAME_2 = "GENE_2";
    public static final String GENE_NAME_3 = "GENE_3";

    public static final int TRANS_ID_1 = 1;
    public static final int TRANS_ID_2 = 2;
    public static final int TRANS_ID_3 = 3;

    public static BreakendGeneData createGeneAnnotation(int svId, boolean isStart, final String geneName, String stableId, byte strand,
            final String chromosome, int position, int orientation)
    {
        String karyotypeBand = "";

        BreakendGeneData gene = new BreakendGeneData(svId, isStart, geneName, stableId, strand, karyotypeBand);
        gene.setPositionalData(chromosome, position, (byte)orientation);

        return gene;
    }

    // Ensembl data types
    public static EnsemblDataCache createGeneDataCache()
    {
        return new EnsemblDataCache("", RefGenomeVersion.V37);
    }

    public static GeneData createEnsemblGeneData(String geneId, String geneName, String chromosome, int strand, int geneStart, int geneEnd)
    {
        return new GeneData(geneId, geneName, chromosome, (byte)strand, geneStart, geneEnd,  "");
    }

    public static void addTransExonData(EnsemblDataCache geneTransCache, final String geneId, List<TranscriptData> transDataList)
    {
        geneTransCache.getTranscriptDataMap().put(geneId, transDataList);
    }

    public static void addGeneData(EnsemblDataCache geneTransCache, final String chromosome, List<GeneData> geneDataList)
    {
        geneTransCache.getChrGeneDataMap().put(chromosome, geneDataList);
    }

    public static int getCodingBases(final Integer start, final Integer end)
    {
        if(start != null && end != null)
            return (end - start) + 1;
        return 0;
    }

    public static int[] generateExonStarts(int startBase, int exonCount, int exonLength, int intronLength)
    {
        int[] exonStarts = new int[exonCount];

        for(int i = 0; i < exonCount; ++i)
        {
            exonStarts[i] = startBase + i * (exonLength + intronLength);
        }

        return exonStarts;
    }

    public static String generateTransName(int transId) { return String.format("TRAN%04d", transId); }

    public static TranscriptData createTransExons(
            final String geneId, int transId, byte strand,
            int[] exonStarts, int exonLength, Integer codingStart, Integer codingEnd, boolean isCanonical, final String biotype)
    {
        if(exonStarts.length == 0 || exonLength <= 0)
            return null;

        int exonCount = exonStarts.length;
        int transStart = exonStarts[0];
        int transEnd = exonStarts[exonCount-1] + exonLength;

        final List<ExonData> exons = Lists.newArrayList();

        boolean hasCodingBases = codingStart != null && codingEnd != null;

        // work out phases based on coding start & end
        boolean inCoding = false;
        boolean finishedCoding = false;
        int lastExonEndPhase = PHASE_NONE;

        if(strand == POS_STRAND)
        {
            for (int i = 0; i < exonCount; ++i)
            {
                int exonStart = exonStarts[i];
                int exonEnd = exonStart + exonLength;
                int exonRank = i + 1;

                int exonPhase = PHASE_0;
                int exonStartPhase = PHASE_NONE;
                int exonEndPhase = PHASE_NONE;
                int exonCodingStart = 0;

                if (hasCodingBases && !finishedCoding)
                {
                    if (!inCoding)
                    {
                        if (codingStart <= exonEnd)
                        {
                            // coding starts in this exon
                            inCoding = true;
                            exonPhase = PHASE_1; // pre-exon base is 1 less so coding start with 1
                            exonStartPhase = codingStart == exonStart ? PHASE_0 : PHASE_NONE;
                            exonCodingStart = codingStart;
                        }
                    }
                    else
                    {
                        exonPhase = tickPhaseForward(lastExonEndPhase); // tick forward
                        exonStartPhase = lastExonEndPhase; // by convention
                        exonCodingStart = exonStart;
                    }

                    if (inCoding)
                    {
                        int exonCodingBases = min(exonEnd, codingEnd) - exonCodingStart + 1;
                        exonEndPhase = codingBasesToPhase(exonCodingBases, exonPhase);
                        lastExonEndPhase = exonEndPhase;

                        if (codingEnd <= exonEnd)
                        {
                            finishedCoding = true;
                            inCoding = false;
                            exonEndPhase = PHASE_NONE;
                        }
                    }
                }

                exons.add(new ExonData(transId, exonStart, exonEnd, exonRank, exonStartPhase, exonEndPhase));
            }
        }
        else
        {
            for (int i = exonCount-1; i >= 0; --i)
            {
                int exonStart = exonStarts[i];
                int exonEnd = exonStart + exonLength;
                int exonRank = exonCount - i;

                int exonPhase = PHASE_0;
                int exonStartPhase = PHASE_NONE;
                int exonEndPhase = PHASE_NONE;
                int exonCodingEnd = 0;

                if (hasCodingBases && !finishedCoding)
                {
                    if (!inCoding)
                    {
                        if (codingEnd >= exonStart)
                        {
                            inCoding = true;
                            exonPhase = PHASE_1;
                            exonStartPhase = codingEnd == exonEnd ? PHASE_0 : PHASE_NONE;
                            exonCodingEnd = codingEnd;
                        }
                    }
                    else
                    {
                        exonPhase = tickPhaseForward(lastExonEndPhase);
                        exonStartPhase = lastExonEndPhase;
                        exonCodingEnd = exonEnd;
                    }

                    if (inCoding)
                    {
                        int exonCodingBases = exonCodingEnd - max(exonStart, codingStart) + 1;
                        exonEndPhase = codingBasesToPhase(exonCodingBases, exonPhase);
                        lastExonEndPhase = exonEndPhase;

                        if (codingStart >= exonStart)
                        {
                            finishedCoding = false;
                            inCoding = false;
                            exonEndPhase = PHASE_NONE;
                        }
                    }
                }

                exons.add(0, new ExonData(transId, exonStart, exonEnd, exonRank, exonStartPhase, exonEndPhase));
            }
        }

        TranscriptData transData = new TranscriptData(transId, generateTransName(transId), geneId, isCanonical, strand, transStart, transEnd,
                codingStart, codingEnd, biotype);

        transData.setExons(exons);

        return transData;

    }

}

package com.hartwig.hmftools.linx.fusion;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.fusion.GeneRangeData.GENE_PHASING_REGION_5P_UTR;
import static com.hartwig.hmftools.linx.fusion.GeneRangeData.GENE_PHASING_REGION_CODING_0;
import static com.hartwig.hmftools.linx.fusion.GeneRangeData.GENE_PHASING_REGION_CODING_1;
import static com.hartwig.hmftools.linx.fusion.GeneRangeData.GENE_PHASING_REGION_CODING_2;
import static com.hartwig.hmftools.linx.fusion.GeneRangeData.GENE_PHASING_REGION_MAX;
import static com.hartwig.hmftools.linx.fusion.GeneRangeData.GENE_PHASING_REGION_PROMOTOR;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptExonData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class FusionLikelihood
{
    private SvGeneTranscriptCollection mGeneTransCache;
    private String mOutputDir;

    private List<Long> mDelBucketLengths;
    private List<Long> mDupBucketLengths;

    private final Map<String, List<GeneRangeData>> mChrForwardGeneDataMap;
    private final Map<String, List<GeneRangeData>> mChrReverseGeneDataMap;

    private static final String DEL_BUCKET_LENGTHS = "fl_del_bucket_lengths";
    private static final String DUP_BUCKET_LENGTHS = "fl_dup_bucket_lengths";

    private static final Logger LOGGER = LogManager.getLogger(SvGeneTranscriptCollection.class);

    public FusionLikelihood()
    {
        mDelBucketLengths = Lists.newArrayList();
        mDupBucketLengths= Lists.newArrayList();
        mChrForwardGeneDataMap = Maps.newHashMap();
        mChrReverseGeneDataMap = Maps.newHashMap();
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(DEL_BUCKET_LENGTHS, true, "Semi-colon separated DEL bucket lengths");
        options.addOption(DUP_BUCKET_LENGTHS, true, "Semi-colon separated DUP bucket lengths");
    }

    public void initialise(final CommandLine cmdLineArgs, final String outputDir, final SvGeneTranscriptCollection geneTransCache)
    {
        mOutputDir = outputDir;
        mGeneTransCache = geneTransCache;

        if(cmdLineArgs.hasOption(DEL_BUCKET_LENGTHS))
        {
            final String delLengthData = cmdLineArgs.getOptionValue(DEL_BUCKET_LENGTHS);
            Arrays.stream(delLengthData.split(";")).forEach(x -> mDelBucketLengths.add(Long.parseLong(x)));
        }

        if(cmdLineArgs.hasOption(DUP_BUCKET_LENGTHS))
        {
            final String dupLengthData = cmdLineArgs.getOptionValue(DUP_BUCKET_LENGTHS);
            Arrays.stream(dupLengthData.split(";")).forEach(x -> mDupBucketLengths.add(Long.parseLong(x)));
        }
    }

    public void generateGeneRangeData()
    {
        // first walk forwards through the genes on the positive strand:
        // for DELs taking the first as upstream and the second as downstream for each of the bucket ranges
        // for DUPs taking the second as upstream and the first as downstream
        //


        final Map<String, List<EnsemblGeneData>> chrGeneDataMap = mGeneTransCache.getChrGeneDataMap();

        for(Map.Entry<String, List<EnsemblGeneData>> entry : chrGeneDataMap.entrySet())
        {
            for(final EnsemblGeneData geneData :entry.getValue())
            {
                final List<TranscriptExonData> transExonDataList = mGeneTransCache.getTransExonData(geneData.GeneId);

            }
        }

    }


    public void generateGenePhasingCounts()
    {
        final Map<String, List<EnsemblGeneData>> chrGeneDataMap = mGeneTransCache.getChrGeneDataMap();

        for(Map.Entry<String, List<EnsemblGeneData>> entry : chrGeneDataMap.entrySet())
        {
            final String chromosome = entry.getKey();

            List<GeneRangeData> geneList = Lists.newArrayList();

            for(final EnsemblGeneData geneData :entry.getValue())
            {
                final List<TranscriptExonData> transExonDataList = mGeneTransCache.getTransExonData(geneData.GeneId);

                final int[] phasingCounts = getGenePhasingCounts(geneData, transExonDataList);

                GeneRangeData geneRangeData = new GeneRangeData(geneData);
                geneRangeData.setPhasingCounts(phasingCounts);
                geneList.add(geneRangeData);
            }

            mChrForwardGeneDataMap.put(chromosome, geneList);
        }
    }

    @NotNull
    public static int[] getGenePhasingCounts(final EnsemblGeneData geneData, final List<TranscriptExonData> transExonDataList)
    {
        long geneStart = geneData.GeneStart;
        long geneEnd = geneData.GeneEnd;
        int geneLength = (int) (geneEnd - geneStart + 1);

        boolean[][] geneBases = new boolean[geneLength][GENE_PHASING_REGION_MAX];

        boolean hasCodingExons = false;
        int teIndex = 0;
        List<TranscriptExonData> transcriptExons = SvGeneTranscriptCollection.nextTranscriptExons(transExonDataList, teIndex);

        while (!transcriptExons.isEmpty())
        {
            if (geneData.Strand == 1)
            {
                for (int i = 0; i < transcriptExons.size(); ++i)
                {
                    final TranscriptExonData exonData = transcriptExons.get(i);

                    if (exonData.CodingStart == null)
                        break;

                    hasCodingExons = true;

                    if (exonData.ExonStart > exonData.CodingEnd) // past end of coding region
                        break;

                    // first mark the entire pre-coding region as 5' UTR
                    if (i == 0)
                    {
                        for (long j = geneStart; j <= exonData.CodingStart; ++j)
                        {
                            int gbPos = (int) (j - geneStart);
                            geneBases[gbPos][GENE_PHASING_REGION_5P_UTR] = true;
                        }
                    }

                    if (exonData.ExonEnd < exonData.CodingStart) // already accounted for
                        continue;

                    // now handle the exon's phasing
                    long codingStart = max(exonData.ExonStart, exonData.CodingStart);
                    for (long j = codingStart; j <= exonData.ExonEnd; ++j)
                    {
                        int gbPos = (int) (j - geneStart);

                        if (j > exonData.CodingEnd)
                            break;

                        long adjustedPhase = exonData.ExonPhase + (j - codingStart);
                        int calcPhase = (int) (adjustedPhase % 3);
                        geneBases[gbPos][phaseToRegion(calcPhase)] = true;
                    }

                    // fill in the intronic phasing between coding exons
                    if (i < transcriptExons.size() - 1)
                    {
                        final TranscriptExonData nextExon = transcriptExons.get(i + 1);

                        if (nextExon.ExonStart <= nextExon.CodingEnd)
                        {
                            int regionType = phaseToRegion(exonData.ExonPhaseEnd);

                            for (long j = exonData.ExonEnd + 1; j < nextExon.ExonStart; ++j)
                            {
                                int gbPos = (int) (j - geneStart);
                                geneBases[gbPos][regionType] = true;
                            }
                        }
                    }
                }
            }
            else
            {
                // navigate through as per the exon rank
                for (int i = transcriptExons.size() - 1; i >= 0; --i)
                {
                    final TranscriptExonData exonData = transcriptExons.get(i);

                    if (exonData.CodingStart == null)
                        break;

                    hasCodingExons = true;

                    if (exonData.ExonEnd < exonData.CodingStart) // past end of coding region
                        break;

                    // first mark the entire pre-coding region as 5' UTR
                    if (i == transcriptExons.size() - 1)
                    {
                        for (long j = exonData.CodingEnd; j <= geneEnd; ++j)
                        {
                            int gbPos = (int) (geneEnd - j);
                            geneBases[gbPos][GENE_PHASING_REGION_5P_UTR] = true;
                        }
                    }

                    if (exonData.ExonStart > exonData.CodingEnd) // already accounted for
                        continue;

                    // allocate the exon's phasing - working backwwards this time
                    long codingEnd = min(exonData.ExonEnd, exonData.CodingEnd);
                    for (long j = codingEnd; j >= exonData.ExonStart; --j)
                    {
                        int gbPos = (int) (geneEnd - j);

                        if (j < exonData.CodingStart)
                            break;

                        long adjustedPhase = exonData.ExonPhase + (codingEnd - j);
                        int calcPhase = (int) (adjustedPhase % 3);
                        geneBases[gbPos][phaseToRegion(calcPhase)] = true;
                    }

                    // fill in the intronic phasing between coding exons
                    if (i > 0)
                    {
                        final TranscriptExonData nextExon = transcriptExons.get(i - 1);

                        if (nextExon.ExonEnd >= nextExon.CodingStart)
                        {
                            int regionType = phaseToRegion(exonData.ExonPhaseEnd);

                            for (long j = nextExon.ExonEnd + 1; j < exonData.ExonStart; ++j)
                            {
                                int gbPos = (int) (geneEnd - j);
                                geneBases[gbPos][regionType] = true;
                            }
                        }
                    }
                }
            }

            teIndex += transcriptExons.size();
            transcriptExons = SvGeneTranscriptCollection.nextTranscriptExons(transExonDataList, teIndex);
        }

        // now compute the number of bases for each phasing region
        int[] regionTotals = new int[GENE_PHASING_REGION_MAX];

        for(int i = 0; i < geneLength; ++i)
        {
            for(int j = 0; j < GENE_PHASING_REGION_MAX; ++j)
            {
                if(geneBases[i][j])
                    ++regionTotals[j];
            }
        }
        if(hasCodingExons)
        {
            LOGGER.debug("gene({}) length({}) region counts: pre-coding({}) phases(0={} 1={} 2={})",
                    geneData.GeneId, geneData.GeneName, geneLength, regionTotals[GENE_PHASING_REGION_5P_UTR],
                    regionTotals[GENE_PHASING_REGION_CODING_0], regionTotals[GENE_PHASING_REGION_CODING_1],
                    regionTotals[GENE_PHASING_REGION_CODING_2]);
        }

        return regionTotals;
    }

    public static int phaseToRegion(int phase)
    {
        switch(phase)
        {
            case -1: return GENE_PHASING_REGION_5P_UTR;
            case 0: return GENE_PHASING_REGION_CODING_0;
            case 1: return GENE_PHASING_REGION_CODING_1;
            case 2: return GENE_PHASING_REGION_CODING_2;
        }

        return GENE_PHASING_REGION_5P_UTR;
    }

    public static int regionToPhase(int region)
    {
        if(region == GENE_PHASING_REGION_5P_UTR) return -1;
        if(region == GENE_PHASING_REGION_CODING_0) return 0;
        if(region == GENE_PHASING_REGION_CODING_1) return 1;
        if(region == GENE_PHASING_REGION_CODING_2) return 2;
        if(region == GENE_PHASING_REGION_PROMOTOR) return -1;

        return -1;
    }

}

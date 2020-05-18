package com.hartwig.hmftools.isofox.novel.cohort;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RnaUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunction.ASJ_TRANS_NONE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunction;

public class SpliceVariantMatching
{
    private final CohortConfig mConfig;
    private final Map<String,Integer> mFieldsMap;
    private final EnsemblDataCache mGeneTransCache;

    private Map<String,List<SpliceVariant>> mSampleSpliceVariants;

    private BufferedWriter mWriter;

    public SpliceVariantMatching(final CohortConfig config, final String spliceVariantFile)
    {
        mConfig = config;
        mFieldsMap = Maps.newHashMap();

        mGeneTransCache = new EnsemblDataCache(mConfig.EnsemblDataCache, RefGenomeVersion.HG37);
        mGeneTransCache.setRequiredData(true, false, false, false);
        mGeneTransCache.load(false);
        mSampleSpliceVariants = Maps.newHashMap();
        mWriter = null;

        initialiseWriter();
        loadSpliceVariants(spliceVariantFile);
    }

    public void close() { closeBufferedWriter(mWriter); }

    public void evaluateSpliceVariants(final String sampleId, final List<AltSpliceJunction> altSpliceJunctions)
    {
        if(mGeneTransCache == null || mSampleSpliceVariants.isEmpty())
            return;

        final List<SpliceVariant> spliceVariants = mSampleSpliceVariants.get(sampleId);

        if(spliceVariants == null || spliceVariants.isEmpty())
            return;

        ISF_LOGGER.debug("sampleId({}) evaluating {} splice variants",
                sampleId, spliceVariants.size());

        spliceVariants.forEach(x -> evaluateSpliceVariant(sampleId, altSpliceJunctions, x));
    }

    private void evaluateSpliceVariant(final String sampleId, final List<AltSpliceJunction> altSpliceJunctions, final SpliceVariant variant)
    {
        // search the specified transcript for the next splice junction - alt or otherwise
        /*
        EnsemblGeneData geneData = mGeneTransCache.getGeneDataByName(variant.GeneName);

        if(geneData == null)
        {
            ISF_LOGGER.warn("sampleId({}) variant({}:{}) cannot match gene({})",
                    sampleId, variant.Chromosome, variant.Position, variant.GeneName);
            return;
        }

        TranscriptData transData = mGeneTransCache.getTranscriptData(geneData.GeneId, variant.TransName);

        if(transData == null)
        {
            for(final List<TranscriptData> transList : mGeneTransCache.getTranscriptDataMap().values())
            {
                transData = transList.stream().filter(x -> x.TransName.equals(variant.TransName)).findFirst().orElse(null);
                if(transData != null)
                    break;
            }

            if(transData == null)
            {
                ISF_LOGGER.warn("sampleId({}) variant({}:{}) cannot match transcript({})",
                        sampleId, variant.Chromosome, variant.Position, variant.TransName);
                return;
            }

            geneData = mGeneTransCache.getGeneDataById(transData.GeneId);
        }
        */

        final Map<EnsemblGeneData,List<TranscriptData>> candidateGeneTrans = findCandidateTranscripts(variant);

        /*
        The splice region variant will be within 10 base of an exon.
        It may cause a novel splice junction anywhere in that exon or the intron.
        It could also cause the exon to be skipped entirely. eg. the splice variant is 7 bases after exon 3.
        Then you should check for novel splice junctions that are wholly contained within the entire region from end of exon 2 to start of exon 4.
        */

        boolean foundMatch = false;

        for(final Map.Entry<EnsemblGeneData,List<TranscriptData>> entry : candidateGeneTrans.entrySet())
        {
            final EnsemblGeneData geneData = entry.getKey();

            TranscriptData topTransData = null;
            int[] topExonMatchData = null;

            for(final TranscriptData transData : entry.getValue())
            {
                int[] exonMatchData = getTranscriptExonMatchData(transData, variant.Position, geneData.Strand);

                if(exonMatchData[ED_REGION_START] == -1 || exonMatchData[ED_REGION_END] == -1)
                    continue;

                if(topTransData == null || exonMatchData[ED_MIN_DISTANCE] < topExonMatchData[ED_MIN_DISTANCE])
                {
                    topTransData = transData;
                    topExonMatchData = exonMatchData;
                }
            }

            final List<AltSpliceJunction> candidateAltSJs = Lists.newArrayList();

            for(final AltSpliceJunction altSJ : altSpliceJunctions)
            {
                if(!altSJ.getGeneId().equals(geneData.GeneId))
                    continue;

                for(int se = SE_START; se <= SE_END; ++se)
                {
                    if (positionWithin(altSJ.SpliceJunction[se], topExonMatchData[ED_REGION_START], topExonMatchData[ED_REGION_END])
                    && (altSJ.getTranscriptNames()[se].equals(ASJ_TRANS_NONE) || altSJ.getTranscriptNames()[se].contains(topTransData.TransName)))
                    {
                        candidateAltSJs.add(altSJ);

                        ISF_LOGGER.debug("sampleId({}) variant({}:{}) gene({}) transcript({}) matched altSJ({})",
                                sampleId, variant.Chromosome, variant.Position, geneData.GeneName, topTransData.TransName, altSJ.toString());

                        break;
                    }
                }
            }

            if(!candidateAltSJs.isEmpty())
            {
                foundMatch = true;
                writeMatchData(sampleId, variant, geneData, topTransData, topExonMatchData, candidateAltSJs);
            }
        }

        if(!foundMatch)
        {
            final EnsemblGeneData geneData = mGeneTransCache.getGeneDataByName(variant.GeneName);
            writeMatchData(sampleId, variant, geneData, null, null, Lists.newArrayList());
        }
    }

    private static final int MAX_EXON_SPLICE_DISTANCE = 20;

    private Map<EnsemblGeneData,List<TranscriptData>> findCandidateTranscripts(final SpliceVariant variant)
    {
        final List<EnsemblGeneData> geneDataList = mGeneTransCache.getChrGeneDataMap().get(variant.Chromosome);
        final Map<EnsemblGeneData,List<TranscriptData>> matchedTrans = Maps.newHashMap();

        for(final EnsemblGeneData geneData : geneDataList)
        {
            if(geneData.Strand == 1)
            {
                if (!positionWithin(variant.Position, geneData.GeneStart - MAX_EXON_SPLICE_DISTANCE, geneData.GeneEnd))
                    continue;
            }
            else
            {
                if(!positionWithin(variant.Position, geneData.GeneStart, geneData.GeneEnd + MAX_EXON_SPLICE_DISTANCE))
                    continue;
            }

            final List<TranscriptData> transList = mGeneTransCache.getTranscripts(geneData.GeneId);

            final List<TranscriptData> candidateTrans = Lists.newArrayList();

            for(final TranscriptData transData : transList)
            {
                if(transData.exons().stream().anyMatch(x -> isCandidateExon(x, variant.Position)))
                {
                    candidateTrans.add(transData);
                }
            }

            if(!candidateTrans.isEmpty())
                matchedTrans.put(geneData, candidateTrans);
        }

        return matchedTrans;
    }

    private boolean isCandidateExon(final ExonData exon, int position)
    {
        if(positionWithin(position, exon.ExonStart, exon.ExonEnd))
            return true;

        int distanceFromExon = min(abs(position - exon.ExonStart), abs(position - exon.ExonEnd));
        return distanceFromExon <= MAX_EXON_SPLICE_DISTANCE;
    }

    private static final int ED_REGION_START = 0;
    private static final int ED_REGION_END = 1;
    private static final int ED_RANK_START = 2;
    private static final int ED_RANK_END = 3;
    private static final int ED_CONTEXT = 4;
    private static final int ED_MIN_DISTANCE = 5;

    private static final int ED_CONTEXT_PRE = 0;
    private static final int ED_CONTEXT_POST = 1;
    private static final int ED_CONTEXT_INTRONIC = 2;
    private static final int ED_CONTEXT_EXONIC = 3;

    private int[] getTranscriptExonMatchData(final TranscriptData transData, int splicePosition, int geneStrand)
    {
        final int[] exonMatchData = new int[ED_MIN_DISTANCE+1];

        int exonCount = transData.exons().size();
        if(splicePosition < transData.exons().get(0).ExonStart)
        {
            exonMatchData[ED_CONTEXT] = geneStrand == 1 ? ED_CONTEXT_PRE : ED_CONTEXT_POST;
            exonMatchData[ED_RANK_START] = exonMatchData[ED_RANK_END] = transData.exons().get(0).ExonRank;
            exonMatchData[ED_REGION_START] = splicePosition;
            exonMatchData[ED_REGION_END] = exonCount > 1 ? transData.exons().get(1).ExonStart : transData.exons().get(0).ExonEnd;
            exonMatchData[ED_MIN_DISTANCE] = splicePosition - transData.exons().get(0).ExonStart;
        }
        else if(splicePosition > transData.exons().get(exonCount - 1).ExonEnd)
        {
            exonMatchData[ED_CONTEXT] = geneStrand == -1 ? ED_CONTEXT_PRE : ED_CONTEXT_POST;
            exonMatchData[ED_RANK_START] = exonMatchData[ED_RANK_END] = transData.exons().get(exonCount - 1).ExonRank;
            exonMatchData[ED_REGION_START] = exonCount > 1 ? transData.exons().get(exonCount - 1).ExonEnd : transData.exons().get(exonCount - 1).ExonStart;
            exonMatchData[ED_REGION_END] = splicePosition;
            exonMatchData[ED_MIN_DISTANCE] = splicePosition - transData.exons().get(exonCount - 1).ExonEnd;
        }
        else
        {
            for (int i = 0; i < exonCount; ++i)
            {
                final ExonData prevExon = i > 0 ? transData.exons().get(i - 1) : null;
                final ExonData exon = transData.exons().get(i);
                final ExonData nextExon = i < exonCount - 1 ? transData.exons().get(i + 1) : null;

                if (positionWithin(splicePosition, exon.ExonStart, exon.ExonEnd))
                {
                    exonMatchData[ED_RANK_START] = exonMatchData[ED_RANK_END] = exon.ExonRank;
                    exonMatchData[ED_REGION_START] = prevExon != null ? prevExon.ExonEnd : exon.ExonStart;
                    exonMatchData[ED_REGION_END] = nextExon != null ? nextExon.ExonStart : exon.ExonEnd;
                    exonMatchData[ED_CONTEXT] = ED_CONTEXT_EXONIC;
                    exonMatchData[ED_MIN_DISTANCE] = min(abs(splicePosition - exon.ExonStart), abs(splicePosition - exon.ExonEnd));
                    break;
                }

                if (exon.ExonEnd < splicePosition && nextExon != null && nextExon.ExonStart > splicePosition)
                {
                    exonMatchData[ED_RANK_START] = exon.ExonRank;
                    exonMatchData[ED_RANK_END] = nextExon.ExonRank;
                    exonMatchData[ED_REGION_START] = prevExon != null ? prevExon.ExonEnd : exon.ExonStart;
                    exonMatchData[ED_REGION_END] = nextExon.ExonStart;
                    exonMatchData[ED_CONTEXT] = ED_CONTEXT_INTRONIC;
                    exonMatchData[ED_MIN_DISTANCE] = min(abs(splicePosition - exon.ExonEnd), abs(splicePosition - nextExon.ExonStart));
                    break;
                }
            }
        }

        return exonMatchData;
    }

    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("splice_variant_matching.csv");
            mWriter = createBufferedWriter(outputFileName, false);

            mWriter.write("SampleId,GeneId,GeneName,Strand,TransName,SnpeffTransMatch");

            mWriter.write(",Chromosome,Position,CodingEffect,Ref,Alt,ExonRankStart,ExonRankEnd,Context");
            mWriter.write(",AsjType,AsjContextStart,AsjContextEnd,AsjStart,AsjEnd,AsjFragmentCount,AsjDepthStart,AsjDepthEnd");

            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write splice variant file: {}", e.toString());
        }
    }

    private void writeMatchData(
            final String sampleId, final SpliceVariant variant, final EnsemblGeneData geneData,
            final TranscriptData transData, final int[] exonMatchData, final List<AltSpliceJunction> candidateAltSJs)
    {
        try
        {
            if(candidateAltSJs.isEmpty() || transData == null || exonMatchData == null)
            {
                ISF_LOGGER.debug("sampleId({}) variant({}:{}) gene({}) transcript({}) finds no candidate alt-SJs",
                        sampleId, variant.Chromosome, variant.Position, variant.GeneName, variant.TransName);

                mWriter.write(String.format("%s,%s,%s,%d,%s,%s",
                        sampleId, geneData.GeneId, geneData.GeneName, geneData.Strand,
                        variant.TransName, "true"));

                mWriter.write(String.format(",%s,%d,%s,%s,%s",
                        variant.Chromosome, variant.Position, variant.CodingEffect, variant.Ref, variant.Alt));

                mWriter.write(",-1,-1,NONE");
                mWriter.write(",NONE,NONE,NONE,-1,-1,0,0,0");
                mWriter.newLine();

                return;
            }

            for(final AltSpliceJunction altSJ : candidateAltSJs)
            {
                mWriter.write(String.format("%s,%s,%s,%d,%s,%s",
                        sampleId, geneData.GeneId, geneData.GeneName, geneData.Strand,
                        transData.TransName, transData.TransName.equals(variant.TransName)));

                mWriter.write(String.format(",%s,%d,%s,%s,%s",
                        variant.Chromosome, variant.Position, variant.CodingEffect, variant.Ref, variant.Alt));

                String variantContext = "";

                if(exonMatchData[ED_CONTEXT] == ED_CONTEXT_PRE)
                    variantContext = "PRE_TRANS";
                else if(exonMatchData[ED_CONTEXT] == ED_CONTEXT_POST)
                    variantContext = "POST_TRANS";
                else if(exonMatchData[ED_CONTEXT] == ED_CONTEXT_INTRONIC)
                    variantContext = "INTRONIC";
                else
                    variantContext = "EXONIC";

                mWriter.write(String.format(",%d,%d,%s",
                        exonMatchData[ED_RANK_START], exonMatchData[ED_RANK_END], variantContext));

                mWriter.write(String.format(",%s,%s,%s,%d,%d,%d,%d,%d",
                        altSJ.type(), altSJ.RegionContexts[SE_START], altSJ.RegionContexts[SE_END],
                        altSJ.SpliceJunction[SE_START], altSJ.SpliceJunction[SE_END],
                        altSJ.getFragmentCount(), altSJ.getPositionCount(SE_START), altSJ.getPositionCount(SE_END)));

                mWriter.newLine();
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write splice variant file: {}", e.toString());
        }
    }

    private void loadSpliceVariants(final String filename)
    {
        if(!Files.exists(Paths.get(filename)))
        {
            ISF_LOGGER.error("invalid splice variant file({})", filename);
            return;
        }

        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));

            if(mFieldsMap.isEmpty())
                mFieldsMap.putAll(createFieldsIndexMap(lines.get(0), DELIMITER));

            lines.remove(0);

            int sampleIdIndex = mFieldsMap.get("SampleId");

            List<SpliceVariant> spliceVariants = null;
            String currentSampleId = "";

            for(final String data : lines)
            {
                final String[] items = data.split(DELIMITER);
                final String sampleId = items[sampleIdIndex];

                if(!sampleId.equals(currentSampleId))
                {
                    currentSampleId = sampleId;
                    spliceVariants = Lists.newArrayList();
                    mSampleSpliceVariants.put(sampleId, spliceVariants);
                }

                spliceVariants.add(SpliceVariant.fromCsv(items, mFieldsMap));
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load splice variant data file({}): {}", filename, e.toString());
            return;
        }
    }

}

package com.hartwig.hmftools.isofox.novel.cohort;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionWithin;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunction.ASJ_TRANS_NONE;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunctionType.MIXED_TRANS;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunctionType.NOVEL_3_PRIME;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunctionType.NOVEL_5_PRIME;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunctionType.NOVEL_EXON;
import static com.hartwig.hmftools.isofox.novel.AltSpliceJunctionType.SKIPPED_EXONS;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunction;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunctionType;

public class SpliceVariantMatching
{
    private final CohortConfig mConfig;
    private final Map<String,Integer> mFieldsMap;
    private final EnsemblDataCache mGeneTransCache;

    private Map<String,List<SpliceVariant>> mSampleSpliceVariants;

    private BufferedWriter mWriter;

    private final Map<String,EnsemblGeneData> mGeneDataMap;

    public SpliceVariantMatching(final CohortConfig config, final String spliceVariantFile)
    {
        mConfig = config;
        mFieldsMap = Maps.newHashMap();
        mGeneDataMap = Maps.newHashMap();

        mGeneTransCache = new EnsemblDataCache(mConfig.EnsemblDataCache, RefGenomeVersion.HG37);
        mGeneTransCache.setRequiredData(true, false, false, false);
        mGeneTransCache.load(false);
        mSampleSpliceVariants = Maps.newHashMap();
        mWriter = null;

        initialiseWriter();

        if(spliceVariantFile != null)
            loadSpliceVariants(spliceVariantFile);

        for(String geneId : mConfig.RestrictedGeneIds)
        {
            final EnsemblGeneData geneData = mGeneTransCache.getGeneDataById(geneId);

            if(geneData != null)
            {
                mGeneDataMap.put(geneData.GeneName, geneData);
            }
        }
    }

    public void close() { closeBufferedWriter(mWriter); }

    public void evaluateSpliceVariants(final String sampleId, final List<AltSpliceJunction> altSpliceJunctions)
    {
        if(mGeneTransCache == null || (mSampleSpliceVariants.isEmpty() && mConfig.DbAccess == null))
            return;

        final List<SpliceVariant> spliceVariants = getSomaticVariants(sampleId);

        if(spliceVariants == null || spliceVariants.isEmpty())
            return;

        ISF_LOGGER.debug("sampleId({}) evaluating {} splice variants",
                sampleId, spliceVariants.size());

        spliceVariants.forEach(x -> evaluateSpliceVariant(sampleId, altSpliceJunctions, x));
    }

    private final List<SpliceVariant> getSomaticVariants(final String sampleId)
    {
        if(!mSampleSpliceVariants.isEmpty())
            return mSampleSpliceVariants.get(sampleId);

        final List<SpliceVariant> spliceVariants = Lists.newArrayList();

        final List<SomaticVariant> somaticVariants = mConfig.DbAccess.readSomaticVariants(sampleId, VariantType.UNDEFINED);

        // filter to specific gene list
        for(final SomaticVariant variant : somaticVariants)
        {
            if(!mGeneDataMap.isEmpty() && !mGeneDataMap.containsKey(variant.gene()))
                continue;

            spliceVariants.add(
                    new SpliceVariant(variant.gene(), variant.chromosome(), (int)variant.position(), variant.type(),
                            variant.ref(), variant.alt(), variant.worstEffectTranscript(), variant.worstEffect()));
        }

        return spliceVariants;
    }

    private void evaluateSpliceVariant(final String sampleId, final List<AltSpliceJunction> altSpliceJunctions, final SpliceVariant variant)
    {
        findCloseAltSpliceJunctions(sampleId, altSpliceJunctions, variant);
        findRelatedAltSpliceJunctions(sampleId, altSpliceJunctions, variant);
    }

    private static final int CLOSE_ALT_SJ_DISTANCE = 5;

    private void findCloseAltSpliceJunctions(
            final String sampleId, final List<AltSpliceJunction> altSpliceJunctions, final SpliceVariant variant)
    {
        final List<AltSpliceJunction> closeAltSJs = altSpliceJunctions.stream()
                .filter(x -> withinRange(x.SpliceJunction[SE_START], variant.Position, CLOSE_ALT_SJ_DISTANCE)
                        || withinRange(x.SpliceJunction[SE_END], variant.Position, CLOSE_ALT_SJ_DISTANCE))
                .collect(Collectors.toList());

        if(closeAltSJs.isEmpty())
            return;

        final EnsemblGeneData geneData = mGeneTransCache.getGeneDataByName(variant.GeneName);

        for(AltSpliceJunction altSJ : closeAltSJs)
        {
            List<String> allTransNames = Lists.newArrayList();

            for(int se = SE_START; se <= SE_END; ++se)
            {
                final String[] transNames = altSJ.getTranscriptNames()[se].split(ITEM_DELIM);
                for(String transName : transNames)
                {
                    if(!transName.equals("NONE") && !allTransNames.contains(transName))
                        allTransNames.add(transName);
                }
            }

            writeMatchData(sampleId, variant, geneData, altSJ, "PROXIMATE", appendStrList(allTransNames, ITEM_DELIM.charAt(0)));
        }
    }

    private static final int SPLICE_REGION_DISTANCE = 10;

    private static boolean withinRange(int pos1, int pos2, int distance) { return abs(pos1 - pos2) <= distance; }

    private static boolean validRelatedType(AltSpliceJunctionType type)
    {
        return type == SKIPPED_EXONS || type == NOVEL_5_PRIME || type == NOVEL_3_PRIME || type == NOVEL_EXON || type == MIXED_TRANS;
    }

    private void findRelatedAltSpliceJunctions(final String sampleId, final List<AltSpliceJunction> altSpliceJunctions, final SpliceVariant variant)
    {
        /* Candidate alt SJs must be either

        The splice region variant will be within 10 base of an exon.
        It may cause a novel splice junction anywhere in that exon or the intron.
        It could also cause the exon to be skipped entirely. eg. the splice variant is 7 bases after exon 3.
        Then you should check for novel splice junctions that are wholly contained within the entire region from end of exon 2 to start of exon 4.
        */

        final EnsemblGeneData geneData = mGeneDataMap.get(variant.GeneName);

        if(geneData == null)
        {
            ISF_LOGGER.error("variant({}) gene data not found");
            return;
        }

        final List<TranscriptData> transDataList = mGeneTransCache.getTranscripts(geneData.GeneId);

        for(final TranscriptData transData : transDataList)
        {
            for(int i = 0; i < transData.exons().size(); ++i)
            {
                final ExonData exon = transData.exons().get(i);

                if(!withinRange(exon.ExonStart, variant.Position, SPLICE_REGION_DISTANCE)
                && !withinRange(exon.ExonEnd, variant.Position, SPLICE_REGION_DISTANCE))
                {
                    continue;
                }

                final ExonData prevExon = i > 0 ? transData.exons().get(i - 1) : null;
                final ExonData nextExon = i < transData.exons().size() - 1 ? transData.exons().get(i + 1) : null;

                // look for any alt SJs

                for(final AltSpliceJunction altSJ : altSpliceJunctions)
                {
                    if(!altSJ.getGeneId().equals(geneData.GeneId))
                        continue;

                    if(!validRelatedType(altSJ.type()))
                        continue;

                    // check for a position within the exon boundaries
                    int minAsjPos = min(altSJ.SpliceJunction[SE_START], altSJ.SpliceJunction[SE_END]);
                    int maxAsjPos = max(altSJ.SpliceJunction[SE_START], altSJ.SpliceJunction[SE_END]);

                    if(prevExon != null && minAsjPos < prevExon.ExonEnd)
                        continue;

                    if(nextExon != null && maxAsjPos > nextExon.ExonEnd)
                        continue;

                    ISF_LOGGER.trace("sampleId({}) variant({}:{}) gene({}) transcript({}) matched altSJ({})",
                            sampleId, variant.Chromosome, variant.Position, geneData.GeneName, transData.TransName, altSJ.toString());

                    final String transDataStr = String.format("%s;%d", transData.TransName, exon.ExonRank);
                    writeMatchData(sampleId, variant, geneData, altSJ, "RELATED", transDataStr);
                }
            }
        }
    }

    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("splice_variant_matching.csv");
            mWriter = createBufferedWriter(outputFileName, false);

            mWriter.write("SampleId,MatchType,GeneId,GeneName,Chromosome,Strand");
            mWriter.write(",Position,Type,CodingEffect,Ref,Alt,SnpeffTransName");
            mWriter.write(",AsjType,AsjContextStart,AsjContextEnd,AsjPosStart,AsjPosEnd,AsjFragmentCount,AsjDepthStart,AsjDepthEnd");
            mWriter.write(",AsjTransData");

            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write splice variant file: {}", e.toString());
        }
    }

    private void writeMatchData(
            final String sampleId, final SpliceVariant variant, final EnsemblGeneData geneData,
            final AltSpliceJunction altSJ, final String matchType, final String transDataStr)
    {
        try
        {
            mWriter.write(String.format("%s,%s,%s,%s,%s,%d",
                    sampleId, matchType, geneData.GeneId, geneData.GeneName, geneData.Chromosome, geneData.Strand));

            mWriter.write(String.format(",%d,%s,%s,%s,%s,%s",
                    variant.Position, variant.Type, variant.CodingEffect, variant.Ref, variant.Alt, variant.TransName));

            mWriter.write(String.format(",%s,%s,%s,%d,%d,%d,%d,%d,%s",
                    altSJ.type(), altSJ.RegionContexts[SE_START], altSJ.RegionContexts[SE_END],
                    altSJ.SpliceJunction[SE_START], altSJ.SpliceJunction[SE_END],
                    altSJ.getFragmentCount(), altSJ.getPositionCount(SE_START), altSJ.getPositionCount(SE_END),
                    transDataStr));

            mWriter.newLine();
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

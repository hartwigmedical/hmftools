package com.hartwig.hmftools.isofox.data_loaders;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.common.RnaUtils.positionWithin;
import static com.hartwig.hmftools.isofox.data_loaders.DataLoader.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunction;

public class SpliceVariantMatching
{
    private final DataLoaderConfig mConfig;
    private final Map<String,Integer> mFieldsMap;
    private final EnsemblDataCache mGeneTransCache;

    private Map<String,List<SpliceVariant>> mSampleSpliceVariants;

    public SpliceVariantMatching(final DataLoaderConfig config)
    {
        mConfig = config;
        mFieldsMap = Maps.newHashMap();

        mGeneTransCache = new EnsemblDataCache(mConfig.EnsemblDataCache, RefGenomeVersion.HG37);
        mGeneTransCache.setRequiredData(true, false, false, false);
        mGeneTransCache.load(false);
        mSampleSpliceVariants = Maps.newHashMap();

        loadSpliceVariants(mConfig.SpliceVariantFile);
    }

    public void evaluateSpliceVariants(final String sampleId, final List<AltSpliceJunction> altSpliceJunctions)
    {
        if(mGeneTransCache == null || mSampleSpliceVariants.isEmpty())
            return;

        final List<SpliceVariant> spliceVariants = mSampleSpliceVariants.get(sampleId);

        if(spliceVariants == null || spliceVariants.isEmpty())
            return;

        spliceVariants.forEach(x -> evaluateSpliceVariant(sampleId, altSpliceJunctions, x));
    }

    private void evaluateSpliceVariant(final String sampleId, final List<AltSpliceJunction> altSpliceJunctions, final SpliceVariant variant)
    {
        // search the specified transcript for the next splice junction - alt or otherwise
        final EnsemblGeneData geneData = mGeneTransCache.getGeneDataByName(variant.GeneName);
        final TranscriptData transData = geneData != null ? mGeneTransCache.getTranscriptData(geneData.GeneId, variant.TransName) : null;

        if(geneData == null || transData == null)
        {
            ISF_LOGGER.warn("sampleId({}) variant({}:{}) cannot match gene({}) or transcript({})",
                    sampleId, variant.Chromosome, variant.Position, variant.GeneName, variant.TransName);
            return;
        }

        /*
        The splice region variant will be within 10 base of an exon.
        It may cause a novel splice junction anywhere in that exon or the intron.
        It could also cause the exon to be skipped entirely. eg. the splice variant is 7 bases after exon 3.
        Then you should check for novel splice junctions that are wholly contained within the entire region from end of exon 2 to start of exon 4.
        */

        final long[] evalRegion = { -1, -1 };

        for(int i = 0; i < transData.exons().size(); ++i)
        {
            final ExonData prevExon = i > 0 ? transData.exons().get(i - 1) : null;
            final ExonData exon = transData.exons().get(i);
            final ExonData nextExon = i < transData.exons().size() - 1 ? transData.exons().get(i + 1) : null;

            if(positionWithin(variant.Position, exon.ExonStart, exon.ExonEnd))
            {
                evalRegion[SE_START] = prevExon != null ? prevExon.ExonEnd : exon.ExonStart;
                evalRegion[SE_END] = nextExon != null ?  nextExon.ExonStart : exon.ExonEnd;
                break;
            }

            if(exon.ExonEnd < variant.Position && nextExon != null && nextExon.ExonStart > variant.Position)
            {
                evalRegion[SE_START] = prevExon != null ? prevExon.ExonEnd : exon.ExonStart;
                evalRegion[SE_END] = nextExon.ExonStart;
                break;
            }
        }

        if(evalRegion[SE_START] == -1 || evalRegion[SE_END] == -1)
        {
            ISF_LOGGER.warn("sampleId({}) variant({}:{}) gene({}) cannot match transcript({}) exons",
                    sampleId, variant.Chromosome, variant.Position, variant.GeneName, variant.TransName);
            return;
        }


        final List<AltSpliceJunction> candidateAltSJs = altSpliceJunctions.stream()
                .filter(x -> x.getGeneId().equals(geneData.GeneId))
                .filter(x -> positionWithin(x.SpliceJunction[SE_START], evalRegion[SE_START], evalRegion[SE_END]))
                .collect(Collectors.toList());

        if(candidateAltSJs.isEmpty())
        {
            ISF_LOGGER.debug("sampleId({}) variant({}:{}) gene({}) transcript({}) finds no candidate alt-SJs",
                    sampleId, variant.Chromosome, variant.Position, variant.GeneName, variant.TransName);
            return;
        }

        for(final AltSpliceJunction altSJ : candidateAltSJs)
        {
            ISF_LOGGER.info("sampleId({}) variant({}:{}) gene({}) transcript({}) matched altSJ({})",
                    sampleId, variant.Chromosome, variant.Position, variant.GeneName, variant.TransName, altSJ.toString());
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

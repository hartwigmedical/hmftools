package com.hartwig.hmftools.purple.targeted;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.TaggedRegion;

public class TargetRegionsData
{
    private final Map<String,List<TaggedRegion>> mTargetRegions;

    private int mTotalBases;
    private int mCodingBases;

    private double mTmlRatio;
    private double mTmbRatio;

    private int mCodingBaseFactor;

    private boolean mIsValid;

    public static final List<String> TMB_GENE_EXCLUSIONS = Lists.newArrayList("HLA-A", "HLA-B", "HLA-C", "PIM1", "BCL2");

    // target-region TML and TMB
    public static final int DEFAULT_CODING_BASE_FACTOR = 150000;
    public static final double PANEL_SOMATIC_LIKELIHOOD_DIFF_LOW = 0.08;
    public static final double PANEL_SOMATIC_LIKELIHOOD_DIFF_HIGH = -0.05;

    public TargetRegionsData(final String ratiosFile)
    {
        mTotalBases = 0;
        mCodingBases = 0;
        mTmlRatio = 1;
        mTmbRatio = 1;
        mCodingBaseFactor = DEFAULT_CODING_BASE_FACTOR;

        mTargetRegions = Maps.newHashMap();
        mIsValid = true;

        loadTargetRegionsRatios(ratiosFile);
    }

    public Map<String,List<TaggedRegion>> targetRegions() { return mTargetRegions; }
    public boolean hasTargetRegions()
    {
        return !mTargetRegions.isEmpty();
    }

    public boolean isValid()
    {
        return mIsValid;
    }

    public boolean inTargetRegions(final String chromosome, int position)
    {
        final List<TaggedRegion> chrRegions = mTargetRegions.get(chromosome);

        if(chrRegions == null)
        {
            return false;
        }

        return chrRegions.stream().anyMatch(x -> x.containsPosition(position));
    }

    public int codingBases()
    {
        return mCodingBases;
    }
    public double tmlRatio()
    {
        return mTmlRatio;
    }
    public double tmbRatio()
    {
        return mTmbRatio;
    }
    public int codingBaseFactor()
    {
        return mCodingBaseFactor;
    }

    public void loadTargetRegionsBed(final String targetRegionsBed, final EnsemblDataCache ensemblDataCache)
    {
        if(targetRegionsBed == null)
            return;

        Map<Chromosome, List<TaggedRegion>> chrRegionsMap = TaggedRegion.loadRegionsFromBedFile(targetRegionsBed);

        if(chrRegionsMap == null)
        {
            System.exit(1);
        }

        for(Map.Entry<Chromosome, List<TaggedRegion>> entry : chrRegionsMap.entrySet())
        {
            String chromosome = ensemblDataCache.refGenomeVersion().versionedChromosome(entry.getKey().toString());

            List<TaggedRegion> chrRegions = entry.getValue();

            mTargetRegions.put(chromosome, chrRegions);

            List<GeneData> geneList = ensemblDataCache.getChrGeneDataMap().get(chromosome);
            List<TranscriptData> overlappedTranscripts = Lists.newArrayList();

            // find the genes and then coding transcripts which overlap with these entries
            for(GeneData geneData : geneList)
            {
                if(TMB_GENE_EXCLUSIONS.contains(geneData.GeneName))
                {
                    continue;
                }

                if(chrRegions.stream().anyMatch(x -> positionsOverlap(x.start(), x.end(), geneData.GeneStart, geneData.GeneEnd)))
                {
                    TranscriptData transcriptData = ensemblDataCache.getCanonicalTranscriptData(geneData.GeneId);

                    if(transcriptData != null && !transcriptData.nonCoding())
                    {
                        overlappedTranscripts.add(transcriptData);
                    }
                }
            }

            // now compute the coding bases that overlap the target regions
            for(TaggedRegion region : chrRegions)
            {
                mTotalBases += region.baseLength();

                List<BaseRegion> exonicRegions = Lists.newArrayList();

                for(TranscriptData transcriptData : overlappedTranscripts)
                {
                    for(ExonData exon : transcriptData.exons())
                    {
                        if(transcriptData.CodingStart > exon.End)
                            continue;

                        if(exon.Start > transcriptData.CodingEnd)
                            break;

                        int maxCodingStart = max(transcriptData.CodingStart, exon.Start);
                        int minCodingEnd = min(transcriptData.CodingEnd, exon.End);

                        if(positionsOverlap(maxCodingStart, minCodingEnd, region.start(), region.end()))
                        {
                            exonicRegions.add(new BaseRegion(max(maxCodingStart, region.start()), min(minCodingEnd, region.end())));
                        }
                    }
                }

                BaseRegion.checkMergeOverlaps(exonicRegions);
                int regionCodingBases = exonicRegions.stream().mapToInt(x -> x.baseLength()).sum();

                mCodingBases += regionCodingBases;
            }
        }

        PPL_LOGGER.info("loaded {} target regions bases(total={} coding={}) from file({})",
                mTargetRegions.values().stream().mapToInt(List::size).sum(), mTotalBases, mCodingBases, targetRegionsBed);
    }

    private void loadTargetRegionsRatios(final String filename)
    {
        if(filename != null)
        {
            try
            {
                List<String> lines = Files.readAllLines(Paths.get(filename));

                String header = lines.get(0);
                Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
                String[] values = lines.get(1).split(TSV_DELIM, -1);

                mTmbRatio = Double.parseDouble(values[fieldsIndexMap.get("TmbRatio")]);
                mTmlRatio = Double.parseDouble(values[fieldsIndexMap.get("TmlRatio")]);
                if(fieldsIndexMap.containsKey("CodingBaseFactor"))
                {
                    mCodingBaseFactor = Integer.parseInt(values[fieldsIndexMap.get("CodingBaseFactor")]);
                }

                PPL_LOGGER.info("target regions: tml({}) tmb({}) codingBaseFactor({})",
                        mTmlRatio, mTmbRatio, mCodingBaseFactor);
            }
            catch(IOException e)
            {
                mIsValid = false;
                PPL_LOGGER.error("failed to load target regions ratios file: {}", e.toString());
            }
        }
    }
}

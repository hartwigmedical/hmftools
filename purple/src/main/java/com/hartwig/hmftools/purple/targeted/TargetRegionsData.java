package com.hartwig.hmftools.purple.targeted;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.purple.PurpleConstants.DEFAULT_CODING_BASE_FACTOR;
import static com.hartwig.hmftools.purple.PurpleConstants.DEFAULT_TARGETED_TMB_RATIO;
import static com.hartwig.hmftools.purple.PurpleConstants.DEFAULT_TARGETED_TML_RATIO;
import static com.hartwig.hmftools.purple.PurpleConstants.TARGETED_TMB_GENE_EXCLUSIONS;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

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
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class TargetRegionsData
{
    private final Map<String,List<TaggedRegion>> mTargetRegions;

    private int mTotalBases;
    private int mCodingBases;

    private int mCodingBaseFactor;
    private double mTmbRatio;
    private double mTmlRatio;

    private boolean mIsValid;

    private static final String TARGETED_TMB_RATIO = "target_tmb_ratio";
    private static final String TARGETED_TML_RATIO = "target_tml_ratio";
    private static final String TARGETED_CODING_FACTOR = "target_coding_factor";

    public TargetRegionsData(final int codingBaseFactor, final double tmbRatio, final double tmlRatio)
    {
        mTotalBases = 0;
        mCodingBases = 0;
        mTmlRatio = tmlRatio;
        mTmbRatio = tmbRatio;
        mCodingBaseFactor = codingBaseFactor;

        mTargetRegions = Maps.newHashMap();
        mIsValid = true;
    }

    public TargetRegionsData(final ConfigBuilder configBuilder)
    {
        this(configBuilder.getInteger(TARGETED_CODING_FACTOR),
                configBuilder.getDecimal(TARGETED_TMB_RATIO),
                configBuilder.getDecimal(TARGETED_TML_RATIO));
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addInteger(TARGETED_CODING_FACTOR, "Targeted panel coding base factor" ,DEFAULT_CODING_BASE_FACTOR);
        configBuilder.addDecimal(TARGETED_TMB_RATIO, "Targeted panel TMB adjustment factor", DEFAULT_TARGETED_TMB_RATIO);
        configBuilder.addDecimal(TARGETED_TML_RATIO, "Targeted panel TML adjustment factor", DEFAULT_TARGETED_TML_RATIO);
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
                if(TARGETED_TMB_GENE_EXCLUSIONS.contains(geneData.GeneName))
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
}

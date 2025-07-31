package com.hartwig.hmftools.purple;

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
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.TaggedRegion;

public class TargetRegionsData
{
    private final Map<String, List<TaggedRegion>> mTargetRegions;
    private final Map<String, List<Integer>> mTargetRegionsMsiIndels;

    private int mTotalBases;
    private int mCodingBases;

    private double mTmlRatio;
    private double mTmbRatio;
    private double mMsiIndelRatio;
    private double mMsi23BaseAF;
    private double mMsi4BaseAF;
    private int mCodingBaseFactor;

    private boolean mIsValid;

    public static final List<String> TMB_GENE_EXCLUSIONS = Lists.newArrayList("HLA-A", "HLA-B", "HLA-C", "PIM1", "BCL2");

    // target-region TML, TMB and MSI-Indels
    public static final double DEFAULT_MSI_2_3_BASE_AF = 0.15;
    public static final double DEFAULT_MSI_4_BASE_AF = 0.08;
    public static final int DEFAULT_CODING_BASE_FACTOR = 150000;
    public static final double PANEL_SOMATIC_LIKELIHOOD_DIFF_LOW = 0.08;
    public static final double PANEL_SOMATIC_LIKELIHOOD_DIFF_HIGH = -0.05;

    public TargetRegionsData(final String ratiosFile, final String msiIndelsFile)
    {
        mTotalBases = 0;
        mCodingBases = 0;
        mTmlRatio = 1;
        mTmbRatio = 1;
        mMsiIndelRatio = 1;
        mMsi23BaseAF = DEFAULT_MSI_2_3_BASE_AF;
        mMsi4BaseAF = DEFAULT_MSI_4_BASE_AF;
        mCodingBaseFactor = DEFAULT_CODING_BASE_FACTOR;

        mTargetRegions = Maps.newHashMap();
        mTargetRegionsMsiIndels = Maps.newHashMap();
        mIsValid = true;

        loadTargetRegionsMsiIndels(msiIndelsFile);
        loadTargetRegionsRatios(ratiosFile);
    }

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

    public boolean isTargetRegionsMsiIndel(final String chromsome, int position)
    {
        final List<Integer> chrRegions = mTargetRegionsMsiIndels.get(chromsome);

        if(chrRegions == null)
        {
            return false;
        }

        return chrRegions.stream().anyMatch(x -> position == x);
    }

    public int codingBases()
    {
        return mCodingBases;
    }

    public int msiIndelSiteCount()
    {
        return mTargetRegionsMsiIndels.values().stream().mapToInt(List::size).sum();
    }

    public double tmlRatio()
    {
        return mTmlRatio;
    }

    public double tmbRatio()
    {
        return mTmbRatio;
    }

    public double msiIndelRatio()
    {
        return mMsiIndelRatio;
    }

    public double msi23BaseAF()
    {
        return mMsi23BaseAF;
    }

    public double msi4BaseAF()
    {
        return mMsi4BaseAF;
    }

    public int codingBaseFactor()
    {
        return mCodingBaseFactor;
    }

    public List<TaggedRegion> targetRegions(String chromosome)
    {
        return mTargetRegions.get(chromosome);
    }

    public void loadTargetRegionsBed(final String targetRegionsBed, final EnsemblDataCache ensemblDataCache)
    {
        if(targetRegionsBed == null)
        {
            return;
        }

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

            // now compute the coding bases
            for(BaseRegion region : chrRegions)
            {
                mTotalBases += region.baseLength();

                int codingMinStart = -1;
                int codingMaxEnd = -1;

                for(TranscriptData transcriptData : overlappedTranscripts)
                {
                    if(positionsOverlap(transcriptData.CodingStart, transcriptData.CodingEnd, region.start(), region.end()))
                    {
                        int maxStart = max(transcriptData.CodingStart, region.start());
                        codingMinStart = codingMinStart > 0 ? min(codingMinStart, maxStart) : maxStart;

                        int minEnd = min(transcriptData.CodingEnd, region.end());
                        codingMaxEnd = codingMaxEnd > 0 ? max(codingMaxEnd, minEnd) : minEnd;

                        if(codingMinStart == region.start() && codingMaxEnd == region.end())
                        {
                            break;
                        }
                    }
                }

                if(codingMinStart > 0 && codingMaxEnd > 0)
                {
                    mCodingBases += codingMaxEnd - codingMinStart + 1;
                }
            }
        }

        PPL_LOGGER.info("loaded {} target regions bases(total={} coding={}) from file({})",
                mTargetRegions.values().stream().mapToInt(List::size).sum(), mTotalBases, mCodingBases, targetRegionsBed);
    }

    private void loadTargetRegionsMsiIndels(final String filename)
    {
        if(filename != null)
        {
            try
            {
                List<String> lines = Files.readAllLines(Paths.get(filename));

                String header = lines.get(0);
                lines.remove(0);
                Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
                int chrIndex = fieldsIndexMap.get("Chromosome");
                int posIndex = fieldsIndexMap.get("Position");

                for(String line : lines)
                {
                    String[] values = line.split(TSV_DELIM, -1);

                    String chromosome = values[chrIndex];
                    int position = Integer.parseInt(values[posIndex]);
                    List<Integer> positions = mTargetRegionsMsiIndels.computeIfAbsent(chromosome, k -> Lists.newArrayList());
                    positions.add(position);
                }

                PPL_LOGGER.info("loaded {} MSI indels from file({})",
                        mTargetRegionsMsiIndels.values().stream().mapToInt(List::size).sum(), filename);
            }
            catch(IOException e)
            {
                mIsValid = false;
                PPL_LOGGER.error("failed to load target regions ratios file: {}", e.toString());
            }
        }
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
                mMsiIndelRatio = Double.parseDouble(values[fieldsIndexMap.get("MsiIndelRatio")]);

                if(fieldsIndexMap.containsKey("Msi23BaseAF"))
                {
                    mMsi23BaseAF = Double.parseDouble(values[fieldsIndexMap.get("Msi23BaseAF")]);
                }

                if(fieldsIndexMap.containsKey("Msi4BaseAF"))
                {
                    mMsi4BaseAF = Double.parseDouble(values[fieldsIndexMap.get("Msi4BaseAF")]);
                }

                if(fieldsIndexMap.containsKey("CodingBaseFactor"))
                {
                    mCodingBaseFactor = Integer.parseInt(values[fieldsIndexMap.get("CodingBaseFactor")]);
                }

                PPL_LOGGER.info("target regions: tml({}) tmb({}) msiIndels({}) msiAF(2-3 base={} 4 base={}) codingBaseFactor({})",
                        mTmlRatio, mTmbRatio, mMsiIndelRatio, mMsi23BaseAF, mMsi4BaseAF, mCodingBaseFactor);
            }
            catch(IOException e)
            {
                mIsValid = false;
                PPL_LOGGER.error("failed to load target regions ratios file: {}", e.toString());
            }
        }
    }
}

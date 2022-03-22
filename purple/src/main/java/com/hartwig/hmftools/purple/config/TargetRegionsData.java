package com.hartwig.hmftools.purple.config;

import static java.lang.Math.round;

import static com.hartwig.hmftools.common.genome.bed.NamedBedFile.readBedFile;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;
import static com.hartwig.hmftools.purple.config.PurpleConstants.CODING_BASES_PER_GENOME;
import static com.hartwig.hmftools.purple.config.PurpleConstants.MB_PER_GENOME;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public class TargetRegionsData
{
    private final Map<Chromosome,List<NamedBed>> mTargetRegions;
    private int mTotalBases;
    private int mCodingBases;
    private double mTmlRatio;
    private double mTmbRatio;
    private double mMsiIndelRatio;

    private boolean mIsValid;

    private static final String CODING_REGION_ID = "CODING";

    public TargetRegionsData(final String bedFile, final String ratiosFile)
    {
        mTotalBases = 0;
        mCodingBases = 0;
        mTmlRatio = 0;
        mTmbRatio = 0;
        mMsiIndelRatio = 0;
        mTargetRegions = Maps.newHashMap();
        mIsValid = true;

        loadTargetRegionsBed(bedFile);
        loadTargetRegionsRatios(ratiosFile);

        if(!mTargetRegions.isEmpty())
        {
            PPL_LOGGER.info("loaded {} target regions from file({}) bases(total={} coding={}) ratios(tml={} tmb={} misIndels={})",
                    mTargetRegions.values().stream().mapToInt(x -> x.size()).sum(), bedFile,
                    mTotalBases, mCodingBases, mTmlRatio, mTmbRatio, mMsiIndelRatio);
        }
    }

    public boolean hasTargetRegions() { return !mTargetRegions.isEmpty(); }
    public boolean isValid() { return mIsValid; }

    public boolean inTargetRegions(final String chromsome, int position)
    {
        final List<NamedBed> chrRegions = mTargetRegions.get(HumanChromosome.fromString(chromsome));

        if(chrRegions == null)
            return false;

        return chrRegions.stream().anyMatch(x -> positionWithin(position, x.start(), x.end()));
    }

    public int calcTml(int rawTml)
    {
        if(mTargetRegions.isEmpty())
            return rawTml;

        // # of missense in targetedRegions * RefGenomeCodingBases / TargetRegionCodingBases
        return (int)round(rawTml * CODING_BASES_PER_GENOME / mCodingBases);
    }

    public double calcTmb(int rawTmb)
    {
        if(mTargetRegions.isEmpty())
            return rawTmb / MB_PER_GENOME;

        // # of variants in targetedRegions * RefGenomeSize / TargetRegionSize * Constant
        return rawTmb * MB_PER_GENOME * 1e6 / mTotalBases * mTmbRatio;
    }

    public double calcMsiIndels(int rawCount)
    {
        if(mTargetRegions.isEmpty())
            return (double) rawCount / MB_PER_GENOME;

        // # of MSI indels in targetedRegions  * RefGenomeSize / TargetRegionSize * Constant
        return rawCount * MB_PER_GENOME * 1e6 / mTotalBases * mMsiIndelRatio;
    }

    private void loadTargetRegionsBed(final String bedFile)
    {
        if(bedFile == null)
            return;

        try
        {
            List<NamedBed> namedBedRecords = readBedFile(bedFile);
            int totalRegionSize = 0;

            for(NamedBed namedBed : namedBedRecords)
            {
                HumanChromosome chromosome = HumanChromosome.fromString(namedBed.chromosome());
                List<NamedBed> chrRegions = mTargetRegions.get(chromosome);

                if(chrRegions == null)
                {
                    chrRegions = Lists.newArrayList();
                    mTargetRegions.put(chromosome, chrRegions);
                }

                chrRegions.add(namedBed);
                mTotalBases += namedBed.bases();

                if(namedBed.name().contains(CODING_REGION_ID))
                    mCodingBases += namedBed.bases();
            }
        }
        catch (IOException e)
        {
            mIsValid = false;
            PPL_LOGGER.error("failed to load target regions BED file: {}", e.toString());
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
                Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, ",");
                String[] values = lines.get(1).split(",", -1);

                mTmbRatio = Double.parseDouble(values[fieldsIndexMap.get("TmbRatio")]);
                mTmlRatio = Double.parseDouble(values[fieldsIndexMap.get("TmlRatio")]);
                mMsiIndelRatio = Double.parseDouble(values[fieldsIndexMap.get("MsiIndelRatio")]);
            }
            catch(IOException e)
            {
                mIsValid = false;
                PPL_LOGGER.error("failed to load target regions ratios file: {}", e.toString());
            }
        }
    }
}

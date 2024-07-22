package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.genome.bed.NamedBedFile.readBedFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.bed.NamedBed;

public class TargetRegionsData
{
    private final Map<String,List<NamedBed>> mTargetRegions;
    private final Map<String,List<Integer>> mTargetRegionsMsiIndels;

    private int mTotalBases;
    private int mCodingBases;

    private double mTmlRatio;
    private double mTmbRatio;
    private double mMsiIndelRatio;
    private double mMsi23BaseAF;
    private double mMsi4BaseAF;
    private int mCodingBaseFactor;

    private boolean mIsValid;

    private static final String CODING_REGION_ID = "CODING";

    public static final List<String> TMB_GENE_EXCLUSIONS = Lists.newArrayList("HLA-A","HLA-B","HLA-C","PIM1","BCL2");

    // target-region TML, TMB and MSI-Indels
    public static final double DEFAULT_MSI_2_3_BASE_AF = 0.15;
    public static final double DEFAULT_MSI_4_BASE_AF = 0.08;
    public static final int DEFAULT_CODING_BASE_FACTOR = 150000;
    public static final double PANEL_SOMATIC_LIKELIHOOD_DIFF_LOW = 0.08;
    public static final double PANEL_SOMATIC_LIKELIHOOD_DIFF_HIGH = -0.05;

    public TargetRegionsData(final String bedFile, final String ratiosFile, final String msiIndelsFile)
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

        loadTargetRegionsBed(bedFile);
        loadTargetRegionsMsiIndels(msiIndelsFile);
        loadTargetRegionsRatios(ratiosFile);
    }

    public boolean hasTargetRegions() { return !mTargetRegions.isEmpty(); }
    public boolean isValid() { return mIsValid; }

    public boolean inTargetRegions(final String chromsome, int position)
    {
        final List<NamedBed> chrRegions = mTargetRegions.get(chromsome);

        if(chrRegions == null)
            return false;

        return chrRegions.stream().anyMatch(x -> positionWithin(position, x.start(), x.end()));
    }

    public boolean isTargetRegionsMsiIndel(final String chromsome, int position)
    {
        final List<Integer> chrRegions = mTargetRegionsMsiIndels.get(chromsome);

        if(chrRegions == null)
            return false;

        return chrRegions.stream().anyMatch(x -> position == x);
    }

    public int codingBases() { return mCodingBases; }
    public int msiIndelSiteCount() { return mTargetRegionsMsiIndels.values().stream().mapToInt(x -> x.size()).sum(); }
    public double tmlRatio() { return mTmlRatio; }
    public double tmbRatio() { return mTmbRatio; }
    public double msiIndelRatio() { return mMsiIndelRatio; }
    public double msi23BaseAF() { return mMsi23BaseAF; }
    public double msi4BaseAF() { return mMsi4BaseAF; }
    public int codingBaseFactor() { return mCodingBaseFactor; }

    private void loadTargetRegionsBed(final String bedFile)
    {
        if(bedFile == null)
            return;

        try
        {
            List<NamedBed> namedBedRecords = readBedFile(bedFile);

            for(NamedBed namedBed : namedBedRecords)
            {
                List<NamedBed> chrRegions = mTargetRegions.get(namedBed.chromosome());

                if(chrRegions == null)
                {
                    chrRegions = Lists.newArrayList();
                    mTargetRegions.put(namedBed.chromosome(), chrRegions);
                }

                chrRegions.add(namedBed);
                mTotalBases += namedBed.bases();

                if(namedBed.name().contains(CODING_REGION_ID))
                    mCodingBases += namedBed.bases();
            }

            PPL_LOGGER.info("loaded {} target regions bases(total={} coding={}) from file({})",
                    mTargetRegions.values().stream().mapToInt(x -> x.size()).sum(), mTotalBases, mCodingBases, bedFile);
        }
        catch (IOException e)
        {
            mIsValid = false;
            PPL_LOGGER.error("failed to load target regions BED file: {}", e.toString());
        }
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
                Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
                int chrIndex = fieldsIndexMap.get("Chromosome");
                int posIndex = fieldsIndexMap.get("Position");

                for(String line : lines)
                {
                    String[] values = line.split(TSV_DELIM, -1);

                    String chromosome = values[chrIndex];
                    int position = Integer.parseInt(values[posIndex]);

                    List<Integer> positions = mTargetRegionsMsiIndels.get(chromosome);

                    if(positions == null)
                    {
                        positions = Lists.newArrayList();
                        mTargetRegionsMsiIndels.put(chromosome, positions);
                    }

                    positions.add(position);
                }

                PPL_LOGGER.info("loaded {} MSI INDELs from file({})",
                        mTargetRegionsMsiIndels.values().stream().mapToInt(x -> x.size()).sum(), filename);
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
                Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
                String[] values = lines.get(1).split(TSV_DELIM, -1);

                mTmbRatio = Double.parseDouble(values[fieldsIndexMap.get("TmbRatio")]);
                mTmlRatio = Double.parseDouble(values[fieldsIndexMap.get("TmlRatio")]);
                mMsiIndelRatio = Double.parseDouble(values[fieldsIndexMap.get("MsiIndelRatio")]);

                if(fieldsIndexMap.containsKey("Msi23BaseAF"))
                    mMsi23BaseAF = Double.parseDouble(values[fieldsIndexMap.get("Msi23BaseAF")]);

                if(fieldsIndexMap.containsKey("Msi4BaseAF"))
                    mMsi4BaseAF = Double.parseDouble(values[fieldsIndexMap.get("Msi4BaseAF")]);

                if(fieldsIndexMap.containsKey("CodingBaseFactor"))
                    mCodingBaseFactor = Integer.parseInt(values[fieldsIndexMap.get("CodingBaseFactor")]);

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

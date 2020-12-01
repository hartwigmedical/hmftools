package com.hartwig.hmftools.isofox.novel.cohort;

import static com.hartwig.hmftools.common.sigs.DataUtils.convertList;
import static com.hartwig.hmftools.common.sigs.Percentiles.calcPercentileValues;
import static com.hartwig.hmftools.common.sigs.Percentiles.getPercentile;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.SPLICE_SITE_PERCENTILES;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.expression.cohort.TransExpressionDistribution.DISTRIBUTION_SIZE;
import static com.hartwig.hmftools.isofox.novel.cohort.SpliceVariantMatcher.COHORT_ALT_SJ_FILE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunctionType;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class SpliceSiteCache
{
    private final CohortConfig mConfig;

    // for building a cache
    private final Map<String, Set<Integer>> mCohortAltSJPositions;
    private final Map<String,Map<Integer,List<Double>>> mSpliceSiteRates;

    // when used as a cache
    private final Map<String,Map<Integer,double[]>> mSpliceSitePercentilesMap;
    private final Map<String,Map<Integer,Double>> mCurrentSampleSpliceSiteMap;

    protected static final String COHORT_SPLICE_SITE_FILE = "cohort_splice_site_file";

    public SpliceSiteCache(final CohortConfig config, final CommandLine cmd)
    {
        mConfig = config;
        mCohortAltSJPositions = Maps.newHashMap();
        mSpliceSiteRates = Maps.newHashMap();

        if(cmd.hasOption(COHORT_ALT_SJ_FILE) && config.AnalysisTypes.contains(SPLICE_SITE_PERCENTILES))
        {
            loadCohortAltSJs(cmd.getOptionValue(COHORT_ALT_SJ_FILE));
        }

        mSpliceSitePercentilesMap = Maps.newHashMap();
        mCurrentSampleSpliceSiteMap = Maps.newHashMap();

        if(cmd.hasOption(COHORT_SPLICE_SITE_FILE))
        {
            loadSpliceSitePercentiles(cmd.getOptionValue(COHORT_SPLICE_SITE_FILE));
        }
    }

    public static void addCmdLineOptions(final Options options)
    {
        options.addOption(COHORT_ALT_SJ_FILE, true, "Cohort frequency for alt SJs");
        options.addOption(COHORT_SPLICE_SITE_FILE, true, "Cohort splice site percentiles");
    }

    public void createPercentiles()
    {
        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, SPLICE_SITE_PERCENTILES, filenames))
            return;

        final Map<String,Integer> fieldsIndexMap = Maps.newHashMap();

        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path sampleFile = filenames.get(i);

            final Map<String,Map<Integer,Double>> spliceSiteMap = loadSpliceSiteFile(sampleFile, fieldsIndexMap);

            ISF_LOGGER.debug("{}: sample({}) loaded {} splice-site records",
                    i, sampleId, spliceSiteMap.values().stream().mapToInt(x -> x.size()).sum());

            for(Map.Entry<String,Map<Integer,Double>> chrEntry : spliceSiteMap.entrySet())
            {
                final String chromosome = chrEntry.getKey();

                for(Map.Entry<Integer,Double> posEntry : chrEntry.getValue().entrySet())
                {
                    int splicePosition = posEntry.getKey();

                    if(!hasAltSpliceSitePosition(chromosome, splicePosition))
                        continue;

                    double supportRate = posEntry.getValue();

                    addSpliceSiteSupport(chromosome, splicePosition, supportRate);
                }
            }
        }

        ISF_LOGGER.info("loaded {} sample splice-site files", mConfig.SampleData.SampleIds.size());

        // write a cohort file
        writePercentiles();
    }

    public static final Map<String,Map<Integer,Double>> loadSpliceSiteFile(final Path filename, final Map<String,Integer> fieldsIndexMap)
    {
        final Map<String,Map<Integer,Double>> spliceSiteMap = Maps.newHashMap();

        try
        {
            final List<String> lines = Files.readAllLines(filename);

            if(fieldsIndexMap.isEmpty())
                fieldsIndexMap.putAll(createFieldsIndexMap(lines.get(0), DELIMITER));

            int chromosomeIndex = fieldsIndexMap.get("Chromosome");
            int splicePosIndex = fieldsIndexMap.get("SpliceSitePosition");
            int traverseFragsIndex = fieldsIndexMap.get("TraverseFrags");
            int supportFragsIndex = fieldsIndexMap.get("SupportFrags");

            lines.remove(0);

            int spliceSiteCount = 0;

            for(final String line : lines)
            {
                // GeneSetId,Chromosome,SpliceSitePosition,TraverseFrags,SupportFrags,SkipFrags
                final String[] items = line.split(DELIMITER);
                String chromosome = items[chromosomeIndex];
                int splicePosition = Integer.parseInt(items[splicePosIndex]);

                int traverseFrags = Integer.parseInt(items[traverseFragsIndex]);
                int supportFrags = Integer.parseInt(items[supportFragsIndex]);

                if(traverseFrags == 0 && supportFrags == 0)
                    continue;

                double supportRate = supportFrags / (double)(traverseFrags + supportFrags);

                Map<Integer,Double> positionsMap = spliceSiteMap.get(chromosome);

                if(positionsMap == null)
                {
                    positionsMap = Maps.newHashMap();
                    spliceSiteMap.put(chromosome, positionsMap);
                }

                positionsMap.put(splicePosition, supportRate);
                ++spliceSiteCount;
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load splice site file({}): {}", filename.toString(), e.toString());
        }

        return spliceSiteMap;
    }

    private void addSpliceSiteSupport(final String chromosome, int splicePosition, double supportRate)
    {
        Map<Integer,List<Double>> positionsMap = mSpliceSiteRates.get(chromosome);

        if(positionsMap == null)
        {
            positionsMap = Maps.newHashMap();
            mSpliceSiteRates.put(chromosome, positionsMap);
        }

        List<Double> supportRates = positionsMap.get(splicePosition);

        if(supportRates == null)
        {
            positionsMap.put(splicePosition, Lists.newArrayList(supportRate));
            return;
        }

        int index = 0;

        while(index < supportRates.size())
        {
            if(supportRate < supportRates.get(index))
                break;

            ++index;
        }

        supportRates.add(index, supportRate);
    }

    private boolean processType(AltSpliceJunctionType type)
    {
        switch(type)
        {
            case SKIPPED_EXONS:
            case NOVEL_3_PRIME:
            case NOVEL_5_PRIME:
            case NOVEL_EXON:
            case MIXED_TRANS:
            case EXON_INTRON:
                return true;
            default:
                return false;
        }
    }

    private void loadCohortAltSJs(final String filename)
    {
        if(!Files.exists(Paths.get(filename)))
        {
            ISF_LOGGER.error("invalid cohort alt-SJ file({})", filename);
            return;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null)
                return;

            final Map<String,Integer> fieldsMap = createFieldsIndexMap(line, DELIMITER);

            int typeIndex = fieldsMap.get("Type");
            int chromosomeIndex = fieldsMap.get("Chromosome");
            int sjStartPosIndex = fieldsMap.get("SjStart");
            int sjEndPosIndex = fieldsMap.get("SjEnd");

            int asjCount = 0;

            while ((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(DELIMITER);

                AltSpliceJunctionType type = AltSpliceJunctionType.valueOf(items[typeIndex]);

                if(!processType(type))
                    continue;

                String chromosome = items[chromosomeIndex];
                int asjPosStart = Integer.parseInt(items[sjStartPosIndex]);
                int asjPosEnd = Integer.parseInt(items[sjEndPosIndex]);

                Set<Integer> positions = mCohortAltSJPositions.get(chromosome);

                if(positions == null)
                {
                    positions = Sets.newHashSet();
                    mCohortAltSJPositions.put(chromosome, positions);
                }

                positions.add(asjPosStart);
                positions.add(asjPosEnd);

                ++asjCount;
            }

            ISF_LOGGER.info("loaded {} cohort alt-SJ records", asjCount);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load cohort alt-SJ data file({}): {}", filename, e.toString());
            return;
        }
    }

    private boolean hasAltSpliceSitePosition(final String chromosome, int splicePosition)
    {
        final Set<Integer> positions = mCohortAltSJPositions.get(chromosome);

        if(positions == null)
            return false;

        return positions.contains(splicePosition);
    }

    private void writePercentiles()
    {
        final double[] percentileValues = new double[DISTRIBUTION_SIZE];

        try
        {
            final String outputFileName = mConfig.formCohortFilename("splice_site_support_perc.csv");
            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("Chromosome,SplicePosition");

            for(int i = 0; i < DISTRIBUTION_SIZE; ++i)
            {
                writer.write(String.format(",Pct_%d", i));
            }

            writer.newLine();

            for(Map.Entry<String,Map<Integer,List<Double>>> chrEntry : mSpliceSiteRates.entrySet())
            {
                final String chromosome = chrEntry.getKey();

                for(Map.Entry<Integer,List<Double>> posEntry : chrEntry.getValue().entrySet())
                {
                    final int splicePosition = posEntry.getKey();

                    final double[] supportRates = convertList(posEntry.getValue());
                    calcPercentileValues(supportRates, percentileValues);

                    writer.write(String.format("%s,%d", chromosome, splicePosition));

                    for(int i = 0; i < DISTRIBUTION_SIZE; ++i)
                    {
                        writer.write(String.format(",%.3f", percentileValues[i]));
                    }

                    writer.newLine();
                }
            }

            closeBufferedWriter(writer);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write splice site percentiles file: {}", e.toString());
        }
    }

    private void loadSpliceSitePercentiles(final String filename)
    {
        if(!Files.exists(Paths.get(filename)))
        {
            ISF_LOGGER.error("invalid cohort splice site percentiles file({})", filename);
            return;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            // skip field names
            String line = fileReader.readLine();

            if (line == null)
            {
                ISF_LOGGER.error("empty cohort splice site percentiles file({})", filename);
                return;
            }

            final Map<String,Integer> fieldsMap = createFieldsIndexMap(line, DELIMITER);

            int chromosomeIndex = fieldsMap.get("Chromosome");
            int positionIndex = fieldsMap.get("SplicePosition");
            String currentChromosome = "";
            Map<Integer,double[]> postionsMap = null;

            int expectedColCount = 2 + DISTRIBUTION_SIZE;
            int splicePosCount = 0;

            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(DELIMITER, -1);

                if (items.length != expectedColCount)
                {
                    ISF_LOGGER.error("invalid splice site percentile data length({}) vs expected({}): {}",
                            items.length, expectedColCount, line);
                    return;
                }

                final String chromosome = items[chromosomeIndex];
                final int splicePosition = Integer.parseInt(items[positionIndex]);

                double[] percentileData = new double[DISTRIBUTION_SIZE];

                int startIndex = 2;
                for(int i = startIndex; i < items.length; ++i)
                {
                    double supportRate = Double.parseDouble(items[i]);
                    percentileData[i - startIndex] = supportRate;
                }

                if(!currentChromosome.equals(chromosome))
                {
                    currentChromosome = chromosome;
                    postionsMap = Maps.newHashMap();
                    mSpliceSitePercentilesMap.put(chromosome, postionsMap);
                }

                postionsMap.put(splicePosition, percentileData);
                ++splicePosCount;
            }

            ISF_LOGGER.info("loaded {} splice-position percentiles from file({})", splicePosCount, filename);
        }
        catch (IOException e)
        {
            ISF_LOGGER.warn("failed to load splice-position percentiles file({}): {}", filename, e.toString());
        }
    }

    public void loadSampleSpliceSites(final String sampleId)
    {
        mCurrentSampleSpliceSiteMap.clear();

        if(mSpliceSitePercentilesMap.isEmpty())
            return;

        final Map<String,Integer> ssFieldsIndexMap = Maps.newHashMap();
        final Path ssFilename = Paths.get(CohortConfig.formSampleFilename(mConfig, sampleId, SPLICE_SITE_PERCENTILES));
        mCurrentSampleSpliceSiteMap.putAll(SpliceSiteCache.loadSpliceSiteFile(ssFilename, ssFieldsIndexMap));
    }

    public static final double PSI_NO_RATE = -1;

    public double getSampleSpliceSitePsi(final String chromosome, int splicePosition)
    {
        final Map<Integer,Double> positions = mCurrentSampleSpliceSiteMap.get(chromosome);

        if(positions == null)
            return PSI_NO_RATE;

        Double rate = positions.get(splicePosition);

        return rate != null ? rate : PSI_NO_RATE;
    }

    public double getCohortSpliceSitePsi(final String chromosome, int splicePosition, double sampleRate)
    {
        final Map<Integer,double[]> percentilesMap = mSpliceSitePercentilesMap.get(chromosome);

        if(percentilesMap == null)
            return PSI_NO_RATE;

        final double[] percentiles = percentilesMap.get(splicePosition);

        if(percentiles == null)
            return PSI_NO_RATE;

        return getPercentile(percentiles, sampleRate);
    }

}

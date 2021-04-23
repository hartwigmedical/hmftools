package com.hartwig.hmftools.isofox.novel.cohort;

import static com.hartwig.hmftools.common.sigs.DataUtils.convertList;
import static com.hartwig.hmftools.common.stats.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.stats.Percentiles.calcPercentileValues;
import static com.hartwig.hmftools.common.stats.Percentiles.getPercentile;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.SPLICE_SITE_PERCENTILES;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
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

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class SpliceSiteCache
{
    private final CohortConfig mConfig;

    // for building a cache
    private final Map<String,Map<Integer,List<Double>>> mSpliceSiteRates;
    private final Map<String,Map<Integer,int[]>> mSpliceSiteTotals;

    // when used as a cache
    private final Map<String,Map<Integer,double[]>> mCohortSpliceSitePercentilesMap;
    private final Map<String,Map<Integer,int[]>> mCohortSpliceSiteTotalsMap;
    private final Map<String,Map<Integer,int[]>> mSampleSpliceSiteMap;

    private final EnsemblDataCache mGeneTransCache;

    protected static final String COHORT_SPLICE_SITE_FILE = "cohort_splice_site_file";

    public static final int SS_TRAVERSED = 0;
    public static final int SS_SUPPORT = 1;

    public SpliceSiteCache(final CohortConfig config, final CommandLine cmd)
    {
        mConfig = config;
        mSpliceSiteRates = Maps.newHashMap();
        mSpliceSiteTotals = Maps.newHashMap();

        mCohortSpliceSitePercentilesMap = Maps.newHashMap();
        mCohortSpliceSiteTotalsMap = Maps.newHashMap();
        mSampleSpliceSiteMap = Maps.newHashMap();

        if(config.EnsemblDataCache != null)
        {
            mGeneTransCache = new EnsemblDataCache(config.EnsemblDataCache, RefGenomeVersion.V37);
            mGeneTransCache.setRequiredData(false, false, false, false);
            mGeneTransCache.setRestrictedGeneIdList(config.RestrictedGeneIds);
            mGeneTransCache.load(true); // only need  genes, not transcript data
        }
        else
        {
            mGeneTransCache = null;
        }

        if(cmd.hasOption(COHORT_SPLICE_SITE_FILE))
        {
            loadSpliceSitePercentiles(cmd.getOptionValue(COHORT_SPLICE_SITE_FILE));
        }
    }

    public static void addCmdLineOptions(final Options options)
    {
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

            final Map<String,Map<Integer,int[]>> spliceSiteMap = loadSpliceSiteFile(sampleFile, fieldsIndexMap);

            ISF_LOGGER.debug("{}: sample({}) loaded {} splice-site records",
                    i, sampleId, spliceSiteMap.values().stream().mapToInt(x -> x.size()).sum());

            for(Map.Entry<String,Map<Integer,int[]>> chrEntry : spliceSiteMap.entrySet())
            {
                final String chromosome = chrEntry.getKey();

                final List<EnsemblGeneData> geneDataList = mGeneTransCache != null ?
                        mGeneTransCache.getChrGeneDataMap().get(chromosome) : null;

                for(Map.Entry<Integer,int[]> posEntry : chrEntry.getValue().entrySet())
                {
                    int splicePosition = posEntry.getKey();

                    if(mGeneTransCache != null && !inRequiredGenes(geneDataList, splicePosition))
                        continue;

                    final int[] spliceSiteData = posEntry.getValue();

                    addSpliceSiteSupport(chromosome, splicePosition, spliceSiteData);
                }
            }
        }

        ISF_LOGGER.info("loaded {} sample splice-site files", mConfig.SampleData.SampleIds.size());

        // write a cohort file
        writePercentiles();
    }

    public static Map<String,Map<Integer,int[]>> loadSpliceSiteFile(final Path filename, final Map<String,Integer> fieldsIndexMap)
    {
        final Map<String,Map<Integer,int[]>> spliceSiteMap = Maps.newHashMap();

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

                final int[] spliceSiteData = new int[2];
                spliceSiteData[SS_SUPPORT] = supportFrags;
                spliceSiteData[SS_TRAVERSED] = traverseFrags;

                Map<Integer,int[]> positionsMap = spliceSiteMap.get(chromosome);

                if(positionsMap == null)
                {
                    positionsMap = Maps.newHashMap();
                    spliceSiteMap.put(chromosome, positionsMap);
                }

                positionsMap.put(splicePosition, spliceSiteData);
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load splice site file({}): {}", filename.toString(), e.toString());
        }

        return spliceSiteMap;
    }

    private static boolean inRequiredGenes(final List<EnsemblGeneData> geneDataList, int position)
    {
        if(geneDataList == null)
            return false;

        for(final EnsemblGeneData geneData : geneDataList)
        {
            if(positionWithin(position, geneData.GeneStart, geneData.GeneEnd))
                return true;
        }

        return false;
    }

    public static double calcSupportRate(final int[] spliceSiteData)
    {
        if(spliceSiteData == null || spliceSiteData.length != 2)
            return -1;

        double denom = spliceSiteData[SS_SUPPORT] + spliceSiteData[SS_TRAVERSED];

        if(denom == 0)
            return -1;

        return spliceSiteData[SS_SUPPORT] / denom;
    }

    private void addSpliceSiteSupport(final String chromosome, int splicePosition, final int[] spliceSiteData)
    {
        Map<Integer,List<Double>> positionsMap = mSpliceSiteRates.get(chromosome);

        if(positionsMap == null)
        {
            positionsMap = Maps.newHashMap();
            mSpliceSiteRates.put(chromosome, positionsMap);
        }

        double supportRate = calcSupportRate(spliceSiteData);

        List<Double> supportRates = positionsMap.get(splicePosition);

        if(supportRates == null)
        {
            positionsMap.put(splicePosition, Lists.newArrayList(supportRate));
        }
        else
        {
            int index = 0;

            while(index < supportRates.size())
            {
                if(supportRate < supportRates.get(index))
                    break;

                ++index;
            }

            supportRates.add(index, supportRate);
        }

        Map<Integer,int[]> positionTotalsMap = mSpliceSiteTotals.get(chromosome);

        if(positionTotalsMap == null)
        {
            positionTotalsMap = Maps.newHashMap();
            mSpliceSiteTotals.put(chromosome, positionTotalsMap);
        }

        int[] siteTotals = positionTotalsMap.get(splicePosition);

        if(siteTotals == null)
        {
            positionTotalsMap.put(splicePosition, spliceSiteData);
        }
        else
        {
            siteTotals[SS_SUPPORT] += spliceSiteData[SS_SUPPORT];
            siteTotals[SS_TRAVERSED] += spliceSiteData[SS_TRAVERSED];
        }
    }

    private void writePercentiles()
    {
        final double[] percentileValues = new double[PERCENTILE_COUNT];

        try
        {
            final String outputFileName = mConfig.formCohortFilename("splice_site_support_perc.csv");
            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("Chromosome,SplicePosition,FragTotal,SupportTotal");

            for(int i = 0; i < PERCENTILE_COUNT; ++i)
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

                    final int[] spliceSiteTotals = mSpliceSiteTotals.get(chromosome).get(splicePosition);
                    int totalSupport = spliceSiteTotals[SS_TRAVERSED] + spliceSiteTotals[SS_SUPPORT];

                    writer.write(String.format("%s,%d,%d,%d",
                            chromosome, splicePosition, totalSupport, spliceSiteTotals[SS_SUPPORT]));

                    for(int i = 0; i < PERCENTILE_COUNT; ++i)
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
            int fragsIndex = fieldsMap.get("FragTotal");
            int supportIndex = fieldsMap.get("SupportTotal");
            String currentChromosome = "";
            Map<Integer,double[]> percentilesMap = null;
            Map<Integer,int[]> totalsMap = null;

            int expectedColCount = 4 + PERCENTILE_COUNT;
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

                int totalFrags = Integer.parseInt(items[fragsIndex]);
                int supportFrags = Integer.parseInt(items[supportIndex]);

                final int[] totals = new int[2];
                totals[SS_TRAVERSED] = totalFrags - supportFrags;
                totals[SS_SUPPORT] = supportFrags;

                double[] percentileData = new double[PERCENTILE_COUNT];

                int startIndex = 4;
                for(int i = startIndex; i < items.length; ++i)
                {
                    double supportRate = Double.parseDouble(items[i]);
                    percentileData[i - startIndex] = supportRate;
                }

                if(!currentChromosome.equals(chromosome))
                {
                    currentChromosome = chromosome;
                    percentilesMap = Maps.newHashMap();
                    totalsMap = Maps.newHashMap();
                    mCohortSpliceSitePercentilesMap.put(chromosome, percentilesMap);
                    mCohortSpliceSiteTotalsMap.put(chromosome, totalsMap);
                }

                percentilesMap.put(splicePosition, percentileData);
                totalsMap.put(splicePosition, totals);
                ++splicePosCount;
            }

            ISF_LOGGER.info("loaded {} cohort splice-positions from file({})", splicePosCount, filename);
        }
        catch (IOException e)
        {
            ISF_LOGGER.warn("failed to load splice-position percentiles file({}): {}", filename, e.toString());
        }
    }

    public void loadSampleSpliceSites(final String sampleId)
    {
        mSampleSpliceSiteMap.clear();

        if(mCohortSpliceSitePercentilesMap.isEmpty())
            return;

        final Map<String,Integer> ssFieldsIndexMap = Maps.newHashMap();
        final Path ssFilename = Paths.get(CohortConfig.formSampleFilename(mConfig, sampleId, SPLICE_SITE_PERCENTILES));
        mSampleSpliceSiteMap.putAll(SpliceSiteCache.loadSpliceSiteFile(ssFilename, ssFieldsIndexMap));
    }

    public static final double PSI_NO_RATE = -1;

    public int[] getSampleSpliceSiteData(final String chromosome, int splicePosition)
    {
        final Map<Integer,int[]> positions = mSampleSpliceSiteMap.get(chromosome);

        if(positions == null)
            return null;

        return positions.get(splicePosition);
    }

    public double getCohortSpliceSitePsiPercentile(final String chromosome, int splicePosition, double sampleRate)
    {
        final Map<Integer,double[]> percentilesMap = mCohortSpliceSitePercentilesMap.get(chromosome);

        if(percentilesMap == null)
            return PSI_NO_RATE;

        final double[] percentiles = percentilesMap.get(splicePosition);

        if(percentiles == null)
            return PSI_NO_RATE;

        return getPercentile(percentiles, sampleRate);
    }

    public int[] getCohortSpliceSiteData (final String chromosome, int splicePosition)
    {
        final Map<Integer,int[]> totalsMap = mCohortSpliceSiteTotalsMap.get(chromosome);

        if(totalsMap == null)
            return null;

        return totalsMap.get(splicePosition);
    }
}

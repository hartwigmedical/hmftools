package com.hartwig.hmftools.wisp.purity.cn;

import static java.lang.Math.max;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.FileType.AMBER_LOH;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.AMBER_LOH_MIN_AF;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.AMBER_LOH_MIN_TUMOR_BAF;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.addCommonFields;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.addCommonHeaderFields;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.wisp.purity.PurityConfig;
import com.hartwig.hmftools.wisp.purity.ResultsWriter;
import com.hartwig.hmftools.wisp.purity.SampleData;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class AmberLohCalcs
{
    private final PurityConfig mConfig;

    private final SampleData mSample;
    private final List<PurpleCopyNumber> mCopyNumbers;

    private final BufferedWriter mWriter;

    private final Multimap<Chromosome,AmberBAF> mTumorChromosomeBafs;

    public AmberLohCalcs(final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sample)
    {
        mConfig = config;
        mWriter = resultsWriter.getAmberLohWriter();
        mSample = sample;

        mCopyNumbers = Lists.newArrayList();
        mTumorChromosomeBafs = ArrayListMultimap.create();

        try
        {
            String cnFile = PurpleCopyNumberFile.generateFilenameForReading(mConfig.getPurpleDir(sample.TumorId), mSample.TumorId);
            mCopyNumbers.addAll(PurpleCopyNumberFile.read(cnFile));

            String tumorAmberFile = AmberBAFFile.generateAmberFilenameForReading(mConfig.getAmberDir(sample.TumorId), mSample.TumorId);
            mTumorChromosomeBafs.putAll(AmberBAFFile.read(tumorAmberFile, true));
        }
        catch(Exception e)
        {
            CT_LOGGER.error("sample({}) failed to load Purple and Amber data: {}", mSample.TumorId, e.toString());
        }
    }

    public AmberLohResult processSample(final String sampleId, final PurityContext purityContext)
    {
        try
        {
            String sampleAmberFile = AmberBAFFile.generateAmberFilenameForReading(mConfig.getAmberDir(sampleId), sampleId);
            Multimap<Chromosome,AmberBAF> sampleChromosomeBafs = AmberBAFFile.read(sampleAmberFile, true);

            // select only regions with an LOH
            Map<String,List<PurpleCopyNumber>> chrCopyNumbers = Maps.newHashMap();

            for(PurpleCopyNumber copyNumber : mCopyNumbers)
            {
                if(copyNumber.minorAlleleCopyNumber() > 0.5)
                    continue;

                List<PurpleCopyNumber> copyNumbers = chrCopyNumbers.get(copyNumber.chromosome());

                if(copyNumbers == null)
                {
                    copyNumbers = Lists.newArrayList();
                    chrCopyNumbers.put(copyNumber.chromosome(), copyNumbers);
                }

                copyNumbers.add(copyNumber);
            }

            List<SiteData> siteDataList = Lists.newArrayList();

            int totalAmberSites = 0;

            for(Map.Entry<String,List<PurpleCopyNumber>> entry : chrCopyNumbers.entrySet())
            {
                String chrStr = entry.getKey();
                Chromosome chromosome = HumanChromosome.fromString(chrStr);

                Collection<AmberBAF> tumorChrSites = mTumorChromosomeBafs.get(chromosome);
                Collection<AmberBAF> sampleChrites = sampleChromosomeBafs.get(chromosome);

                totalAmberSites += tumorChrSites.size();
                totalAmberSites += sampleChrites.size();

                for(PurpleCopyNumber copyNumber : entry.getValue())
                {
                    List<AmberBAF> tumorLohSites = tumorChrSites.stream()
                            .filter(x -> positionWithin(x.position(), copyNumber.start(), copyNumber.end())).collect(Collectors.toList());

                    List<AmberBAF> sampleLohSites = sampleChrites.stream()
                            .filter(x -> positionWithin(x.position(), copyNumber.start(), copyNumber.end())).collect(Collectors.toList());

                    SiteData siteData = calculateAmberSupport(copyNumber, tumorLohSites, sampleLohSites);

                    if(siteData == null)
                        continue;

                    siteDataList.add(siteData);

                    writeLohData(mWriter, mConfig, mSample, sampleId, siteData);
                }
            }

            List<SiteData> filteredSiteData = siteDataList.stream().filter(x -> x.averageAF() < AMBER_LOH_MIN_AF).collect(Collectors.toList());

            if(filteredSiteData.isEmpty())
                return AmberLohResult.INVALID_RESULT;

            double lohPercent = filteredSiteData.size() / (double)totalAmberSites;

            Collections.sort(filteredSiteData, Comparator.comparingDouble(x -> x.impliedPurity()));
            int medianIndex = filteredSiteData.size() / 2;
            double lohEstimatedPurity = filteredSiteData.get(medianIndex).impliedPurity();

            double lohTotalCopyNumber = 0;
            double lohTotalAf = 0;
            int lohTotalSupport = 0;

            for(SiteData siteData : filteredSiteData)
            {
                lohTotalCopyNumber += siteData.CopyNumber.averageTumorCopyNumber();
                lohTotalSupport += siteData.SupportCount;
                lohTotalAf += siteData.averageAF();
            }

            double filteredSites = filteredSiteData.size();

            double lohMeanCN = lohTotalCopyNumber / filteredSites;
            double lohMeanAf = lohTotalAf / filteredSites;

            PoissonDistribution poissonDistribution = new PoissonDistribution(0.5 * lohTotalSupport);
            int observed = (int)round(lohMeanAf * lohTotalSupport);
            double lohProbability = poissonDistribution.cumulativeProbability(observed);

            // LOHpValue = ppois(medianAF*fragments,0.5*fragments,TRUE)

            return new AmberLohResult(
                    filteredSiteData.size(), lohEstimatedPurity, lohPercent, lohMeanCN, lohMeanAf, lohProbability, lohTotalSupport);
        }
        catch(Exception e)
        {
            CT_LOGGER.error("sample({}) failed to load Amber data: {}", sampleId, e.toString());
            e.printStackTrace();
            return null;
        }
    }

    private class SiteData
    {
        public final PurpleCopyNumber CopyNumber;
        public int SiteCount;
        public int SupportCount;
        public int SampleDepthTotal;
        public double AfTotal;

        public SiteData(final PurpleCopyNumber copyNumber)
        {
            CopyNumber = copyNumber;
            SiteCount = 0;
            SupportCount = 0;
            SampleDepthTotal = 0;
            AfTotal = 0;
        }

        public void merge(final SiteData other)
        {
            SiteCount += other.SiteCount;
            SupportCount += other.SupportCount;
            SampleDepthTotal += other.SampleDepthTotal;
        }

        // public double meanAF() { return SampleDepthTotal > 0 ? SupportCount / (double)SampleDepthTotal : 0; }
        public double averageAF() { return SiteCount > 0 ? AfTotal / (double)SiteCount : 0; }

        public double impliedPurity()
        {
            // impliedPurity=max(0,(1-2*avgAF) / (avgAF*(copyNumber-2) + 1)

            double avgAf = averageAF();
            double denom = avgAf * (CopyNumber.averageTumorCopyNumber() - 2) + 1;
            return denom > 0 ? max((1 - 2 * avgAf) / denom, 0) : 0;
        }
    }

    private SiteData calculateAmberSupport(
            final PurpleCopyNumber copyNumber, final List<AmberBAF> tumorLohSites, final List<AmberBAF> sampleLohSites)
    {
        // find matching sites and accumulate VAF counts
        if(tumorLohSites.isEmpty() || sampleLohSites.isEmpty())
            return null;

        int sampleSiteIndex = 0;
        AmberBAF sampleSite = sampleLohSites.get(sampleSiteIndex);

        SiteData siteData = new SiteData(copyNumber);

        for(AmberBAF tumorSite : tumorLohSites)
        {
            if(tumorSite.tumorBAF() < AMBER_LOH_MIN_TUMOR_BAF)
                continue;

            while(sampleSite.position() < tumorSite.position())
            {
                ++sampleSiteIndex;

                if(sampleSiteIndex >= sampleLohSites.size())
                {
                    sampleSite = null;
                    break;
                }

                sampleSite = sampleLohSites.get(sampleSiteIndex);
            }

            if(sampleSite == null)
                break;

            if(sampleSite.position() > tumorSite.position())
                continue;

            if(sampleSite.position() == tumorSite.position())
            {
                ++siteData.SiteCount;

                double sampleBaf = tumorSite.tumorBAF() < 0.5 ? sampleSite.tumorBAF() : 1 - sampleSite.tumorBAF();
                int sampleSupport = (int)round(sampleSite.tumorDepth() * sampleBaf);

                siteData.SupportCount += sampleSupport;
                siteData.SampleDepthTotal += sampleSite.tumorDepth();

                double af = sampleSite.tumorDepth() > 0 ? sampleSupport / (double)sampleSite.tumorDepth() : 0;
                siteData.AfTotal += af;
            }
        }

        return siteData;
    }

    public static BufferedWriter initialiseAmberLohWriter(final PurityConfig config)
    {
        try
        {
            String fileName = config.formFilename(AMBER_LOH);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            addCommonHeaderFields(sj, config);

            sj.add("Chromosome").add("CnSegmentStart").add("CnSegmentEnd").add("CopyNumber");
            sj.add("SiteCount").add("SupportCount").add("AvgAF").add("ImpliedPurity");

            writer.write(sj.toString());
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to initialise Amber LOH file: {}", e.toString());
            return null;
        }
    }

    private static synchronized void writeLohData(
            final BufferedWriter writer, final PurityConfig config,
            final SampleData sampleData, final String sampleId, final SiteData siteData)
    {
        if(writer == null)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);
            addCommonFields(sj, config, sampleData, sampleId);

            sj.add(siteData.CopyNumber.chromosome());
            sj.add(String.valueOf(siteData.CopyNumber.start()));
            sj.add(String.valueOf(siteData.CopyNumber.end()));
            sj.add(format("%.2f", siteData.CopyNumber.averageTumorCopyNumber()));
            sj.add(String.valueOf(siteData.SiteCount));
            sj.add(String.valueOf(siteData.SupportCount));
            sj.add(format("%.2f", siteData.averageAF()));
            sj.add(format("%.2f", siteData.impliedPurity()));

            writer.write(sj.toString());
            writer.newLine();
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to write Amber LOH file: {}", e.toString());
        }
    }

}
package com.hartwig.hmftools.wisp.purity.cn;

import static java.lang.Math.round;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;

import java.io.File;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
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

public class AmberLohCalcs
{
    private final PurityConfig mConfig;

    private final SampleData mSample;
    private final List<PurpleCopyNumber> mCopyNumbers;

    private final Multimap<Chromosome,AmberBAF> mTumorChromosomeBafs;

    public AmberLohCalcs(final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sample)
    {
        mConfig = config;

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

            double cnTotal = 0;
            int lohDistanceTotal = 0;

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

            SiteData combinedSiteData = new SiteData();

            for(Map.Entry<String,List<PurpleCopyNumber>> entry : chrCopyNumbers.entrySet())
            {
                String chrStr = entry.getKey();
                Chromosome chromosome = HumanChromosome.fromString(chrStr);

                Collection<AmberBAF> tumorChrSites = mTumorChromosomeBafs.get(chromosome);
                Collection<AmberBAF> sampleChrites = sampleChromosomeBafs.get(chromosome);

                for(PurpleCopyNumber copyNumber : entry.getValue())
                {
                    List<AmberBAF> tumorLohSites = tumorChrSites.stream()
                            .filter(x -> positionWithin(x.position(), copyNumber.start(), copyNumber.end())).collect(Collectors.toList());

                    List<AmberBAF> sampleLohSites = sampleChrites.stream()
                            .filter(x -> positionWithin(x.position(), copyNumber.start(), copyNumber.end())).collect(Collectors.toList());

                    SiteData siteData = calculateAmberSupport(copyNumber, tumorLohSites, sampleLohSites);

                    if(siteData == null)
                        continue;

                    combinedSiteData.merge(siteData);

                    lohDistanceTotal += copyNumber.length();
                    cnTotal += copyNumber.length() * copyNumber.averageTumorCopyNumber();
                }
            }

            double avgCopyNumber = cnTotal / lohDistanceTotal;
            double avgAf = combinedSiteData.SupportCount / (double)combinedSiteData.SampleDepthTotal;

            double estimatedPurity = (1 - 2 * avgAf) / (avgAf * (avgCopyNumber - 2) + 1);

            return new AmberLohResult(combinedSiteData.SiteCount, lohDistanceTotal, avgCopyNumber, avgAf, estimatedPurity);

            /*
            AF = sum(temp$support) / sum(temp$tumorDepthP)
            avgCN = mean(temp$copyNumber)
            purity = (1-2*AF) / (AF*(avgCN-2) + 1)
            */
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
        public int SiteCount;
        public int SupportCount;
        public int SampleDepthTotal;

        public SiteData()
        {
            SiteCount = 0;
            SupportCount = 0;
            SampleDepthTotal = 0;
        }

        public void merge(final SiteData other)
        {
            SiteCount += other.SiteCount;
            SupportCount += other.SupportCount;
            SampleDepthTotal += other.SampleDepthTotal;
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

        SiteData siteData = new SiteData();

        for(AmberBAF tumorSite : tumorLohSites)
        {
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

                // temp = temp %>% mutate(support=round(ifelse(tumorBAFT<0.5,tumorBAFP,1-tumorBAFP)*tumorDepthP,0))
                double sampleBaf = tumorSite.tumorBAF() < 0.5 ? sampleSite.tumorBAF() : 1 - sampleSite.tumorBAF();
                int sampleSupport = (int)round(sampleSite.tumorDepth() * sampleBaf);

                siteData.SupportCount += sampleSupport;

                siteData.SampleDepthTotal += sampleSite.tumorDepth();
            }
        }

        return siteData;
    }
}
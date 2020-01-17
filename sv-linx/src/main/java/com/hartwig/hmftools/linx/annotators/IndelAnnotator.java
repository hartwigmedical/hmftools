package com.hartwig.hmftools.linx.annotators;

import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.types.ResolvedType.COMPLEX;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_INV;
import static com.hartwig.hmftools.linx.types.ResolvedType.SIMPLE_GRP;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.linx.analysis.ClusterMetrics;
import com.hartwig.hmftools.linx.types.ResolvedType;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class IndelAnnotator
{
    private final DatabaseAccess mDbAccess;

    private final Map<String,List<IndelData>> mSampleIndelData;

    private final Map<String,List<IndelData>> mSampleChrIndelData;
    private int mIndelCount;
    private int mMsiIndelCount;

    private static final int INDEL_THRESHOLD = 500;
    private static final int MSI_THRESHOLD = 250;

    private static final Logger LOGGER = LogManager.getLogger(IndelAnnotator.class);

    public IndelAnnotator(final DatabaseAccess dbAccess, final String indelFile)
    {
        mDbAccess = dbAccess;
        mSampleIndelData = Maps.newHashMap();
        mSampleChrIndelData = Maps.newHashMap();

        if(indelFile != null && !indelFile.isEmpty())
        {
            loadIndelsFromFile(indelFile);
        }
    }

    public final Map<String,List<IndelData>> getSampleIndelData() { return mSampleIndelData; }

    public void loadIndels(final String sampleId)
    {
        mSampleChrIndelData.clear();
        mIndelCount = 0;
        mMsiIndelCount = 0;

        final List<IndelData> sampleIndelData = Lists.newArrayList();

        if(mSampleIndelData.containsKey(sampleId))
        {
            sampleIndelData.addAll(mSampleIndelData.get(sampleId));
        }
        else if(mDbAccess != null)
        {
            sampleIndelData.addAll(loadIndelsFromDatabase(sampleId));
        }

        if(sampleIndelData.isEmpty())
            return;

        mIndelCount = sampleIndelData.size();

        for(IndelData indel : sampleIndelData)
        {
            // only cache non-MSI indels and record the MSI count
            if (indel.highRepeatCount())
            {
                ++mMsiIndelCount;
                continue;
            }

            List<IndelData> indelDataList = mSampleChrIndelData.get(indel.Chromosome);

            if (indelDataList == null)
            {
                indelDataList = Lists.newArrayList();
                mSampleChrIndelData.put(indel.Chromosome, indelDataList);
            }

            indelDataList.add(indel);
        }
    }

    private List<IndelData> findIndels(final String chromosome, long posStart, long posEnd)
    {
        final List<IndelData> indelDataList = mSampleChrIndelData.get(chromosome);

        if(indelDataList == null)
            return Lists.newArrayList();

        return indelDataList.stream().filter(x -> x.Position >= posStart && x.Position <= posEnd).collect(Collectors.toList());
    }

    public boolean exceedsThresholds()
    {
        // ignore samples with high INDEL count or MSI evidence
        if(mIndelCount >= INDEL_THRESHOLD || mMsiIndelCount >= MSI_THRESHOLD)
        {
            LOGGER.debug("skipped indel annotation, high indel({}) msi({}) counts", mIndelCount, mMsiIndelCount);
            return true;
        }

        return false;
    }

    private static final double GENOME_LENGTH = 3e9;

    public void annotateCluster(final SvCluster cluster)
    {
        if(cluster.getResolvedType() != COMPLEX && cluster.getResolvedType() != SIMPLE_GRP && cluster.getResolvedType() != RECIP_INV)
            return;

        int totalIndels = 0;

        for(SvLinkedPair pair : cluster.getLinkedPairs())
        {
            long lowerPos = pair.getBreakend(true).position();
            long upperPos = pair.getBreakend(false).position();
            pair.setIndelCount(findIndels(pair.chromosome(), lowerPos, upperPos).size());

            if(pair.getIndelCount() > 0)
            {
                totalIndels += pair.getIndelCount();
                LOGGER.debug("cluster({}) pair({} len={}) indelCount({})",
                        cluster.id(), pair, pair.length(), pair.getIndelCount());
            }
        }

        if(totalIndels > 1 && !cluster.hasVariedPloidy())
        {
            ClusterMetrics metrics = cluster.getMetrics();

            metrics.IndelCount = totalIndels;
            metrics.IndelProbability = 1;

            long clusterRange = metrics.TotalRange;

            if (clusterRange > 0)
            {
                double expectedIndelCount = min(clusterRange / GENOME_LENGTH, 1) * (mIndelCount - mMsiIndelCount);

                if (totalIndels > expectedIndelCount)
                {
                    PoissonDistribution poisson = new PoissonDistribution(expectedIndelCount);
                    double indelProbability = 1 - poisson.cumulativeProbability(totalIndels - 1);
                    metrics.IndelProbability = indelProbability;

                    if(indelProbability < 0.01)
                    {
                        LOGGER.debug(String.format("cluster(%d) indelCount(%d) range(%d) expected(%.1f) probability(%.9f)",
                                cluster.id(), totalIndels, clusterRange, expectedIndelCount, indelProbability));
                    }
                }
            }
        }
    }

    private static final int CSV_REQUIRED_FIELDS = 8;
    private final static int INDEL_COL_SAMPLE = 0;
    private final static int INDEL_COL_CHR = 1;
    private final static int INDEL_COL_POS = 2;
    private final static int INDEL_COL_REF = 3;
    private final static int INDEL_COL_ALT = 4;
    private final static int INDEL_COL_MH = 5;
    private final static int INDEL_COL_RC = 6;
    private final static int INDEL_COL_PLOIDY = 7;

    private void loadIndelsFromFile(final String filename)
    {
        if(filename.isEmpty())
            return;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            int recordCount = 0;

            // SampleId,Chromosome,Position,Ref,Alt,Microhomology,RepeatCount,Ploidy

            String line;
            String currentSample = "";
            List<IndelData> indelDataList = null;

            while ((line = fileReader.readLine()) != null)
            {
                if(line.contains("SampleId"))
                    continue;

                // parse CSV data
                String[] items = line.split(",");

                if(items.length < CSV_REQUIRED_FIELDS)
                    continue;

                final String sampleId = items[INDEL_COL_SAMPLE];

                IndelData indel = new IndelData(items[INDEL_COL_CHR], Long.parseLong(items[INDEL_COL_POS]),
                        items[INDEL_COL_REF], items[INDEL_COL_ALT], items[INDEL_COL_MH],
                        Integer.parseInt(items[INDEL_COL_RC]), Double.parseDouble(items[INDEL_COL_PLOIDY]));

                ++recordCount;

                if(!currentSample.equals(sampleId))
                {
                    currentSample = sampleId;
                    indelDataList = Lists.newArrayList();
                    mSampleIndelData.put(sampleId, indelDataList);
                }

                indelDataList.add(indel);
            }

            LOGGER.debug("loaded {} indel data records, samples({})", recordCount, mSampleIndelData.size());
        }
        catch(IOException exception)
        {
            LOGGER.error("Failed to read indel CSV file({})", filename);
        }
    }

    private List<IndelData> loadIndelsFromDatabase(final String sampleId)
    {
        final List<SomaticVariant> variants = mDbAccess.readSomaticVariants(sampleId, VariantType.INDEL);
        List<IndelData> sampleIndelList = Lists.newArrayList();

        for(final SomaticVariant variant : variants)
        {
            if (variant.isFiltered())
                continue;

            sampleIndelList.add(new IndelData(variant.chromosome(), variant.position(), variant.ref(), variant.alt(),
                    variant.microhomology(), variant.repeatCount(), variant.ploidy()));
        }

        LOGGER.debug("sample({}) retrieved {} indels", sampleId, sampleIndelList.size());

        return sampleIndelList;
    }

}

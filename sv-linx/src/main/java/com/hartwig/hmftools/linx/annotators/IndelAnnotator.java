package com.hartwig.hmftools.linx.annotators;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.annotators.IndelData.CSV_REQUIRED_FIELDS;
import static com.hartwig.hmftools.linx.annotators.IndelData.INDEL_COL_SAMPLE;
import static com.hartwig.hmftools.linx.annotators.IndelData.fromString;
import static com.hartwig.hmftools.linx.types.ResolvedType.COMPLEX;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_INV;
import static com.hartwig.hmftools.linx.types.ResolvedType.SIMPLE_GRP;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.analysis.ClusterMetrics;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.math3.distribution.PoissonDistribution;

public class IndelAnnotator
{
    private final DatabaseAccess mDbAccess;

    private final Map<String,List<IndelData>> mSampleIndelData;

    private final Map<String,List<IndelData>> mSampleChrIndelData;
    private int mIndelCount;
    private int mMsiIndelCount;
    private BufferedWriter mFileWriter;

    private static final int INDEL_THRESHOLD = 500;
    private static final int MSI_THRESHOLD = 250;

    public IndelAnnotator(final DatabaseAccess dbAccess, final LinxConfig config)
    {
        mDbAccess = dbAccess;
        mSampleIndelData = Maps.newHashMap();
        mSampleChrIndelData = Maps.newHashMap();

        if(config.IndelFile != null && !config.IndelFile.isEmpty())
        {
            loadIndelsFromFile(config.IndelFile);
        }
        else if(mDbAccess != null && config.Output.WriteAll)
        {
            createOutputFile(config.OutputDataPath);
        }
    }

    public final Map<String,List<IndelData>> getSampleIndelData() { return mSampleIndelData; }

    public void close() { closeBufferedWriter(mFileWriter);}

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
            writeIndexData(sampleId, sampleIndelData);
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

    private List<IndelData> findIndels(final String chromosome, int posStart, int posEnd)
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
            LNX_LOGGER.debug("skipped indel annotation, high indel({}) msi({}) counts", mIndelCount, mMsiIndelCount);
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

        for(LinkedPair pair : cluster.getLinkedPairs())
        {
            int lowerPos = pair.getBreakend(true).position();
            int upperPos = pair.getBreakend(false).position();
            pair.setIndelCount(findIndels(pair.chromosome(), lowerPos, upperPos).size());

            if(pair.getIndelCount() > 0)
            {
                totalIndels += pair.getIndelCount();
                LNX_LOGGER.debug("cluster({}) pair({} len={}) indelCount({})",
                        cluster.id(), pair, pair.length(), pair.getIndelCount());
            }
        }

        ClusterMetrics metrics = cluster.getMetrics();
        metrics.IndelCount = totalIndels;
        metrics.IndelProbability = 1; // default indicates nothing unexpected

        if(!cluster.hasVariedJcn() && metrics.TotalRange > 0 && metrics.TotalDeleted < metrics.TotalRange)
        {
            long clusterRange = metrics.TotalRange - metrics.TotalDeleted;
            double expectedIndelCount = min(clusterRange / GENOME_LENGTH, 1) * (mIndelCount - mMsiIndelCount);

            if (totalIndels > expectedIndelCount)
            {
                PoissonDistribution poisson = new PoissonDistribution(expectedIndelCount);
                double indelProbability = 1 - poisson.cumulativeProbability(totalIndels - 1);
                metrics.IndelProbability = indelProbability;

                if(indelProbability < 0.01)
                {
                    LNX_LOGGER.debug(String.format("cluster(%d) indelCount(%d) range(%d) expected(%.1f) probability(%.9f)",
                            cluster.id(), totalIndels, clusterRange, expectedIndelCount, indelProbability));
                }
            }
        }
    }


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

                IndelData indel = fromString(items);

                ++recordCount;

                if(!currentSample.equals(sampleId))
                {
                    currentSample = sampleId;
                    indelDataList = Lists.newArrayList();
                    mSampleIndelData.put(sampleId, indelDataList);
                }

                indelDataList.add(indel);
            }

            LNX_LOGGER.debug("loaded {} indel data records, samples({})", recordCount, mSampleIndelData.size());
        }
        catch(IOException exception)
        {
            LNX_LOGGER.error("Failed to read indel CSV file({})", filename);
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

            sampleIndelList.add(new IndelData(variant.chromosome(), (int)variant.position(), variant.ref(), variant.alt(),
                    variant.microhomology(), variant.repeatCount(), variant.variantCopyNumber()));
        }

        LNX_LOGGER.debug("sample({}) retrieved {} indels", sampleId, sampleIndelList.size());

        return sampleIndelList;
    }

    private void createOutputFile(final String outputDir)
    {
        try
        {
            String outputFileName = outputDir + "LNX_INDELS.csv";

            mFileWriter = createBufferedWriter(outputFileName, false);
            mFileWriter.write("SampleId,Chromosome,Position,Ref,Alt,Microhomology,RepeatCount,Ploidy");
            mFileWriter.newLine();
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing indel output file: {}", e.toString());
        }
    }

    private void writeIndexData(final String sampleId, final List<IndelData> indelDataList)
    {
        if(mFileWriter == null)
            return;

        try
        {
            for(final IndelData data : indelDataList)
            {
                // SampleId,Chromosome,Position,Ref,Alt,Microhomology,RepeatCount,Ploidy
                mFileWriter.write(String.format("%s,%s,%d,%s,%s,%s,%d,%.4f",
                        sampleId, data.Chromosome, data.Position, data.Ref, data.Alt,
                        data.Microhomology, data.RepeatCount, data.Ploidy));
                mFileWriter.newLine();
            }
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing indel output file: {}", e.toString());
        }
    }


}

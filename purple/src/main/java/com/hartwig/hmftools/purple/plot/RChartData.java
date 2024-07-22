package com.hartwig.hmftools.purple.plot;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.KATAEGIS_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_VARIANT_CN;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.concurrent.atomic.AtomicInteger;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.purple.PurpleConfig;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.VariantContext;

public class RChartData
{
    private static final double COPY_NUMBER_BUCKET_SIZE = 1;
    private static final double VARIANT_COPY_NUMBER_BUCKET_SIZE = 0.05;
    private static final int MAX_SOMATIC_PLOT_COUNT = 100_000;

    private final Map<String, AtomicInteger> mSomaticHistogram = Maps.newHashMap();
    private final String mHistogramFilename;
    private final BufferedWriter mSomaticWriter;
    private int mSomaticCount;

    public RChartData(final PurpleConfig config, final String tumorSample)
    {
        mHistogramFilename = config.OutputDir + tumorSample + ".purple.somatic.hist.tsv";

        String somaticFilename = somaticDataFilename(config, tumorSample);
        mSomaticWriter = !config.Charting.Disabled ? initialiseSomaticWriter(somaticFilename) : null;
        mSomaticCount = 0;
    }

    private static String somaticDataFilename(final PurpleConfig config, final String tumorSample)
    {
        return config.Charting.PlotDirectory + tumorSample + ".somatic_data.tsv";
    }

    public void processVariant(final SomaticVariant variant)
    {
        somaticVariantCopyNumberPdf(variant.context());

        if(mSomaticCount >= MAX_SOMATIC_PLOT_COUNT)
            return;

        ++mSomaticCount;

        writeSomaticData(variant);
    }

    public void write() throws IOException
    {
        Files.write(new File(mHistogramFilename).toPath(), variantCopyNumberByCopyNumberString());
        closeBufferedWriter(mSomaticWriter);
    }

    public static void cleanupFiles(final PurpleConfig config, final String tumorSample)
    {
        String somaticFilename = somaticDataFilename(config, tumorSample);

        if(Files.exists(Paths.get(somaticFilename)))
        {
            try { Files.delete(Paths.get(somaticFilename)); } catch(IOException e) {}
        }
    }

    private List<String> variantCopyNumberByCopyNumberString()
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(histogramHeader());
        mSomaticHistogram.entrySet().stream().map(RChartData::histogramEntryToString).sorted().forEach(lines::add);
        return lines;
    }

    private static String histogramHeader()
    {
        return new StringJoiner(TSV_DELIM, "", "")
                .add("variantCopyNumberBucket")
                .add("copyNumberBucket")
                .add("count")
                .toString();
    }

    private BufferedWriter initialiseSomaticWriter(final String somaticFilename)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(somaticFilename);

            // write data required for the rainfall plot
            writer.write("chromosome\tposition\tmutation\tkataegis");
            writer.newLine();
            return writer;

        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to write plotting somatics: {}", e.toString());
            return null;
        }
    }

    private void writeSomaticData(final SomaticVariant variant)
    {
        if(mSomaticWriter == null)
            return;

        if(variant.type() != VariantType.SNP || !variant.isPass())
            return;

        try
        {
            String mutation = format("%s>%s", variant.decorator().ref(), variant.decorator().alt());
            String kataegis = variant.context().getAttributeAsString(KATAEGIS_FLAG, "");

            mSomaticWriter.write(format("%s\t%d\t%s\t%s",
                    variant.chromosome(), variant.position(), mutation, kataegis));
            mSomaticWriter.newLine();
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to write plotting somatics: {}", e.toString());
        }
    }

    private static String histogramEntryToString(final Map.Entry<String,AtomicInteger> entry)
    {
        String[] keys = entry.getKey().split(">");

        return new StringJoiner(TSV_DELIM)
                .add(format("%.2f", Integer.parseInt(keys[0]) * VARIANT_COPY_NUMBER_BUCKET_SIZE))
                .add(format("%.0f", Integer.parseInt(keys[1]) * COPY_NUMBER_BUCKET_SIZE))
                .add(String.valueOf(entry.getValue()))
                .toString();
    }

    private void somaticVariantCopyNumberPdf(final VariantContext somaticVariant)
    {
        CommonInfo commonInfo = somaticVariant.getCommonInfo();
        double copyNumber = commonInfo.getAttributeAsDouble(PURPLE_CN, 0.0);
        double variantCopyNumber = commonInfo.getAttributeAsDouble(PURPLE_VARIANT_CN, 0.0);

        int copyNumberBucket = bucket(copyNumber, COPY_NUMBER_BUCKET_SIZE);
        int variantCopyNumberBucket = bucket(variantCopyNumber, VARIANT_COPY_NUMBER_BUCKET_SIZE);

        final String key = variantCopyNumberBucket + ">" + copyNumberBucket;
        mSomaticHistogram.computeIfAbsent(key, x -> new AtomicInteger()).incrementAndGet();
    }

    static int bucket(double value, double binWidth)
    {
        return (int) Math.round((value) / binWidth);
    }
}

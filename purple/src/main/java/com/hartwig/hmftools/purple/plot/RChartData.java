package com.hartwig.hmftools.purple.plot;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_VARIANT_CN_INFO;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.concurrent.atomic.AtomicInteger;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.VariantContext;

public class RChartData
{
    private static final double COPY_NUMBER_BUCKET_SIZE = 1;
    private static final double VARIANT_COPY_NUMBER_BUCKET_SIZE = 0.05;
    private static final DecimalFormat VCN_FORMAT = new DecimalFormat("0.00");
    private static final DecimalFormat CN_FORMAT = new DecimalFormat("0");

    private static final String DELIMITER = "\t";

    private final Map<String, AtomicInteger> mSomaticHistogram = Maps.newHashMap();
    private final String mFilename;

    public RChartData(String outputDirectory, String tumorSample)
    {
        mFilename = outputDirectory + File.separator + tumorSample + ".purple.somatic.hist.tsv";
    }

    public void processVariant(final SomaticVariant variant)
    {
        somaticVariantCopyNumberPdf(variant.context());
    }

    public void write() throws IOException
    {
        Files.write(new File(mFilename).toPath(), variantCopyNumberByCopyNumberString());
    }

    private List<String> variantCopyNumberByCopyNumberString()
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        mSomaticHistogram.entrySet().stream().map(RChartData::toString).sorted().forEach(lines::add);
        return lines;
    }

    private static String header()
    {
        return new StringJoiner(DELIMITER, "", "").add("variantCopyNumberBucket").add("copyNumberBucket").add("count").toString();
    }

    private static String toString(Map.Entry<String, AtomicInteger> entry)
    {
        String[] keys = entry.getKey().split(">");

        return new StringJoiner(DELIMITER).add(VCN_FORMAT.format(Integer.parseInt(keys[0]) * VARIANT_COPY_NUMBER_BUCKET_SIZE))
                .add(CN_FORMAT.format(Integer.parseInt(keys[1]) * COPY_NUMBER_BUCKET_SIZE))
                .add(String.valueOf(entry.getValue()))
                .toString();
    }

    private void somaticVariantCopyNumberPdf(final VariantContext somaticVariant)
    {
        CommonInfo commonInfo = somaticVariant.getCommonInfo();
        double copyNumber = commonInfo.getAttributeAsDouble(PURPLE_CN_INFO, 0.0);
        double variantCopyNumber = commonInfo.getAttributeAsDouble(PURPLE_VARIANT_CN_INFO, 0.0);

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

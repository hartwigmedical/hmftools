package com.hartwig.hmftools.wisp.purity.variant;

import static com.hartwig.hmftools.common.sage.FragmentLengthCounts.ALT_COUNT;
import static com.hartwig.hmftools.common.sage.FragmentLengthCounts.REF_COUNT;
import static com.hartwig.hmftools.common.sage.VariantFragmentLength.FLD_ALT_COUNT;
import static com.hartwig.hmftools.common.sage.VariantFragmentLength.FLD_LENGTH;
import static com.hartwig.hmftools.common.sage.VariantFragmentLength.FLD_REF_COUNT;
import static com.hartwig.hmftools.common.sage.VariantFragmentLength.VARIANT_FRAG_LENGTHS_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.FileType.FRAGMENT_LENGTHS;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.addCommonFields;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.addCommonHeaderFields;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sage.FragmentLengthCounts;
import com.hartwig.hmftools.common.sage.VariantFragmentLength;
import com.hartwig.hmftools.wisp.purity.PurityConfig;
import com.hartwig.hmftools.wisp.purity.ResultsWriter;
import com.hartwig.hmftools.wisp.purity.SampleData;

public class SampleFragmentLengths
{
    private final SampleData mSample;
    private final PurityConfig mConfig;
    private final ResultsWriter mResultsWriter;

    private final Map<String,Boolean> mVariantStatusMap;

    public SampleFragmentLengths(final PurityConfig config, final ResultsWriter resultsWriter, final SampleData sample)
    {
        mConfig = config;
        mResultsWriter = resultsWriter;
        mSample = sample;

        mVariantStatusMap = Maps.newHashMap();
    }

    public void processSample(final String variantVcf, final List<SomaticVariant> variants)
    {
        String fragmentLengthFile = variantVcf.replace(".vcf.gz", VARIANT_FRAG_LENGTHS_FILE_ID);

        List<VariantFragmentLength> fragmentLengths = VariantFragmentLength.read(fragmentLengthFile);

        if(fragmentLengths == null || fragmentLengths.isEmpty())
            return;

        Map<String,FragmentLengthCounts> sampleLengthDistributions = Maps.newHashMap();

        // form a combined length distribution for non-filtered variants
        String currentVariantInfo = "";
        boolean currentUseVariant = false;

        for(VariantFragmentLength varFragLength : fragmentLengths)
        {
            boolean useVariant;

            if(!varFragLength.VariantInfo.equals(currentVariantInfo))
            {
                currentVariantInfo = varFragLength.VariantInfo;
                currentUseVariant = useVariant(varFragLength.VariantInfo, variants);
            }

            useVariant = currentUseVariant;

            if(!useVariant)
                continue;

            FragmentLengthCounts distribution = sampleLengthDistributions.get(varFragLength.SampleId);

            if(distribution == null)
            {
                distribution = new FragmentLengthCounts();
                sampleLengthDistributions.put(varFragLength.SampleId, distribution);
            }

            int[] counts = distribution.lengthCounts().get(varFragLength.Length);

            if(counts == null)
            {
                counts = new int[2];
                distribution.lengthCounts().put(varFragLength.Length, counts);
            }

            counts[REF_COUNT] += varFragLength.RefCount;
            counts[ALT_COUNT] += varFragLength.AltCount;
        }

        for(Map.Entry<String,FragmentLengthCounts> entry : sampleLengthDistributions.entrySet())
        {
            String sampleId = entry.getKey();
            FragmentLengthCounts counts = entry.getValue();

            writeSampleFragLengths(mResultsWriter.getFragLengthWriter(), mConfig, mSample, sampleId, counts);
        }
    }

    private boolean useVariant(final String variantInfo, final List<SomaticVariant> variants)
    {
        Boolean useVariant = mVariantStatusMap.get(variantInfo);

        if(useVariant != null)
            return useVariant;

        // chr1:34610987 G>T
        String[] parts = variantInfo.split(" ");
        String[] coords = parts[0].split(":");
        String[] refAlt = parts[1].split(">");
        String chromosome = coords[0];
        int position = Integer.parseInt(coords[1]);
        String ref = refAlt[0];
        String alt = refAlt[1];

        for(SomaticVariant variant : variants)
        {
            if(variant.Chromosome.equals(chromosome) && variant.Position == position && variant.Ref.equals(ref) && variant.Alt.equals(alt))
            {
                useVariant = !variant.isFiltered();
                mVariantStatusMap.put(variantInfo, useVariant);

                return useVariant;
            }
        }

        return false;
    }

    public static BufferedWriter initialiseFragmentLengthWriter(final PurityConfig config)
    {
        try
        {
            String fileName = config.formFilename(FRAGMENT_LENGTHS);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            addCommonHeaderFields(sj, config);

            sj.add(FLD_LENGTH).add(FLD_REF_COUNT).add(FLD_ALT_COUNT);

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to initialise fragment length output file: {}", e.toString());
            return null;
        }
    }

    private static synchronized void writeSampleFragLengths(
            final BufferedWriter writer, final PurityConfig config,
            final SampleData sampleData, final String sampleId, final FragmentLengthCounts fragLengthCounts)
    {
        try
        {

            for(Map.Entry<Integer,int[]> entry : fragLengthCounts.lengthCounts().entrySet())
            {
                StringJoiner sj = new StringJoiner(TSV_DELIM);

                addCommonFields(sj, config, sampleData, sampleId);

                sj.add(String.valueOf(entry.getKey()));

                final int[] counts = entry.getValue();
                sj.add(String.valueOf(counts[REF_COUNT]));
                sj.add(String.valueOf(counts[ALT_COUNT]));

                writer.write(sj.toString());
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to write fragment length file: {}", e.toString());
        }
    }
}

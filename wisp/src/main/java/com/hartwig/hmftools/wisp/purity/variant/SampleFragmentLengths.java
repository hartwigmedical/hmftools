package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sage.FragmentLengthCounts.ALT_COUNT;
import static com.hartwig.hmftools.common.sage.FragmentLengthCounts.REF_COUNT;
import static com.hartwig.hmftools.common.sage.VariantFragmentLength.FLD_ALT_COUNT;
import static com.hartwig.hmftools.common.sage.VariantFragmentLength.FLD_LENGTH;
import static com.hartwig.hmftools.common.sage.VariantFragmentLength.FLD_REF_COUNT;
import static com.hartwig.hmftools.common.sage.VariantFragmentLength.VARIANT_FRAG_LENGTHS_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.pathFromFile;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.FileType.FRAGMENT_LENGTHS;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.addCommonFields;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.addCommonHeaderFields;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.BaseRegion;
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

    private List<VariantFragmentLength> loadFragmentLengths(final String variantVcf)
    {
        List<VariantFragmentLength> fragmentLengths = Lists.newArrayList();

        // either load from a file named similarly to the VCF or a patient file incremented with multiples of samples
        String somaticVcf = filenamePart(variantVcf);
        String somaticVcfDir = pathFromFile(variantVcf);
        String fragLengthsDir = mConfig.FragmentLengthDir != null ? mConfig.FragmentLengthDir : somaticVcfDir;
        String fragmentLengthFile = fragLengthsDir + somaticVcf.replace(".vcf.gz", VARIANT_FRAG_LENGTHS_FILE_ID);

        boolean foundFiles = false;

        if(Files.exists(Paths.get(fragmentLengthFile)))
        {
            foundFiles = true;
            fragmentLengths.addAll(VariantFragmentLength.read(fragmentLengthFile));
        }
        else
        {
            int index = 0;
            // P014501, like: P014501.sage.frag_lengths.0.tsv.gz
            fragmentLengthFile = fragLengthsDir + format("%s.sage.frag_lengths.%d.tsv.gz", mSample.PatientId, index);

            while(Files.exists(Paths.get(fragmentLengthFile)))
            {
                foundFiles = true;
                fragmentLengths.addAll(VariantFragmentLength.read(fragmentLengthFile));
                index += 3;
                fragmentLengthFile = fragLengthsDir + format("%s.sage.frag_lengths.%d.tsv.gz", mSample.PatientId, index);
            }

            if(!fragmentLengths.isEmpty())
            {
                CT_LOGGER.debug("patient({}) loaded {} fragment entries from {} files",
                        mSample.PatientId, fragmentLengths.size(), index / 3);
            }
        }

        if(!foundFiles)
        {
            CT_LOGGER.warn("patient({}) found no fragment files",mSample.PatientId);
        }

        return fragmentLengths;
    }

    public void processSample(final String variantVcf, final List<SomaticVariant> variants)
    {
        List<VariantFragmentLength> fragmentLengths = loadFragmentLengths(variantVcf);

        if(fragmentLengths == null || fragmentLengths.isEmpty())
            return;

        Map<String,List<SomaticVariant>> variantsMap = Maps.newHashMap();

        String currentChr = "";
        List<SomaticVariant> currentChrVariants = null;

        for(SomaticVariant variant : variants)
        {
            if(!variant.Chromosome.equals(currentChr))
            {
                currentChr = variant.Chromosome;
                currentChrVariants = Lists.newArrayList();
                variantsMap.put(variant.Chromosome, currentChrVariants);
            }

            currentChrVariants.add(variant);
        }

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
                currentUseVariant = useVariant(varFragLength.VariantInfo, variantsMap);
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

    private boolean useVariant(final String variantInfo, final Map<String,List<SomaticVariant>> variantsMap)
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

        List<SomaticVariant> variants = variantsMap.get(chromosome);

        if(variants == null)
            return false;

        SomaticVariant variant = findMatchingVariant(position, ref, alt, variants);

        if(variant == null)
            return false;

        useVariant = !variant.isFiltered();
        mVariantStatusMap.put(variantInfo, useVariant);
        return useVariant;
    }

    private static final int NO_POSITION_MATCH = -1;

    private SomaticVariant findMatchingVariant(final int position, final String ref, final String alt, final List<SomaticVariant> variants)
    {
        int matchedIndex = NO_POSITION_MATCH;

        if(variants.size() < 10)
        {
            for(int i = 0; i < variants.size(); ++i)
            {
                if(variants.get(i).Position == position)
                {
                    matchedIndex = i;
                    break;
                }
            }
        }
        else
        {
            int lowIndex = 0;
            int highIndex = variants.size() - 1;
            int currentIndex = variants.size() / 2;

            while(true)
            {
                int currentPosition = variants.get(currentIndex).Position;

                if(currentPosition == position)
                {
                    matchedIndex = currentIndex;
                    break;
                }

                if(position < currentPosition)
                {
                    // current index is looking too high in the list
                    if(currentIndex == lowIndex + 1)
                        break;

                    highIndex = currentIndex;
                }
                else
                {
                    if(currentIndex == highIndex - 1)
                        break;

                    lowIndex = currentIndex;
                }

                int newIndex = lowIndex + (highIndex - lowIndex) / 2;

                if(newIndex == currentIndex)
                    break;

                currentIndex = newIndex;
            }
        }

        if(matchedIndex == NO_POSITION_MATCH)
            return null;

        SomaticVariant variant = variants.get(matchedIndex);

        if(variant.Alt.equals(alt) && variant.Ref.equals(ref))
            return variant;

        // check for ref & alt match
        for(int i = 0; i <= 1; ++i)
        {
            boolean searchBack = (i == 0);

            int index = searchBack ? matchedIndex - 1 : matchedIndex + 1;

            while(true)
            {
                if(index < 0 || index >= variants.size())
                    break;

                variant = variants.get(index);

                if(variant.Position != position)
                    break;

                if(variant.Alt.equals(alt) && variant.Ref.equals(ref))
                    return variant;

                if(searchBack)
                    --index;
                else
                    ++index;
            }
        }

        return null;
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

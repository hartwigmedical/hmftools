package com.hartwig.hmftools.sig_analyser.loaders;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.MIN_SAMPLE_PURITY;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.OUTPUT_DIR;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.OUTPUT_FILE_ID;
import static com.hartwig.hmftools.sig_analyser.SigAnalyser.SAMPLE_IDS;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.getNewFile;
import static com.hartwig.hmftools.sig_analyser.common.DataUtils.writeMatrixData;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.sig_analyser.common.SigMatrix;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SigSnvLoader
{
    private String mOutputDir;
    private String mOutputFileId;
    private DatabaseAccess mDbAccess;
    private boolean mApplySampleQC;

    final Map<String,Integer> mBucketStringToIndex;
    SigMatrix mSampleBucketCounts;
    final List<String> mSampleIds;

    private Double mPloidyMin;
    private Double mPloidyMax;
    private Double mSubclonalLikelihoodMin;
    private Double mSubclonalLikelihoodMax;

    private static final int SNV_BUCKET_COUNT = 96;
    private static final String APPLY_SAMPLE_QC = "apply_sample_qc";
    private static final String SUBCLONAL_MIN = "subclonal_min";
    private static final String SUBCLONAL_MAX = "subclonal_max";
    private static final String PLOIDY_MAX = "ploidy_max";
    private static final String PLOIDY_MIN = "ploidy_min";

    private static final Logger LOGGER = LogManager.getLogger(SigSnvLoader.class);

    public SigSnvLoader()
    {
        mOutputDir = "";
        mOutputFileId = "";
        mBucketStringToIndex = Maps.newHashMap();
        mSampleIds = Lists.newArrayList();
        mApplySampleQC = true;
        mPloidyMin = null;
        mPloidyMax = null;
        mSubclonalLikelihoodMin = null;
        mSubclonalLikelihoodMax = null;
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(APPLY_SAMPLE_QC, false, "Check sample QC status etc");
        options.addOption(SUBCLONAL_MAX, true, "Optional: subclonal max threshold");
        options.addOption(SUBCLONAL_MIN, true, "Optional: subclonal min threshold");
        options.addOption(PLOIDY_MAX, true, "Optional: ploidy max threshold");
        options.addOption(PLOIDY_MIN, true, "Optional: ploidy min threshold");
    }

    public boolean initialise(final DatabaseAccess dbAccess, final CommandLine cmd)
    {
        mOutputDir = cmd.getOptionValue(OUTPUT_DIR);
        mOutputFileId = cmd.getOptionValue(OUTPUT_FILE_ID);
        mDbAccess = dbAccess;

        mApplySampleQC = cmd.hasOption(APPLY_SAMPLE_QC);

        if(cmd.hasOption(SUBCLONAL_MIN))
            mSubclonalLikelihoodMin = Double.parseDouble(cmd.getOptionValue(SUBCLONAL_MIN, "0"));

        if(cmd.hasOption(SUBCLONAL_MAX))
            mSubclonalLikelihoodMax = Double.parseDouble(cmd.getOptionValue(SUBCLONAL_MAX, "0"));

        if(cmd.hasOption(PLOIDY_MIN))
            mPloidyMin = Double.parseDouble(cmd.getOptionValue(PLOIDY_MIN, "0"));

        if(cmd.hasOption(PLOIDY_MAX))
            mPloidyMax = Double.parseDouble(cmd.getOptionValue(PLOIDY_MAX, "0"));

        if(cmd.hasOption(SAMPLE_IDS))
        {
            mSampleIds.addAll(Arrays.stream(cmd.getOptionValue(SAMPLE_IDS).split(";")).collect(Collectors.toList()));
        }

        buildBucketMap();

        return true;
    }

    private void buildBucketMap()
    {
        char[] refBases = {'C', 'T'};
        char[] bases = {'A','C', 'G', 'T'};
        int index = 0;

        for(int i = 0; i < refBases.length; ++i)
        {
            char ref = refBases[i];

            for(int j = 0; j < bases.length; ++j)
            {
                char alt = bases[j];

                if(ref != alt)
                {
                    String baseChange = String.format("%c>%c", ref, alt);

                    for (int k = 0; k < bases.length; ++k)
                    {
                        char before = bases[k];

                        for (int l = 0; l < bases.length; ++l)
                        {
                            char after = bases[l];

                            String context = String.format("%c%c%c", before, ref, after);

                            String bucketName = baseChange + "_" + context;

                            mBucketStringToIndex.put(bucketName, index);
                            ++index;
                        }
                    }
                }
            }
        }
    }

    public void loadData()
    {
        List<String> sampleIds = null;

        if(!mSampleIds.isEmpty())
        {
            sampleIds = mSampleIds;
        }
        else
        {
            if(mApplySampleQC)
                sampleIds = mDbAccess.readPurpleSampleListPassingQC(MIN_SAMPLE_PURITY);
            else
                sampleIds = mDbAccess.readSomaticVariantSampleList();
        }

        mSampleBucketCounts = new SigMatrix(SNV_BUCKET_COUNT, sampleIds.size());

        LOGGER.info("retrieving SNV data for {} samples", sampleIds.size());

        for(int sampleIndex = 0; sampleIndex < sampleIds.size(); ++sampleIndex)
        {
            String sampleId = sampleIds.get(sampleIndex);
            final List<SomaticVariant> variants = mDbAccess.readSomaticVariants(sampleId);

            LOGGER.info("sample({}) processing {} variants", sampleId, variants.size());

            processSampleVariants(sampleId, variants, sampleIndex);
        }

        try
        {
            BufferedWriter writer = getNewFile(mOutputDir, mOutputFileId + "_sample_counts.csv");

            int i = 0;
            for(; i < sampleIds.size()-1; ++i)
            {
                writer.write(String.format("%s,", sampleIds.get(i)));
            }
            writer.write(String.format("%s", sampleIds.get(i)));

            writer.newLine();

            writeMatrixData(writer, mSampleBucketCounts, true);

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile: {}", e.toString());
        }
    }

    private void processSampleVariants(final String sampleId, List<SomaticVariant> variants, int sampleIndex)
    {
        double[][] sampleCounts = mSampleBucketCounts.getData();

        for(final SomaticVariant variant : variants)
        {
            if(variant.isFiltered() || !variant.isSnp())
                continue;

            if(variant.alt().length() != 1)
                continue;

            String rawContext = variant.trinucleotideContext();

            if(rawContext.contains("N"))
                continue;

            // check filters
            if(mSubclonalLikelihoodMin != null && variant.subclonalLikelihood() < mSubclonalLikelihoodMin)
                continue;

            if(mSubclonalLikelihoodMax != null && variant.subclonalLikelihood() > mSubclonalLikelihoodMax)
                continue;

            if(mPloidyMin != null && variant.ploidy() < mPloidyMin)
                continue;

            if(mPloidyMax != null && variant.ploidy() > mPloidyMax)
                continue;

            // convert base change to standard set and the context accordingly
            String baseChange;
            String context;
            if(variant.ref().charAt(0) == 'A' || variant.ref().charAt(0) == 'G')
            {
                baseChange = String.format("%c>%c", convertBase(variant.ref().charAt(0)), convertBase(variant.alt().charAt(0)));

                // convert the context as well
                context = String.format("%c%c%c",
                        convertBase(rawContext.charAt(2)), convertBase(rawContext.charAt(1)), convertBase(rawContext.charAt(0)));
            }
            else
            {
                baseChange = variant.ref() + ">" + variant.alt();
                context = rawContext;
            }

            String bucketName = baseChange + "_" + context;
            Integer bucketIndex = mBucketStringToIndex.get(bucketName);

            if(bucketIndex == null)
            {
                LOGGER.error("sample({}) invalid bucketName({}) from baseChange({} raw={}>{}) context({} raw={}",
                        sampleId, bucketName, baseChange, variant.ref(), variant.alt(), context, rawContext);

                return;
            }

            ++sampleCounts[bucketIndex][sampleIndex];
        }
    }

    private static char convertBase(char base)
    {
        if(base == 'A') return 'T';
        if(base == 'T') return 'A';
        if(base == 'C') return 'G';
        if(base == 'G') return 'C';
        return base;
    }

    private static String standardiseSnv(final String snv)
    {
        // convert to equivalent strand's base
        if(snv.equals("G>T")) return "C>A";
        if(snv.equals("G>C")) return "C>G";
        if(snv.equals("G>A")) return "C>T";
        if(snv.equals("A>T")) return "T>A";
        if(snv.equals("A>G")) return "T>C";
        if(snv.equals("A>C")) return "T>G";

        return snv;
    }


}

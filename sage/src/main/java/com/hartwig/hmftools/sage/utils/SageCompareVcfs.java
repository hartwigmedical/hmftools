package com.hartwig.hmftools.sage.utils;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sage.SageCommon.SAGE_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ALT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REF;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.QUAL;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsDouble;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LIST_SEPARATOR;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_MICROHOMOLOGY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_REPEAT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_REPEAT_SEQUENCE;
import static com.hartwig.hmftools.sage.SageCommon.APP_NAME;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.utils.VariantData.comparePositions;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_BASE_QUAL;
import static com.hartwig.hmftools.sage.vcf.VcfTags.MAX_READ_EDGE_DISTANCE;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_CORE;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_INFO;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_JITTER;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_LEFT_FLANK;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_RIGHT_FLANK;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VariantReadSupport;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.sage.vcf.ReadContextVcfInfo;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class SageCompareVcfs
{
    private final String mSampleId;
    private final String mOriginalVcf;
    private final String mNewVcf;
    private final String mOutputDir;
    private final String mOutputId;

    // comparison config
    private final double mDiffAbsPermitted;
    private final double mDiffPercPermitted;
    private final boolean mPassOnly;

    private final List<String> mExcludedFields;

    private final BufferedWriter mWriter;

    // counters and state
    private int mCompareCount;
    private int mDiffCount;
    private int mUnmatchedOrigCount;
    private int mUnmatchedNewCount;

    private static final String ORIGINAL_VCF = "original_vcf";
    private static final String NEW_VCF = "new_vcf";

    private static final String DIFF_ABS = "diff_abs";
    private static final String DIFF_PERC = "diff_perc";
    private static final String PASS_ONLY = "pass_only";
    private static final String EXCLUDE_FIELDS = "exclude_fields";

    private static final double DEFAULT_DIFF_PERC = 0.1;
    private static final double DEFAULT_DIFF_ABS = 2;

    public SageCompareVcfs(final ConfigBuilder configBuilder)
    {
        mSampleId = configBuilder.getValue(SAMPLE);
        mOriginalVcf = configBuilder.getValue(ORIGINAL_VCF);
        mNewVcf = configBuilder.getValue(NEW_VCF);
        mOutputDir = parseOutputDir(configBuilder);
        mOutputId = configBuilder.getValue(OUTPUT_ID);

        mDiffAbsPermitted = Double.parseDouble(configBuilder.getValue(DIFF_ABS, String.valueOf(DEFAULT_DIFF_ABS)));
        mDiffPercPermitted = Double.parseDouble(configBuilder.getValue(DIFF_PERC, String.valueOf(DEFAULT_DIFF_PERC)));
        mPassOnly = configBuilder.hasFlag(PASS_ONLY);

        mExcludedFields = Lists.newArrayList();

        if(configBuilder.hasValue(EXCLUDE_FIELDS))
        {
            Arrays.stream(configBuilder.getValue(EXCLUDE_FIELDS).split(ITEM_DELIM, -1)).forEach(x -> mExcludedFields.add(x));
        }

        mWriter = initialiseWriter();

        mDiffCount = 0;
        mCompareCount = 0;
        mUnmatchedOrigCount = 0;
        mUnmatchedNewCount = 0;
    }

    private static final int LOG_COUNT = 100000;

    public void run()
    {
        if(mOriginalVcf == null || mNewVcf == null)
        {
            SG_LOGGER.error("missing VCFs in config");
            System.exit(1);
        }

        VcfFileReader origReader = new VcfFileReader(mOriginalVcf);
        VcfFileReader newReader = new VcfFileReader(mNewVcf);

        if(!origReader.fileValid() || !newReader.fileValid())
        {
            SG_LOGGER.error("invalid VCF paths: original({}) and new({})", mOriginalVcf, mNewVcf);
            System.exit(1);
        }

        SG_LOGGER.info("comparing VCFs: orig({}) vs new({})", mOriginalVcf, mNewVcf);

        Iterator origIter = origReader.iterator();
        Iterator newIter = newReader.iterator();

        VariantData origVar = VariantData.fromContext((VariantContext)origIter.next());
        VariantData newVar = VariantData.fromContext((VariantContext)newIter.next());

        List<VariantData> origVariants = Lists.newArrayList();
        List<VariantData> newVariants = Lists.newArrayList();

        int nextLog = LOG_COUNT;

        while(newVar != null || origVar != null)
        {
            int totalComparisons = totalComparisons();

            if(totalComparisons >= nextLog)
            {
                nextLog += LOG_COUNT;
                SG_LOGGER.info("processed {} variants", totalComparisons);
            }

            if(!origVariants.isEmpty() && !newVariants.isEmpty())
            {
                boolean positionMatch = false;

                if(origVar != null && comparePositions(origVariants.get(0), origVar) == 0)
                {
                    origVariants.add(origVar);
                    origVar = getNextVariant(origIter);
                    positionMatch = true;
                }

                if(newVar != null && comparePositions(newVariants.get(0), newVar) == 0)
                {
                    newVariants.add(newVar);
                    positionMatch = true;
                    newVar = getNextVariant(newIter);
                }

                if(positionMatch)
                    continue; // keep looking for more at this position

                // run comparisons within this group
                compareVariants(origVariants, newVariants);
                continue;
            }

            if(newVar != null && origVar != null)
            {
                int posCompare = comparePositions(origVar, newVar);

                if(posCompare == 0)
                {
                    origVariants.add(origVar);
                    newVariants.add(newVar);

                    // check for others at this exact position

                    origVar = getNextVariant(origIter);
                    newVar = getNextVariant(newIter);
                }
                else if(posCompare < 0)
                {
                    // original SV has a lower positions
                    writeUnmatchedVariant(origVar, false);
                    origVar = getNextVariant(origIter);
                }
                else
                {
                    writeUnmatchedVariant(newVar, true);
                    newVar = getNextVariant(newIter);
                }
            }
            else if(newVar != null)
            {
                writeUnmatchedVariant(newVar, true);
                newVar = getNextVariant(newIter);
            }
            else if(origVar != null)
            {
                writeUnmatchedVariant(origVar, false);
                origVar = getNextVariant(origIter);
            }
        }

        SG_LOGGER.info("summary: total variants({}) compared({} diffs={}) unmatched(orig={} new={})",
                totalComparisons(), mCompareCount, mDiffCount, mUnmatchedOrigCount, mUnmatchedNewCount);

        closeBufferedWriter(mWriter);

        SG_LOGGER.info("Sage compare VCFs complete");
    }

    private static VariantData getNextVariant(final Iterator iter)
    {
        return VariantData.fromContext((VariantContext)iter.next());
    }

    private int totalComparisons() { return mCompareCount + mUnmatchedNewCount + mUnmatchedOrigCount; }

    private void compareVariants(final List<VariantData> origVariants, final List<VariantData> newVariants)
    {
        int origIndex = 0;
        while(origIndex < origVariants.size())
        {
            VariantData origVar = origVariants.get(origIndex);

            boolean matched = false;
            int newIndex = 0;
            while(newIndex < newVariants.size())
            {
                VariantData newVar = newVariants.get(newIndex);

                if(origVar.matches(newVar))
                {
                    compareVariants(origVar, newVar);
                    newVariants.remove(newIndex);
                    origVariants.remove(origIndex);
                    matched = true;
                    break;
                }

                ++newIndex;
            }

            if(!matched)
                ++origIndex;
        }

        origVariants.forEach(x -> writeUnmatchedVariant(x, false));
        newVariants.forEach(x -> writeUnmatchedVariant(x, true));

        origVariants.clear();
        newVariants.clear();
    }

    private void compareVariants(final VariantData origVar, final VariantData newVar)
    {
        if(mPassOnly)
        {
            if(!origVar.isPassing() && !newVar.isPassing())
                return;

            if(!origVar.isPassing())
            {
                writeUnmatchedVariant(newVar, true);
                return;
            }
            else if(!newVar.isPassing())
            {
                writeUnmatchedVariant(origVar, false);
                return;
            }
        }

        mCompareCount++;

        // compare each field in turn
        compareField(origVar, newVar, QUAL, origVar.qual(), newVar.qual());
        compareReadBases(origVar, newVar);
        compareAttributeField(origVar, newVar, READ_CONTEXT_JITTER, FieldType.STRING); // will still work for int arrays but could change
        compareAttributeField(origVar, newVar, READ_CONTEXT_MICROHOMOLOGY, FieldType.STRING);
        compareAttributeField(origVar, newVar, READ_CONTEXT_REPEAT_SEQUENCE, FieldType.STRING);
        compareAttributeField(origVar, newVar, READ_CONTEXT_REPEAT_COUNT, FieldType.EXACT_INT);

        compareAttributeField(origVar, newVar, AVG_BASE_QUAL, FieldType.DECIMAL);
        compareAttributeField(origVar, newVar, MAX_READ_EDGE_DISTANCE, FieldType.DECIMAL);

        // compare genotype fields
        for(int i = 0; i < origVar.context().getGenotypes().size(); ++i)
        {
            Genotype origGenotype = origVar.context().getGenotypes().get(i);
            Genotype newGenotype = newVar.context().getGenotypes().get(i);
            // compareGenotypeAttributeField(origVar, newVar, origGenotype, newGenotype, RAW_SUPPORT_DEPTH);

            compareMatchCountsField(origVar, newVar, origGenotype, newGenotype);
        }

        // compare filters

        Set<String> origFilters = origVar.filters();
        Set<String> newFilters = newVar.filters();

        // first check for a difference in PASS vs not
        Set<String> origFilterDiffs = origFilters.stream().filter(x -> !newFilters.contains(x)).collect(Collectors.toSet());
        Set<String> newFilterDiffs = newFilters.stream().filter(x -> !origFilters.contains(x)).collect(Collectors.toSet());

        if(!newFilterDiffs.isEmpty() || !origFilterDiffs.isEmpty())
        {
            String origFiltersStr = origVar.filters().isEmpty() ? PASS : filtersStr(origFilterDiffs);
            String newFiltersStr = newVar.filters().isEmpty() ? PASS : filtersStr(newFilterDiffs);
            writeDiffs(origVar, newVar, "FILTERS", origFiltersStr, newFiltersStr);
        }

        // compare allelic depth
        int origAd = origVar.allelicDepth();
        int newAd = newVar.allelicDepth();

        if(hasValueDiff(origAd, newAd, mDiffAbsPermitted, mDiffPercPermitted))
        {
            writeDiffs(origVar, newVar, "ALLELE_DEPTH", String.valueOf(origAd), String.valueOf(newAd));
            return;
        }

        boolean origPhased = origVar.context().hasAttribute(LOCAL_PHASE_SET);
        boolean newPhased = origVar.context().hasAttribute(LOCAL_PHASE_SET);

        if(origPhased != newPhased)
        {
            writeDiffs(origVar, newVar, "PHASED", String.valueOf(origPhased), String.valueOf(newPhased));
            return;
        }
    }

    private enum FieldType
    {
        STRING,
        DECIMAL,
        EXACT_INT;
    }

    private void compareAttributeField(
            final VariantData origVar, final VariantData newVar, final String vcfTag, final FieldType fieldType)
    {
        if(fieldType == FieldType.STRING)
        {
            String origValue = origVar.context().getAttributeAsString(vcfTag, "");
            String newValue = newVar.context().getAttributeAsString(vcfTag, "");
            compareField(origVar, newVar, vcfTag, origValue, newValue);
        }
        else if(fieldType == FieldType.EXACT_INT)
        {
            int origValue = origVar.context().getAttributeAsInt(vcfTag, 0);
            int newValue = newVar.context().getAttributeAsInt(vcfTag, 0);

            compareField(origVar, newVar, vcfTag, origValue, newValue);
        }
        else
        {
            double origValue = origVar.context().getAttributeAsDouble(vcfTag, 0);
            double newValue = origVar.context().getAttributeAsDouble(vcfTag, 0);
            compareField(origVar, newVar, vcfTag, origValue, newValue);
        }
    }

    private void compareReadBases(final VariantData origVar, final VariantData newVar)
    {
        String oldReadBases = extractReadBases(origVar);
        String newReadBases = extractReadBases(newVar);

        compareField(origVar, newVar, READ_CONTEXT_INFO, oldReadBases, newReadBases);
    }

    private static String extractReadBases(final VariantData var)
    {
        if(var.context().hasAttribute(READ_CONTEXT_INFO))
        {
            return ReadContextVcfInfo.fromVcfTag(var.context().getAttributeAsString(READ_CONTEXT_INFO, "")).readBases();
        }
        else
        {
            return var.context().getAttributeAsString(READ_CONTEXT_LEFT_FLANK, "")
                    + var.context().getAttributeAsString(READ_CONTEXT_CORE, "")
                    + var.context().getAttributeAsString(READ_CONTEXT_RIGHT_FLANK, "");
        }
    }

    private void compareGenotypeAttributeField(
            final VariantData origVar, final VariantData newVar, final Genotype origGenotype, final Genotype newGenotype, final String vcfTag)
    {
        double origValue = getGenotypeAttributeAsDouble(origGenotype, vcfTag, 0);
        double newValue = getGenotypeAttributeAsDouble(newGenotype, vcfTag, 0);
        compareField(origVar, newVar, format("%s:%s", origGenotype.getSampleName(), vcfTag), origValue, newValue);
    }

    private void compareMatchCountsField(
            final VariantData origVar, final VariantData newVar, final Genotype origGenotype, final Genotype newGenotype)
    {
        final int[] origQualCounts = parseIntegerList(origGenotype, READ_CONTEXT_COUNT);

        final int[] newQualCounts = parseIntegerList(newGenotype, READ_CONTEXT_COUNT);

        if(origQualCounts.length == newQualCounts.length)
        {
            for(int i = 0; i < VariantReadSupport.values().length; ++i)
            {
                String fieldName = format("%s:%s", origGenotype.getSampleName(), VariantReadSupport.values()[i]);

                compareField(origVar, newVar, format("%s:%s", origGenotype.getSampleName(), fieldName), origQualCounts[i], newQualCounts[i]);
            }
        }
        else if(origQualCounts.length == 7 && newQualCounts.length == VariantReadSupport.values().length)
        {
            compareField(origVar, newVar, format("%s:%s", origGenotype.getSampleName(), VariantReadSupport.FULL),
                    origQualCounts[0] + origQualCounts[1], newQualCounts[0] + newQualCounts[1]);

            compareField(origVar, newVar, format("%s:%s", origGenotype.getSampleName(), VariantReadSupport.CORE), origQualCounts[2], newQualCounts[2]);

            compareField(origVar, newVar, format("%s:%s", origGenotype.getSampleName(), VariantReadSupport.REALIGNED), origQualCounts[3], newQualCounts[3]);

            compareField(origVar, newVar, format("%s:%s", origGenotype.getSampleName(), VariantReadSupport.REF), origQualCounts[5], newQualCounts[4]);
        }
    }

    private static int[] parseIntegerList(final Genotype genotype, final String vcfTag)
    {
        final String[] stringValues = genotype.getExtendedAttribute(vcfTag, 0).toString().split(LIST_SEPARATOR, -1);
        int[] values = new int[stringValues.length];

        for(int i = 0; i < stringValues.length; ++i)
        {
            values[i] = Integer.parseInt(stringValues[i]);
        }

        return values;
    }


    private void compareField(
            final VariantData origVar, final VariantData newVar, final String vcfTag, final String origValue, final String newValue)
    {
        if(!origValue.equals(newValue))
            writeDiffs(origVar, newVar, vcfTag, origValue, newValue);
    }

    private void compareField(
            final VariantData origVar, final VariantData newVar, final String vcfTag, final int origValue, final int newValue)
    {
        if(origValue != newValue)
            writeDiffs(origVar, newVar, vcfTag, String.valueOf(origValue), String.valueOf(newValue));
    }

    private void compareField(
            final VariantData origVar, final VariantData newVar, final String vcfTag, final double origValue, final double newValue)
    {
        if(hasValueDiff(origValue, newValue, mDiffAbsPermitted, mDiffPercPermitted))
            writeDiffs(origVar, newVar, vcfTag, formatNumber(vcfTag, origValue), formatNumber(vcfTag, newValue));
    }

    private static final List<String> DECIMAL_TAGS = List.of("");

    private static String formatNumber(final String vcfTag, final double value)
    {
        if(DECIMAL_TAGS.contains(vcfTag))
            return String.valueOf(value);
        else
            return format("%.0f", value);
    }

    private static boolean hasValueDiff(final double value1, final double value2, final double diffAbs, final double diffPerc)
    {
        if(value1 == 0 && value2 == 0)
            return false;

        double diff = abs(value1 - value2);
        return diff > diffAbs && diff / max(value1, value2) > diffPerc;
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            String fileName = mOutputDir + mSampleId + SAGE_FILE_ID + ".compare";

            if(mOutputId != null)
                fileName += "." + mOutputId;

            fileName += TSV_EXTENSION;

            SG_LOGGER.info("writing comparison file: {}", fileName);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add("VariantInfo");
            sj.add(FLD_CHROMOSOME).add(FLD_POSITION).add(FLD_REF).add(FLD_ALT).add("Tier");

            sj.add("DiffType").add("OrigValue").add("NewValue").add("SharedFilters");
            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    private void writeUnmatchedVariant(final VariantData var, boolean isNew)
    {
        if(mPassOnly && !var.isPassing())
            return;

        if(isNew)
        {
            mUnmatchedNewCount++;
            writeDiffs(null, var, "NO_ORIG", "", "");
        }
        else
        {
            mUnmatchedOrigCount++;
            writeDiffs(var, null, "NO_NEW", "", "");
        }
    }

    private void writeDiffs(
            final VariantData origVar, final VariantData newVar, final String fieldName, final String origValue, final String newValue)
    {
        if(mExcludedFields.contains(fieldName))
            return;

        if(origVar != null && newVar != null)
        {
            ++mDiffCount;
        }

        try
        {
            VariantData var = origVar != null ? origVar : newVar;

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(format("%s:%d %s>%s", var.Chromosome, var.Position, var.Ref, var.Alt));
            sj.add(var.Chromosome).add(String.valueOf(var.Position)).add(var.Ref).add(var.Alt).add(var.tier().toString());

            sj.add(fieldName).add(origValue).add(newValue);

            /*
            Set<String> sharedFilters = Sets.newHashSet();

            if(origVar != null)
            {
                origVar.filters().forEach(x -> sharedFilters.add(x));
            }

            if(newVar != null)
            {
                newVar.filters().forEach(x -> sharedFilters.add(x));
            }

            if(sharedFilters.isEmpty())
                sharedFilters.add(PASS);

            mWriter.write(format("\t%.0f\t%.0f\t%s\t%d",
                    origVar != null ? origVar.qual() : -1, newVar != null ? newVar.qual() : -1, filtersStr(sharedFilters), maxDepth));
             */

            mWriter.write(sj.toString());
            mWriter.newLine();
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to write output file: {}", e.toString());
        }
    }

    private String filtersStr(final Set<String> filters)
    {
        StringJoiner sj = new StringJoiner(";");
        filters.forEach(x -> sj.add(x));
        return sj.toString();
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addPath(ORIGINAL_VCF, true, "Original VCF file");
        configBuilder.addPath(NEW_VCF, true, "New VCF file");
        configBuilder.addFlag(PASS_ONLY, "Only compare passing variants");
        configBuilder.addDecimal(DIFF_ABS, "Absolute value difference", DEFAULT_DIFF_ABS);
        configBuilder.addDecimal(DIFF_PERC, "Percentage value difference", DEFAULT_DIFF_PERC);
        configBuilder.addConfigItem(EXCLUDE_FIELDS, false, "List of VCF fields to ignore, separated by ';'");

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        SageCompareVcfs sageCompare = new SageCompareVcfs(configBuilder);
        sageCompare.run();
    }
}

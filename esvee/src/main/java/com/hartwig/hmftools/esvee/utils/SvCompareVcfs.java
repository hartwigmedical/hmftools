package com.hartwig.hmftools.esvee.utils;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.SvVcfTags.PON_FILTER_PON;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.STRAND_BIAS;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.GenotypeIds.fromVcfHeader;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.esvee.caller.HotspotCache;
import com.hartwig.hmftools.esvee.caller.TargetRegions;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class SvCompareVcfs
{
    private final String mSampleId;
    private final String mReferenceId;
    private final String mOriginalVcf;
    private final String mNewVcf;
    private final String mOutputDir;
    private final String mOutputId;

    private final Map<String,List<VariantContext>> mOrigChrBreakendMap;
    private final Map<String,List<VariantContext>> mNewChrBreakendMap;

    private final BufferedWriter mWriter;

    private final boolean mIgnorePonDiff;
    private final boolean mWriteAllDiffs;

    private final List<VcfCompareField> mVcfCheckFields;

    private static final String ORIGINAL_VCF = "original_vcf";
    private static final String NEW_VCF = "new_vcf";
    private static final String IGNORE_PON_DIFF = "ignore_pon_diff";
    private static final String WRITE_ALL_DIFFS = "write_all_diffs";

    private static final int DEFAULT_MAX_DIFF = 20;
    private static final double DEFAULT_MAX_DIFF_PERC = 0.2;

    public SvCompareVcfs(final ConfigBuilder configBuilder)
    {
        mSampleId = configBuilder.getValue(SAMPLE);
        mReferenceId = configBuilder.getValue(REFERENCE, "");
        mOriginalVcf = configBuilder.getValue(ORIGINAL_VCF);
        mNewVcf = configBuilder.getValue(NEW_VCF);
        mOutputDir = parseOutputDir(configBuilder);
        mOutputId = configBuilder.getValue(OUTPUT_ID);

        mOrigChrBreakendMap = Maps.newHashMap();
        mNewChrBreakendMap = Maps.newHashMap();

        mIgnorePonDiff = configBuilder.hasFlag(IGNORE_PON_DIFF);
        mWriteAllDiffs = configBuilder.hasFlag(WRITE_ALL_DIFFS);

        mVcfCheckFields = Lists.newArrayList();
        addComparisonFields();

        mWriter = initialiseWriter();
    }

    private void addComparisonFields()
    {
        mVcfCheckFields.add(new VcfCompareField(REF_DEPTH, GenotypeScope.BOTH, VariantTypeScope.BREAKEND, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));
        mVcfCheckFields.add(new VcfCompareField(REF_DEPTH_PAIR, GenotypeScope.BOTH, VariantTypeScope.BREAKEND, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));

        mVcfCheckFields.add(new VcfCompareField(STRAND_BIAS, GenotypeScope.COMBINED, VariantTypeScope.BREAKEND, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));

        /*
        mVcfCheckFields.add(new VcfCompareField(SPLIT_READS, GenotypeScope.BOTH, VariantTypeScope.BOTH, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));
        mVcfCheckFields.add(new VcfCompareField(INDEL_COUNT, GenotypeScope.BOTH, VariantTypeScope.BOTH, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));

        // SV only
        mVcfCheckFields.add(new VcfCompareField(QUAL, GenotypeScope.BOTH, VariantTypeScope.SV, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));
        mVcfCheckFields.add(new VcfCompareField(SV_FRAG_COUNT, GenotypeScope.BOTH, VariantTypeScope.SV, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));
        mVcfCheckFields.add(new VcfCompareField(READ_PAIRS, GenotypeScope.BOTH, VariantTypeScope.SV, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));
        mVcfCheckFields.add(new VcfCompareField(GRIDSS_ASRP, GenotypeScope.BOTH, VariantTypeScope.SV, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));

        // assembly
        mVcfCheckFields.add(new VcfCompareField(GRIDSS_AS, GenotypeScope.COMBINED, VariantTypeScope.BOTH, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));
        mVcfCheckFields.add(new VcfCompareField(GRIDSS_CAS, GenotypeScope.COMBINED, VariantTypeScope.BOTH, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));
        mVcfCheckFields.add(new VcfCompareField(GRIDSS_RAS, GenotypeScope.COMBINED, VariantTypeScope.BOTH, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));

        // SGL only
        mVcfCheckFields.add(new VcfCompareField(SGL_FRAG_COUNT, GenotypeScope.BOTH, VariantTypeScope.SGL, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));
        mVcfCheckFields.add(new VcfCompareField(GRIDSS_BQ, GenotypeScope.BOTH, VariantTypeScope.SGL, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));
        mVcfCheckFields.add(new VcfCompareField(GRIDSS_BAQ, GenotypeScope.BOTH, VariantTypeScope.SGL, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));
        mVcfCheckFields.add(new VcfCompareField(GRIDSS_BSC, GenotypeScope.BOTH, VariantTypeScope.SGL, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));
        mVcfCheckFields.add(new VcfCompareField(GRIDSS_BASRP, GenotypeScope.BOTH, VariantTypeScope.SGL, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));
        mVcfCheckFields.add(new VcfCompareField(GRIDSS_BASSR, GenotypeScope.BOTH, VariantTypeScope.SGL, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));
        */
    }

    public void run()
    {
        if(mOriginalVcf == null || mNewVcf == null)
        {
            SV_LOGGER.error("missing VCFs");
            return;
        }

        SV_LOGGER.info("loading original VCF({})", mOriginalVcf);
        loadOriginalVariants(mOriginalVcf);

        compareVariants(mNewVcf);
        closeBufferedWriter(mWriter);

        SV_LOGGER.info("Gripss compare VCFs complete");
    }

    private void loadOriginalVariants(final String vcfFile)
    {
        VcfFileReader reader = new VcfFileReader(vcfFile);

        VCFHeader vcfHeader = reader.vcfHeader();
        GenotypeIds genotypeIds = fromVcfHeader(vcfHeader, mReferenceId, mSampleId);

        if(genotypeIds == null)
        {
            System.exit(1);
        }

        SV_LOGGER.info("genetype info: ref({}: {}) tumor({}: {})",
                genotypeIds.ReferenceOrdinal, genotypeIds.ReferenceId, genotypeIds.TumorOrdinal, genotypeIds.TumorId);


        String currentChr = "";
        List<VariantContext> breakends = null;

        for(VariantContext variantContext : reader.iterator())
        {
            String chromosome = variantContext.getContig();

            if(!currentChr.equals(chromosome))
            {
                currentChr = chromosome;
                breakends = Lists.newArrayList();
                mOrigChrBreakendMap.put(chromosome, breakends);
            }

            breakends.add(variantContext);
        }

        SV_LOGGER.info("loaded {} original SVs, incomplete({})",
                mOrigChrBreakendMap.values().stream().mapToInt(x -> x.size()).sum());
    }

    private void compareVariants(final String newVcfFile)
    {
        SV_LOGGER.info("loading new VCF({})", newVcfFile);

        VcfFileReader reader = new VcfFileReader(newVcfFile);

        VCFHeader vcfHeader = reader.vcfHeader();
        GenotypeIds genotypeIds = fromVcfHeader(vcfHeader, mReferenceId, mSampleId);

        if(genotypeIds == null)
            System.exit(1);

        int diffCount = 0;

        try
        {
            int newSvCount = 0;
            String currentChr = "";
            List<VariantContext> origBreakends = null;
            VariantContext origBreakend = null;
            int origBreakendIndex = 0;

            for(VariantContext newBreakend : reader.iterator())
            {
                ++newSvCount;

                if((newSvCount % 10000 == 0))
                {
                    SV_LOGGER.debug("processed {} variants", newSvCount);
                }

                if(!currentChr.equals(newBreakend.getContig()))
                {
                    currentChr = newBreakend.getContig();
                    origBreakends = mOrigChrBreakendMap.get(currentChr);
                    origBreakendIndex = 0;
                    origBreakend = origBreakends != null ? origBreakends.get(origBreakendIndex) : null;
                }

                // move ahead as required
                while(origBreakend != null)
                {
                    // test positions factoring in homology
                    break;
                }

                if(origBreakend == null)
                {
                    writeDiffs(null, newBreakend, "NO_ORIG", "", "");
                    ++diffCount;
                    continue;
                }

                // mOriginalSvData.remove(origSv.id());

                Set<String> origFilters = origBreakend.getFilters();
                Set<String> newFilters = newBreakend.getFilters();

                // first check for a difference in PASS vs not
                Set<String> origFilterDiffs = origFilters.stream().filter(x -> !newFilters.contains(x)).collect(Collectors.toSet());
                Set<String> newFilterDiffs = newFilters.stream().filter(x -> !origFilters.contains(x)).collect(Collectors.toSet());

                boolean ignorePonDiff = mIgnorePonDiff
                        && ((origFilterDiffs.isEmpty() && newFilterDiffs.size() == 1 && newFilterDiffs.contains(PON_FILTER_PON))
                        || (newFilterDiffs.isEmpty() && origFilterDiffs.size() == 1 && origFilterDiffs.contains(PON_FILTER_PON)));

                /*
                if(!ignorePonDiff && (!newFilterDiffs.isEmpty() || !origFilterDiffs.isEmpty()))
                {
                    boolean origIsPass = origStart.Context.isNotFiltered() || (origFilters.size() == 1 && origFilters.contains(PASS));
                    boolean newIsPass = newStart.Context.isNotFiltered() || (newFilters.size() == 1 && newFilters.contains(PASS));

                    if(origIsPass != newIsPass)
                    {
                        writeDiffs(
                                origSv, newSv, "FILTER_PASS",
                                filtersStr(origFilterDiffs, false), filtersStr(newFilterDiffs, false));
                    }
                    else
                    {
                        writeDiffs(
                                origSv, newSv, "FILTER_DIFF",
                                filtersStr(origFilterDiffs, false), filtersStr(newFilterDiffs, false));
                    }

                 */

                    ++diffCount;
                    continue;
                }

            /*
                // check specified VCF tags
                for(VcfCompareField compareField : mVcfCheckFields)
                {
                    checkVcfFieldDiff(origSv, newSv, compareField);
                }

                if(!mGridssDiffsOnly)
                {
                    // check local and remote linked by for assembled links
                    boolean origHasStartAssembled = origStart.Context.getAttributeAsString(LOCAL_LINKED_BY, "").contains("asm");
                    boolean newHasStartAssembled = newStart.Context.getAttributeAsString(LOCAL_LINKED_BY, "").contains("asm");
                    boolean origHasEndAssembled =
                            !origSv.isSgl() && origStart.Context.getAttributeAsString(REMOTE_LINKED_BY, "").contains("asm");
                    boolean newHasEndAssembled =
                            !origSv.isSgl() && newStart.Context.getAttributeAsString(REMOTE_LINKED_BY, "").contains("asm");

                    if(origHasStartAssembled != newHasStartAssembled || origHasEndAssembled != newHasEndAssembled)
                    {
                        writeDiffs(
                                origSv, newSv, "ASSEMBLY",
                                format("%s_%s", origHasStartAssembled, origHasEndAssembled),
                                format("%s_%s", newHasStartAssembled, newHasEndAssembled));
                        ++diffCount;
                        continue;
                    }
                }
            }
            */

            SV_LOGGER.info("loaded {} new SVs", newSvCount);
        }
        catch(Exception e)
        {
            SV_LOGGER.error("error reading vcf({}): {}", newVcfFile, e.toString());
        }

        /*
        for(SvData origSv : mOriginalSvData.values())
        {
            ++diffCount;
            writeDiffs(origSv, null, "NO_NEW", "", "");
        }
        */

        SV_LOGGER.info("diffTotal({})", diffCount);
    }

    private boolean checkVcfFieldDiff(final VariantContext origBreakend, final VariantContext newBreakend, final VcfCompareField compareField)
    {
        boolean hasDiff = false;

        /*
        if(compareField.TypeScope == VariantTypeScope.SV && origBreakend.isSgl())
            return false;
        if(compareField.TypeScope == VariantTypeScope.SGL && !origBreakend.isSgl())
            return false;

        if(compareField.Scope == GenotypeScope.COMBINED)
        {
            double origValue = origBreakend.contextStart().getAttributeAsDouble(compareField.VcfTag, 0);
            double newValue = newBreakend.contextStart().getAttributeAsDouble(compareField.VcfTag, 0);

            if(hasDiff(origValue, newValue, compareField.DiffAbs, compareField.DiffPerc))
            {
                hasDiff = true;
                writeDiffs(origBreakend, newBreakend, compareField.VcfTag, String.valueOf(origValue), String.valueOf(newValue));
            }
        }
        else
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                if(se == SE_END && origBreakend.isSgl())
                    continue;

                Breakend origBreakend = origBreakend.breakends()[se];
                Breakend newBreakend = newBreakend.breakends()[se];

                if(!mReferenceId.isEmpty() && (compareField.Scope == GenotypeScope.NORMAL || compareField.Scope == GenotypeScope.BOTH))
                {
                    double origValue = getGenotypeAttributeAsDouble(origBreakend.RefGenotype, compareField.VcfTag, 0);
                    double newValue = getGenotypeAttributeAsDouble(newBreakend.RefGenotype, compareField.VcfTag, 0);

                    if(hasDiff(origValue, newValue, compareField.DiffAbs, compareField.DiffPerc))
                    {
                        hasDiff = true;
                        writeDiffs(origBreakend, newBreakend, format("REF_%s", compareField.VcfTag), String.valueOf(origValue), String.valueOf(newValue));
                    }
                }

                if(compareField.Scope == GenotypeScope.TUMOR || compareField.Scope == GenotypeScope.BOTH)
                {
                    double origValue = getGenotypeAttributeAsDouble(origBreakend.TumorGenotype, compareField.VcfTag, 0);
                    double newValue = getGenotypeAttributeAsDouble(newBreakend.TumorGenotype, compareField.VcfTag, 0);

                    if(hasDiff(origValue, newValue, compareField.DiffAbs, compareField.DiffPerc))
                    {
                        hasDiff = true;
                        writeDiffs(origBreakend, newBreakend, format("TUMOR_%s", compareField.VcfTag), String.valueOf(origValue), String.valueOf(newValue));
                    }
                }
            }
        }
        */

        return hasDiff;
    }

    private static boolean hasDiff(double value1, double value2)
    {
        return hasDiff(value1, value2, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC);
    }

    private static boolean hasDiff(double value1, double value2, double maxDiff, double maxDiffPerc)
    {
        double diff = abs(value1 - value2);
        double diffPerc = diff / (double)max(value1, value2);
        return diff > maxDiff && diffPerc > maxDiffPerc;
    }

    private String filtersStr(final Set<String> filters, boolean setPassOnEmpty)
    {
        if(setPassOnEmpty && filters.isEmpty())
            return PASS;

        StringJoiner sj = new StringJoiner(";");
        filters.forEach(x -> sj.add(x));
        return sj.toString();
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            String fileName = mOutputDir + mSampleId + ".sv_compare";

            if(mOutputId != null)
                fileName += "." + mOutputId;

            fileName += TSV_EXTENSION;

            SV_LOGGER.info("writing comparison file: {}", fileName);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("OrigId\tNewId");

            writer.write("\tCoords\tType\tDiffType\tOrigValue\tNewValue\tOrigQual\tNewQual\tOrigFilters\tNewFilters");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    private void writeDiffs(final VariantContext origBreakend, final VariantContext newBreakend, final String diffType, final String origValue, final String newValue)
    {
        try
        {
            /*
            String coords = origBreakend != null ? makeSvCoords(origBreakend) : makeSvCoords(newBreakend);
            StructuralVariantType type = origBreakend != null ? origBreakend.type() : newBreakend.type();

            mWriter.write(String.format("%s\t%s",
                    origBreakend != null ? origBreakend.id() : "", newBreakend != null ? newBreakend.id() : ""));

            mWriter.write(format("\t%s\t%s\t%s\t%s\t%s",
                    coords, type, diffType, origValue, newValue));

            mWriter.write(String.format("\t%.1f\t%.1f\t%s\t%s",
                    origBreakend != null ? origBreakend.breakendStart().Qual : -1, newBreakend != null ? newBreakend.breakendStart().Qual : -1,
                    origBreakend != null ? filtersStr(origBreakend.breakendStart().Context.getFilters(), true) : "",
                    newBreakend != null ? filtersStr(newBreakend.breakendStart().Context.getFilters(), true) : ""));

            */
            mWriter.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write output file: {}", e.toString());
        }
    }

    private static String makeSvCoords(final VariantContext breakend)
    {
        /*
        if(sv.isSgl())
        {
            return String.format("%s:%d:%d", sv.chromosomeStart(), sv.posStart(), sv.orientStart());
        }
        else
        {
            return String.format("%s:%d:%d-%s:%d:%d",
                    sv.chromosomeStart(), sv.posStart(), sv.orientStart(), sv.chromosomeEnd(), sv.posEnd(), sv.orientEnd());
        }
        */

        return null;
    }

    private enum GenotypeScope
    {
        NORMAL,
        TUMOR,
        BOTH,
        COMBINED;
    }

    private enum VariantTypeScope
    {
        BREAKEND,
        VARIANT;
    }

    private class VcfCompareField
    {
        public final String VcfTag;
        public final GenotypeScope Scope;
        public final double DiffAbs;
        public final double DiffPerc;
        public final VariantTypeScope TypeScope;

        public VcfCompareField(
                final String vcfTag, final GenotypeScope scope, final VariantTypeScope typeScope, final double diffAbs, final double diffPerc)
        {
            VcfTag = vcfTag;
            Scope = scope;
            TypeScope = typeScope;
            DiffAbs = diffAbs;
            DiffPerc = diffPerc;
        }

        public String toString() { return format("tag(%s) scope(%s) st(%s)", VcfTag, Scope, TypeScope); }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();

        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addConfigItem(REFERENCE, false, REFERENCE_DESC);
        configBuilder.addPath(ORIGINAL_VCF, true, "Optional, name of the reference sample");
        configBuilder.addPath(NEW_VCF, true, "Path to the GRIDSS structural variant VCF file");
        configBuilder.addFlag(IGNORE_PON_DIFF, "Ignore diffs if just PON filter");
        configBuilder.addFlag(WRITE_ALL_DIFFS, "Write all VCF field diffs, not just the first");

        // CHECK: required?
        HotspotCache.addConfig(configBuilder);
        TargetRegions.addConfig(configBuilder);

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);
        setLogLevel(configBuilder);

        SvCompareVcfs gripssCompare = new SvCompareVcfs(configBuilder);
        gripssCompare.run();
    }
}
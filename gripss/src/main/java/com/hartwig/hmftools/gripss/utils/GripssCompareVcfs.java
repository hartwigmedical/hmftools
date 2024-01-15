package com.hartwig.hmftools.gripss.utils;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.sv.SvVcfTags.GRIDSS_AS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.GRIDSS_ASRP;
import static com.hartwig.hmftools.common.sv.SvVcfTags.GRIDSS_BAQ;
import static com.hartwig.hmftools.common.sv.SvVcfTags.GRIDSS_BASRP;
import static com.hartwig.hmftools.common.sv.SvVcfTags.GRIDSS_BASSR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.GRIDSS_BQ;
import static com.hartwig.hmftools.common.sv.SvVcfTags.GRIDSS_BSC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.GRIDSS_CAS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.GRIDSS_RAS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.INDEL_COUNT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.LOCAL_LINKED_BY;
import static com.hartwig.hmftools.common.sv.SvVcfTags.PON_FILTER_PON;
import static com.hartwig.hmftools.common.sv.SvVcfTags.QUAL;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REMOTE_LINKED_BY;
import static com.hartwig.hmftools.common.sv.SvVcfTags.READ_PAIRS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SGL_FRAG_COUNT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SPLIT_READS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.STRAND_BIAS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SV_FRAG_COUNT;
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
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsDouble;
import static com.hartwig.hmftools.gripss.GripssConfig.GR_LOGGER;
import static com.hartwig.hmftools.common.variant.GenotypeIds.fromVcfHeader;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.gripss.VariantBuilder;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.filters.HotspotCache;
import com.hartwig.hmftools.gripss.filters.TargetRegions;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class GripssCompareVcfs
{
    private final String mSampleId;
    private final String mReferenceId;
    private final String mOriginalVcf;
    private final String mNewVcf;
    private final String mOutputDir;
    private final String mOutputId;

    private final VariantBuilder mVariantBuilder;

    private final Map<String,SvData> mOriginalSvData;
    private final Map<String, List<SvData>> mOriginalCoordsSvData; // keyed by chromosome pair

    private final BufferedWriter mWriter;

    private final boolean mIgnorePonDiff;
    private final boolean mKeyByCoords; // instead of assuming VCF Ids match
    private final boolean mWriteAllDiffs;
    private final boolean mGridssDiffsOnly;
    private final boolean mRefDepthDiffsOnly;

    private final List<VcfCompareField> mVcfCheckFields;

    private static final String ORIGINAL_VCF = "original_vcf";
    private static final String NEW_VCF = "new_vcf";
    private static final String IGNORE_PON_DIFF = "ignore_pon_diff";
    private static final String KEY_BY_COORDS = "key_by_coords";
    private static final String WRITE_ALL_DIFFS = "write_all_diffs";

    private static final String GRIDSS_ONLY = "gridss_only";
    private static final String REF_DEPTH_ONLY = "ref_depth_only";

    private static final int DEFAULT_MAX_DIFF = 20;
    private static final double DEFAULT_MAX_DIFF_PERC = 0.2;

    public GripssCompareVcfs(final ConfigBuilder configBuilder)
    {
        mSampleId = configBuilder.getValue(SAMPLE);
        mReferenceId = configBuilder.getValue(REFERENCE, "");
        mOriginalVcf = configBuilder.getValue(ORIGINAL_VCF);
        mNewVcf = configBuilder.getValue(NEW_VCF);
        mOutputDir = parseOutputDir(configBuilder);
        mOutputId = configBuilder.getValue(OUTPUT_ID);

        mOriginalSvData = Maps.newHashMap();
        mOriginalCoordsSvData = Maps.newHashMap();

        mVariantBuilder = new VariantBuilder(
                null, new HotspotCache(configBuilder), new TargetRegions(configBuilder), false);

        mIgnorePonDiff = configBuilder.hasFlag(IGNORE_PON_DIFF);
        mKeyByCoords = configBuilder.hasFlag(KEY_BY_COORDS);
        mWriteAllDiffs = configBuilder.hasFlag(WRITE_ALL_DIFFS);
        mGridssDiffsOnly = configBuilder.hasFlag(GRIDSS_ONLY);
        mRefDepthDiffsOnly = configBuilder.hasFlag(REF_DEPTH_ONLY);

        mVcfCheckFields = Lists.newArrayList();
        addComparisonFields();

        mWriter = initialiseWriter();
    }

    private void addComparisonFields()
    {
        if(mRefDepthDiffsOnly)
        {
            mVcfCheckFields.add(new VcfCompareField(REF_DEPTH, GenotypeScope.BOTH, VariantTypeScope.BOTH, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));
            mVcfCheckFields.add(new VcfCompareField(REF_DEPTH_PAIR, GenotypeScope.BOTH, VariantTypeScope.BOTH, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));
            return;
        }

        if(!mGridssDiffsOnly)
        {
            mVcfCheckFields.add(new VcfCompareField(REF_DEPTH, GenotypeScope.BOTH, VariantTypeScope.BOTH, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));
            mVcfCheckFields.add(new VcfCompareField(REF_DEPTH_PAIR, GenotypeScope.BOTH, VariantTypeScope.BOTH, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));
        }

        mVcfCheckFields.add(new VcfCompareField(STRAND_BIAS, GenotypeScope.COMBINED, VariantTypeScope.BOTH, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));
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
    }

    public void run()
    {
        if(mOriginalVcf == null || mNewVcf == null)
        {
            GR_LOGGER.error("missing VCFs");
            return;
        }

        GR_LOGGER.info("loading original VCF({})", mOriginalVcf);
        loadOriginalVariants(mOriginalVcf);

        GR_LOGGER.info("loaded {} original SVs, incomplete({})", mOriginalSvData.size(), mVariantBuilder.incompleteSVs());

        compareVariants(mNewVcf);
        closeBufferedWriter(mWriter);

        GR_LOGGER.info("Gripss compare VCFs complete");
    }

    private void loadOriginalVariants(final String vcfFile)
    {
        mVariantBuilder.clearState();

        VcfFileReader reader = new VcfFileReader(vcfFile);

        VCFHeader vcfHeader = reader.vcfHeader();
        GenotypeIds genotypeIds = fromVcfHeader(vcfHeader, mReferenceId, mSampleId);

        if(genotypeIds == null)
        {
            System.exit(1);
        }

        GR_LOGGER.info("genetype info: ref({}: {}) tumor({}: {})",
                genotypeIds.ReferenceOrdinal, genotypeIds.ReferenceId, genotypeIds.TumorOrdinal, genotypeIds.TumorId);

        try
        {
            for(VariantContext variantContext : reader.iterator())
            {
                SvData svData = mVariantBuilder.checkCreateVariant(variantContext, genotypeIds);

                if(svData == null || svData.type() == INF)
                    continue;

                mOriginalSvData.put(svData.id(), svData);

                if(mKeyByCoords)
                {
                    String chrPair = chromosomePair(svData);
                    List<SvData> svList = mOriginalCoordsSvData.get(chrPair);

                    if(svList == null)
                    {
                        svList = Lists.newArrayList();
                        mOriginalCoordsSvData.put(chrPair, svList);
                    }

                    svList.add(svData);
                }
            }
        }
        catch(Exception e)
        {
            GR_LOGGER.error("error reading vcf({}): {}", vcfFile, e.toString());
        }
    }

    private static String chromosomePair(final SvData sv)
    {
        if(sv.isSgl())
            return sv.chromosomeStart();
        else
            return format("%s_%s", sv.chromosomeStart(), sv.chromosomeEnd());
    }

    private void compareVariants(final String newVcfFile)
    {
        mVariantBuilder.clearState();

        GR_LOGGER.info("loading new VCF({})", newVcfFile);

        VcfFileReader reader = new VcfFileReader(newVcfFile);

        VCFHeader vcfHeader = reader.vcfHeader();
        GenotypeIds genotypeIds = fromVcfHeader(vcfHeader, mReferenceId, mSampleId);

        if(genotypeIds == null)
            System.exit(1);

        int diffCount = 0;

        try
        {
            int newSvCount = 0;

            for(VariantContext variantContext : reader.iterator())
            {
                SvData newSv = mVariantBuilder.checkCreateVariant(variantContext, genotypeIds);

                if(newSv == null || newSv.type() == INF)
                    continue;

                ++newSvCount;

                if((newSvCount % 10000 == 0))
                {
                    GR_LOGGER.debug("processed {} variants", newSvCount);
                }

                SvData origSv = findOriginalSv(newSv);

                if(origSv == null)
                {
                    writeDiffs(null, newSv, "NO_ORIG", "", "");
                    ++diffCount;
                    continue;
                }

                mOriginalSvData.remove(origSv.id());

                // types of diffs:
                // filters - where both have filters but different ones
                // filtered vs PASS - likely due to rescue
                // LINKED_BY local and remote -

                if(!mKeyByCoords)
                {
                    if(origSv.posStart() != newSv.posStart() || origSv.posEnd() != newSv.posEnd()
                    || origSv.orientStart() != newSv.orientStart() || origSv.orientEnd() != newSv.orientEnd())
                    {
                        ++diffCount;
                        writeDiffs(origSv, newSv, "COORDS", makeSvCoords(origSv), makeSvCoords(newSv));
                        continue;
                    }
                }

                Breakend origStart = origSv.breakends()[SE_START];
                Breakend newStart = newSv.breakends()[SE_START];

                Set<String> origFilters = origStart.Context.getFilters();
                Set<String> newFilters = newStart.Context.getFilters();

                // first check for a difference in PASS vs not
                Set<String> origFilterDiffs = origFilters.stream().filter(x -> !newFilters.contains(x)).collect(Collectors.toSet());
                Set<String> newFilterDiffs = newFilters.stream().filter(x -> !origFilters.contains(x)).collect(Collectors.toSet());

                boolean ignorePonDiff = mIgnorePonDiff
                        && ((origFilterDiffs.isEmpty() && newFilterDiffs.size() == 1 && newFilterDiffs.contains(PON_FILTER_PON))
                        || (newFilterDiffs.isEmpty() && origFilterDiffs.size() == 1 && origFilterDiffs.contains(PON_FILTER_PON)));

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

                    ++diffCount;
                    continue;
                }

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

            GR_LOGGER.info("loaded {} new SVs", newSvCount);
        }
        catch(Exception e)
        {
            GR_LOGGER.error("error reading vcf({}): {}", newVcfFile, e.toString());
        }

        for(SvData origSv : mOriginalSvData.values())
        {
            ++diffCount;
            writeDiffs(origSv, null, "NO_NEW", "", "");
        }

        GR_LOGGER.info("diffTotal({})", diffCount);
    }

    private boolean checkVcfFieldDiff(final SvData origSv, final SvData newSv, final VcfCompareField compareField)
    {
        boolean hasDiff = false;

        if(compareField.TypeScope == VariantTypeScope.SV && origSv.isSgl())
            return false;
        if(compareField.TypeScope == VariantTypeScope.SGL && !origSv.isSgl())
            return false;

        if(compareField.Scope == GenotypeScope.COMBINED)
        {
            double origValue = origSv.contextStart().getAttributeAsDouble(compareField.VcfTag, 0);
            double newValue = newSv.contextStart().getAttributeAsDouble(compareField.VcfTag, 0);

            if(hasDiff(origValue, newValue, compareField.DiffAbs, compareField.DiffPerc))
            {
                hasDiff = true;
                writeDiffs(origSv, newSv, compareField.VcfTag, String.valueOf(origValue), String.valueOf(newValue));
            }
        }
        else
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                if(se == SE_END && origSv.isSgl())
                    continue;

                Breakend origBreakend = origSv.breakends()[se];
                Breakend newBreakend = newSv.breakends()[se];

                if(!mReferenceId.isEmpty() && (compareField.Scope == GenotypeScope.NORMAL || compareField.Scope == GenotypeScope.BOTH))
                {
                    double origValue = getGenotypeAttributeAsDouble(origBreakend.RefGenotype, compareField.VcfTag, 0);
                    double newValue = getGenotypeAttributeAsDouble(newBreakend.RefGenotype, compareField.VcfTag, 0);

                    if(hasDiff(origValue, newValue, compareField.DiffAbs, compareField.DiffPerc))
                    {
                        hasDiff = true;
                        writeDiffs(origSv, newSv, format("REF_%s", compareField.VcfTag), String.valueOf(origValue), String.valueOf(newValue));
                    }
                }

                if(compareField.Scope == GenotypeScope.TUMOR || compareField.Scope == GenotypeScope.BOTH)
                {
                    double origValue = getGenotypeAttributeAsDouble(origBreakend.TumorGenotype, compareField.VcfTag, 0);
                    double newValue = getGenotypeAttributeAsDouble(newBreakend.TumorGenotype, compareField.VcfTag, 0);

                    if(hasDiff(origValue, newValue, compareField.DiffAbs, compareField.DiffPerc))
                    {
                        hasDiff = true;
                        writeDiffs(origSv, newSv, format("TUMOR_%s", compareField.VcfTag), String.valueOf(origValue), String.valueOf(newValue));
                    }
                }
            }
        }

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

    private SvData findOriginalSv(final SvData newSv)
    {
        if(!mKeyByCoords)
            return mOriginalSvData.get(newSv.id());

        List<SvData> svList = mOriginalCoordsSvData.get(chromosomePair(newSv));

        if(svList == null)
            return null;

        for(int i = 0; i < svList.size(); ++i)
        {
            SvData sv = svList.get(i);

            if(sv.type() != newSv.type())
                continue;

            if(sv.posStart() != newSv.posStart() || sv.orientStart() != newSv.orientStart())
                continue;

            if(!sv.isSgl() && (sv.posEnd() != newSv.posEnd() || sv.orientEnd() != newSv.orientEnd()))
                continue;

            svList.remove(i);
            return sv;
        }

        return null;
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

            GR_LOGGER.info("writing comparison file: {}", fileName);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            if(mKeyByCoords)
                writer.write("OrigId\tNewId");
            else
                writer.write("SvId");

            writer.write("\tCoords\tType\tDiffType\tOrigValue\tNewValue\tOrigQual\tNewQual\tOrigFilters\tNewFilters");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            GR_LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    private void writeDiffs(final SvData origSv, final SvData newSv, final String diffType, final String origValue, final String newValue)
    {
        try
        {
            String coords = origSv != null ? makeSvCoords(origSv) : makeSvCoords(newSv);
            StructuralVariantType type = origSv != null ? origSv.type() : newSv.type();

            if(mKeyByCoords)
            {
                mWriter.write(format("%s\t%s",
                        origSv != null ? origSv.id() : "", newSv != null ? newSv.id() : ""));
            }
            else
            {
                mWriter.write(format("%s", origSv != null ? origSv.id() : newSv.id()));
            }

            mWriter.write(format("\t%s\t%s\t%s\t%s\t%s",
                    coords, type, diffType, origValue, newValue));

            mWriter.write(format("\t%.1f\t%.1f\t%s\t%s",
                    origSv != null ? origSv.breakendStart().Qual : -1, newSv != null ? newSv.breakendStart().Qual : -1,
                    origSv != null ? filtersStr(origSv.breakendStart().Context.getFilters(), true) : "",
                    newSv != null ? filtersStr(newSv.breakendStart().Context.getFilters(), true) : ""));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            GR_LOGGER.error("failed to write output file: {}", e.toString());
        }
    }

    private static String makeSvCoords(final SvData sv)
    {
        if(sv.isSgl())
        {
            return format("%s:%d:%d", sv.chromosomeStart(), sv.posStart(), sv.orientStart());
        }
        else
        {
            return format("%s:%d:%d-%s:%d:%d",
                    sv.chromosomeStart(), sv.posStart(), sv.orientStart(), sv.chromosomeEnd(), sv.posEnd(), sv.orientEnd());
        }
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
        SV,
        SGL,
        BOTH;
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
        configBuilder.addFlag(KEY_BY_COORDS, "Match SVs on coords rather than VcfId");
        configBuilder.addFlag(WRITE_ALL_DIFFS, "Write all VCF field diffs, not just the first");
        configBuilder.addFlag(GRIDSS_ONLY, "Only compare fields written by Grids (ie no Gripss)");
        configBuilder.addFlag(REF_DEPTH_ONLY, "Only compare reference depth fields");

        HotspotCache.addConfig(configBuilder);
        TargetRegions.addConfig(configBuilder);
        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);
        setLogLevel(configBuilder);

        GripssCompareVcfs gripssCompare = new GripssCompareVcfs(configBuilder);
        gripssCompare.run();
    }
}

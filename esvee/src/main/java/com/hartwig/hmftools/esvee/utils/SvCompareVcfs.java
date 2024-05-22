package com.hartwig.hmftools.esvee.utils;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.CHR_PREFIX;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.DISC_FRAGS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.IHOMPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SPLIT_FRAGS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.STRAND_BIAS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.DISCORDANT_READS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.SGL_FRAG_COUNT;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.SPLIT_READS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.SV_FRAG_COUNT;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.QUAL;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsDouble;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.common.variant.GenotypeIds.fromVcfHeader;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.CommonUtils.formSvType;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.sv.VariantAltInsertCoords;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.GenotypeIds;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class SvCompareVcfs
{
    private final String mSampleId;
    private final String mReferenceId;
    private final String mOriginalVcf;
    private final String mOriginalUnfilteredVcf;
    private final String mNewVcf;
    private final String mNewUnfilteredVcf;
    private final String mOutputDir;
    private final String mOutputId;

    private final Map<String,List<VariantBreakend>> mOrigChrBreakendMap;
    private final Map<String,List<VariantBreakend>> mNewChrBreakendMap;

    private final BufferedWriter mWriter;

    private final boolean mCompareFilters; // only relevant for the same caller (ie Esvee or Gridss)
    private final boolean mOriginalIsGridss;
    private final boolean mIgnorePonDiff;

    private final List<VcfCompareField> mVcfCheckFields;
    private RefGenomeVersion mRefGenomeVersion;

    private int mMatchedCount;
    private int mDiffCount;

    private static final String ORIGINAL_VCF = "original_vcf";
    private static final String NEW_VCF = "new_vcf";

    private static final String ORIGINAL_UNFILTERED_VCF = "original_unfiltered_vcf";
    private static final String NEW_UNFILTERED_VCF = "new_unfiltered_vcf";

    private static final String COMPARE_FILTERS = "compare_filters";
    private static final String FIELD_FILTERS = "FILTERS";
    private static final String IGNORE_PON_DIFF = "ignore_pon_diff";

    private static final int DEFAULT_MAX_DIFF = 20;
    private static final double DEFAULT_MAX_DIFF_PERC = 0.2;

    public SvCompareVcfs(final ConfigBuilder configBuilder)
    {
        mSampleId = configBuilder.getValue(SAMPLE);
        mReferenceId = configBuilder.getValue(REFERENCE, "");
        mOriginalVcf = configBuilder.getValue(ORIGINAL_VCF);
        mNewVcf = configBuilder.getValue(NEW_VCF);
        mOriginalUnfilteredVcf = configBuilder.getValue(ORIGINAL_UNFILTERED_VCF);
        mNewUnfilteredVcf = configBuilder.getValue(NEW_UNFILTERED_VCF);
        mOutputDir = parseOutputDir(configBuilder);
        mOutputId = configBuilder.getValue(OUTPUT_ID);

        mOrigChrBreakendMap = Maps.newHashMap();
        mNewChrBreakendMap = Maps.newHashMap();
        mRefGenomeVersion = null;

        mOriginalIsGridss = mOriginalVcf.contains("gridss") || mOriginalVcf.contains("gripss");
        mIgnorePonDiff = configBuilder.hasFlag(IGNORE_PON_DIFF);
        mCompareFilters = configBuilder.hasFlag(COMPARE_FILTERS);

        mVcfCheckFields = Lists.newArrayList();
        addComparisonFields();

        mMatchedCount = 0;
        mDiffCount = 0;

        mWriter = initialiseWriter();
    }

    private void addComparisonFields()
    {
        mVcfCheckFields.add(new VcfCompareField(
                REF_DEPTH, FieldType.INTEGER, GenotypeScope.BOTH, VariantTypeScope.BREAKEND, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));

        mVcfCheckFields.add(new VcfCompareField(
                REF_DEPTH_PAIR, FieldType.INTEGER, GenotypeScope.BOTH, VariantTypeScope.BREAKEND, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));

        mVcfCheckFields.add(new VcfCompareField(
                STRAND_BIAS, FieldType.DECIMAL, GenotypeScope.COMBINED, VariantTypeScope.BREAKEND, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));

        mVcfCheckFields.add(new VcfCompareField(
                IHOMPOS, FieldType.STRING, GenotypeScope.COMBINED, VariantTypeScope.BREAKEND, 0, 0));

        if(!mOriginalIsGridss)
        {
            mVcfCheckFields.add(new VcfCompareField(
                    QUAL, FieldType.INTEGER, GenotypeScope.BOTH, VariantTypeScope.BREAKEND, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));

            mVcfCheckFields.add(new VcfCompareField(
                    SPLIT_FRAGS, FieldType.INTEGER, GenotypeScope.BOTH, VariantTypeScope.BREAKEND, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));

            mVcfCheckFields.add(new VcfCompareField(
                    DISC_FRAGS, FieldType.INTEGER, GenotypeScope.BOTH, VariantTypeScope.BREAKEND, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));
        }

        /*
        mVcfCheckFields.add(new VcfCompareField(INDEL_COUNT, GenotypeScope.BOTH, VariantTypeScope.BOTH, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));

        mVcfCheckFields.add(new VcfCompareField(READ_PAIRS, GenotypeScope.BOTH, VariantTypeScope.SV, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));
        mVcfCheckFields.add(new VcfCompareField(GRIDSS_ASRP, GenotypeScope.BOTH, VariantTypeScope.SV, DEFAULT_MAX_DIFF, DEFAULT_MAX_DIFF_PERC));

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
    }

    public void run()
    {
        if(mOriginalVcf == null || mNewVcf == null)
        {
            SV_LOGGER.error("missing VCFs");
            return;
        }

        SV_LOGGER.info("loading original VCF({})", mOriginalVcf);
        loadVariants(mOriginalVcf, mOrigChrBreakendMap);
        loadVariants(mNewVcf, mNewChrBreakendMap);

        compareVariants();

        checkFilteredVariants();

        writeUnmatchedVariants();

        closeBufferedWriter(mWriter);

        SV_LOGGER.info("Esvee compare VCFs complete");
    }

    private void loadVariants(final String vcfFile, final Map<String,List<VariantBreakend>> chrBreakendMap)
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
        List<VariantBreakend> breakends = null;

        for(VariantContext variantContext : reader.iterator())
        {
            String chromosome = variantContext.getContig();

            if(mRefGenomeVersion == null)
                mRefGenomeVersion = chromosome.startsWith(CHR_PREFIX) ? V38 : V37;

            if(!currentChr.equals(chromosome))
            {
                currentChr = chromosome;
                breakends = Lists.newArrayList();
                chrBreakendMap.put(chromosome, breakends);
            }

            breakends.add(new VariantBreakend(variantContext));
        }

        SV_LOGGER.info("loaded {} SVs from {})",
                chrBreakendMap.values().stream().mapToInt(x -> x.size()).sum(), vcfFile);
    }

    private void compareVariants()
    {
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mRefGenomeVersion.versionedChromosome(chromosome.toString());

            List<VariantBreakend> origBreakends = mOrigChrBreakendMap.get(chrStr);
            List<VariantBreakend> newBreakends = mNewChrBreakendMap.get(chrStr);

            if(origBreakends == null || newBreakends == null)
                continue;

            int origIndex = 0;
            int newIndex = 0;

            int processedCount = 0;

            VariantBreakend origBreakend = !origBreakends.isEmpty() ? origBreakends.get(origIndex) : null;
            VariantBreakend newBreakend = !newBreakends.isEmpty() ? newBreakends.get(newIndex) : null;

            while(origBreakend != null && newBreakend != null)
            {
                if(processedCount > 0 && (processedCount % 1000) == 0)
                {
                    SV_LOGGER.debug("processed {} variants", processedCount);
                }

                // scenarios:
                // variants match exactly or match within homology (no distinction for now)
                // original is ahead, so need to skip past new (or multiple)
                // new is ahead, so need to skip past original

                if(origBreakend.matchesWithinHomology(newBreakend))
                {
                    ++mMatchedCount;
                    compareBreakends(origBreakend, newBreakend);

                    origBreakends.remove(origIndex);
                    newBreakends.remove(newIndex);

                    origBreakend = origIndex < origBreakends.size() ? origBreakends.get(origIndex) : null;
                    newBreakend = newIndex < newBreakends.size() ? newBreakends.get(newIndex) : null;
                    continue;
                }

                if(origBreakend.isLowerPosition(newBreakend))
                {
                    // proceed to the next match or past the other breakend
                    while(origBreakend != null && origBreakend.maxPosition() < newBreakend.minPosition())
                    {
                        ++origIndex;

                        if(origIndex >= origBreakends.size())
                        {
                            origBreakend = null;
                            break;
                        }

                        origBreakend = origBreakends.get(origIndex);
                    }
                }
                else
                {
                    while(newBreakend != null && newBreakend.maxPosition() < origBreakend.minPosition())
                    {
                        ++newIndex;

                        if(newIndex >= newBreakends.size())
                        {
                            newBreakend = null;
                            break;
                        }

                        newBreakend = newBreakends.get(newIndex);
                    }
                }
            }

            // SV_LOGGER.info("loaded {} new SVs", newSvCount);
        }

        // SV_LOGGER.info("diffTotal({})", diffCount);
    }

    private static final String DIFF_NO_ORIG = "NO_ORIG";
    private static final String DIFF_NO_NEW = "NO_NEW";

    private void writeUnmatchedVariants()
    {
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mRefGenomeVersion.versionedChromosome(chromosome.toString());

            for(int i = 0; i <= 1; ++i)
            {
                boolean isOriginal = (i == 0);
                List<VariantBreakend> breakends = isOriginal ? mOrigChrBreakendMap.get(chrStr) : mNewChrBreakendMap.get(chrStr);

                if(breakends == null)
                    continue;

                for(VariantBreakend breakend : breakends)
                {
                    String breakendInfo = format("QUAL(%.0f)", breakend.Context.getPhredScaledQual());

                    if(isOriginal)
                        writeDiffs(breakend, null, DIFF_NO_NEW, breakendInfo, "");
                    else
                        writeDiffs(null, breakend, DIFF_NO_ORIG, "", breakendInfo);
                }
            }
        }
    }

    private void checkFilteredVariants()
    {
        if(mNewUnfilteredVcf == null && mOriginalUnfilteredVcf == null)
            return;

        for(int i = 0; i <= 0; ++i)
        {
            boolean checkOriginal = (i == 0);
            Map<String,List<VariantBreakend>> chrBreakendMap = checkOriginal ? mOrigChrBreakendMap : mNewChrBreakendMap;

            int unmatchedBreakends = chrBreakendMap.values().stream().mapToInt(x -> x.size()).sum();

            if(unmatchedBreakends == 0)
                continue;

            String unfilteredVcf = checkOriginal ? mNewUnfilteredVcf : mOriginalUnfilteredVcf;

            if(unfilteredVcf == null)
                continue;

            SV_LOGGER.debug("checking unmatched {} {} breakends in {} unfiltered VCF",
                    unmatchedBreakends, checkOriginal ? "original" : "new", checkOriginal ? "new" : "original");

            VcfFileReader reader = new VcfFileReader(unfilteredVcf);

            String currentChr = "";
            List<VariantBreakend> breakends = null;

            for(VariantContext variantContext : reader.iterator())
            {
                String chromosome = variantContext.getContig();

                if(mRefGenomeVersion == null)
                    mRefGenomeVersion = chromosome.startsWith(CHR_PREFIX) ? V38 : V37;

                if(!currentChr.equals(chromosome))
                {
                    currentChr = chromosome;
                    breakends = chrBreakendMap.get(chromosome);
                }

                if(breakends == null || breakends.isEmpty())
                    continue;

                VariantBreakend filteredBreakend = new VariantBreakend(variantContext);

                for(int j = 0; j < breakends.size(); ++j)
                {
                    VariantBreakend breakend = breakends.get(j);

                    if(breakend.matchesWithinHomology(filteredBreakend))
                    {
                        String filters = filtersStr(filteredBreakend.Context.getFilters(), true);

                        if(checkOriginal)
                            writeDiffs(breakend, filteredBreakend, FIELD_FILTERS, PASS, filters);
                        else
                            writeDiffs(filteredBreakend, breakend, FIELD_FILTERS, filters, PASS);

                        breakends.remove(j);
                        --unmatchedBreakends;
                        break;
                    }
                    else if(filteredBreakend.isLowerPosition(breakend))
                    {
                        break;
                    }
                }

                if(unmatchedBreakends == 0 || breakends.isEmpty())
                    break;
            }

            if(unmatchedBreakends == 0)
                break;
        }
    }

    private void compareBreakends(final VariantBreakend origBreakend, final VariantBreakend newBreakend)
    {
        // compare each field in turn, building up a set of diffs to write to file
        if(mCompareFilters)
        {
            Set<String> origFilters = origBreakend.Context.getFilters();
            Set<String> newFilters = newBreakend.Context.getFilters();

            // first check for a difference in PASS vs not
            Set<String> origFilterDiffs = origFilters.stream().filter(x -> !newFilters.contains(x)).collect(Collectors.toSet());
            Set<String> newFilterDiffs = newFilters.stream().filter(x -> !origFilters.contains(x)).collect(Collectors.toSet());

            /*
            boolean ignorePonDiff = mIgnorePonDiff
                    && ((origFilterDiffs.isEmpty() && newFilterDiffs.size() == 1 && newFilterDiffs.contains(PON_FILTER_PON))
                        || (newFilterDiffs.isEmpty() && origFilterDiffs.size() == 1 && origFilterDiffs.contains(PON_FILTER_PON)));

            */
            if(!newFilterDiffs.isEmpty() || !origFilterDiffs.isEmpty())
            {
                writeDiffs(
                        origBreakend, newBreakend, FIELD_FILTERS,
                        filtersStr(origFilterDiffs, true), filtersStr(newFilterDiffs, true));
            }

        }

        // compare breakend attributes

        // insert sequence
        if(!origBreakend.Coords.InsertSequence.equals(newBreakend.Coords.InsertSequence))
        {
            writeDiffs(origBreakend, newBreakend, "INSSEQ", origBreakend.Coords.InsertSequence, newBreakend.Coords.InsertSequence);
        }

        // check specified VCF tags
        for(VcfCompareField compareField : mVcfCheckFields)
        {
            checkVcfFieldDiff(origBreakend, newBreakend, compareField);
        }
    }

    private void checkVcfFieldDiff(final VariantBreakend origBreakend, final VariantBreakend newBreakend, final VcfCompareField compareField)
    {
        if(compareField.TypeScope == VariantTypeScope.VARIANT && origBreakend.isEnd())
        {
            // skip evaluating on the end breakend
            return;
        }

        String origVcfTag = mOriginalIsGridss ? mapToGridssVcfTag(compareField.VcfTag, origBreakend.isSingle()) : compareField.VcfTag;

        if(compareField.Scope == GenotypeScope.COMBINED)
        {
            if(compareField.Type == FieldType.STRING)
            {
                checkField(
                        origBreakend, newBreakend, compareField,
                        origBreakend.Context.getAttributeAsString(origVcfTag, ""),
                        newBreakend.Context.getAttributeAsString(compareField.VcfTag, ""));
            }
            else if(compareField.Type == FieldType.INTEGER)
            {
                checkField(
                        origBreakend, newBreakend, compareField,
                        origBreakend.Context.getAttributeAsInt(origVcfTag, 0),
                        newBreakend.Context.getAttributeAsInt(compareField.VcfTag, 0));
            }
            else
            {
                checkField(
                        origBreakend, newBreakend, compareField,
                        origBreakend.Context.getAttributeAsDouble(origVcfTag, 0),
                        newBreakend.Context.getAttributeAsDouble(compareField.VcfTag, 0));
            }
        }
        else
        {
            for(Genotype origGenotype : origBreakend.Context.getGenotypes())
            {
                if(origGenotype.getSampleName().equals(mReferenceId) && compareField.Scope == GenotypeScope.TUMOR)
                    continue;
                else if(origGenotype.getSampleName().equals(mSampleId) && compareField.Scope == GenotypeScope.NORMAL)
                    continue;

                Genotype newGenotype = newBreakend.Context.getGenotype(origGenotype.getSampleName());

                double origValue = getGenotypeAttributeAsDouble(origGenotype, origVcfTag, 0);
                double newValue = getGenotypeAttributeAsDouble(newGenotype, compareField.VcfTag, 0);

                if(compareField.Type == FieldType.INTEGER)
                {
                    checkField(
                            origBreakend, newBreakend, compareField,
                            getGenotypeAttributeAsInt(origGenotype, origVcfTag, 0),
                            getGenotypeAttributeAsInt(newGenotype, compareField.VcfTag, 0));
                }
                else
                {
                    checkField(
                            origBreakend, newBreakend, compareField,
                            getGenotypeAttributeAsDouble(origGenotype, origVcfTag, 0),
                            getGenotypeAttributeAsDouble(newGenotype, compareField.VcfTag, 0));
                }

                if(hasDiff(origValue, newValue, compareField.DiffAbs, compareField.DiffPerc))
                {
                    writeDiffs(origBreakend, newBreakend, compareField.VcfTag, String.valueOf(origValue), String.valueOf(newValue));
                }
            }
        }
    }

    private String mapToGridssVcfTag(final String vcfTag, boolean isSingle)
    {
        if(vcfTag.equals(SPLIT_FRAGS))
            return SPLIT_READS;

        if(vcfTag.equals(DISC_FRAGS))
            return DISCORDANT_READS;

        if(vcfTag.equals(TOTAL_FRAGS))
            return isSingle ? SV_FRAG_COUNT : SGL_FRAG_COUNT;

        return vcfTag;
    }

    private void checkField(
            final VariantBreakend origBreakend, final VariantBreakend newBreakend, final VcfCompareField compareField,
            final String origValue, final String newValue)
    {
        if(!origValue.equals(newValue))
        {
            writeDiffs(origBreakend, newBreakend, compareField.VcfTag, origValue, newValue);
        }
    }

    private void checkField(
            final VariantBreakend origBreakend, final VariantBreakend newBreakend, final VcfCompareField compareField,
            int origValue, int newValue)
    {
        if(hasDiff(origValue, newValue, compareField.DiffAbs, compareField.DiffPerc))
        {
            writeDiffs(origBreakend, newBreakend, compareField.VcfTag, String.valueOf(origValue), String.valueOf(newValue));
        }
    }

    private void checkField(
            final VariantBreakend origBreakend, final VariantBreakend newBreakend, final VcfCompareField compareField,
            double origValue, double newValue)
    {
        if(hasDiff(origValue, newValue, compareField.DiffAbs, compareField.DiffPerc))
        {
            writeDiffs(origBreakend, newBreakend, compareField.VcfTag, format("%.2f", origValue), format("%.2f", newValue));
        }
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

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("OrigId").add("NewId").add("Coords").add("SvType").add("DiffType").add("OrigValue").add("NewValue");
            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    private void writeDiffs(
            final VariantBreakend origBreakend, final VariantBreakend newBreakend, final String diffType,
            final String origValue, final String newValue)
    {
        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(origBreakend != null ? origBreakend.Context.getID() : "-1");
            sj.add(newBreakend != null ? newBreakend.Context.getID() : "-1");

            String breakendCoords = origBreakend != null ? origBreakend.coordStr() : newBreakend.coordStr();
            sj.add(breakendCoords);

            StructuralVariantType svType = origBreakend != null ? origBreakend.svType() : newBreakend.svType();
            sj.add(svType.toString());

            sj.add(diffType);
            sj.add(origValue);
            sj.add(newValue);
            mWriter.write(sj.toString());
            mWriter.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write output file: {}", e.toString());
        }
    }

    private class VariantBreakend
    {
        public final VariantContext Context;
        public final VariantAltInsertCoords Coords;
        public final int Position;
        public final int[] Cipos;

        public VariantBreakend(final VariantContext context)
        {
            Context = context;
            Position = Context.getStart();

            String alt = context.getAlternateAllele(0).getDisplayString();
            Coords = VariantAltInsertCoords.fromRefAlt(alt, alt.substring(0, 1));

            List<Integer> ciposList = context.getAttributeAsIntList(CIPOS, 0);
            Cipos = ciposList.size() == 2 ? new int[] { ciposList.get(0), ciposList.get(1) } : new int[] {0, 0};
        }

        public int minPosition() { return Position + Cipos[0];}
        public int maxPosition() { return Position + Cipos[1];}

        public boolean isSingle() { return Coords.OtherChromsome.isEmpty(); }

        public boolean isEnd()
        {
            if(Coords.OtherChromsome.isEmpty())
                return false;

            if(Context.getContig().equals(Coords.OtherChromsome))
            {
                return Position > Coords.OtherPosition;
            }
            else
            {
                return HumanChromosome.lowerChromosome(Coords.OtherChromsome, Context.getContig());
            }
        }

        public boolean exactMatch(final VariantBreakend other)
        {
            return other.Position == Position && other.Coords.Orient == Coords.Orient;
        }

        public boolean matchesWithinHomology(final VariantBreakend other)
        {
            if(other.Coords.Orient != Coords.Orient)
                return false;

            return positionsOverlap(minPosition(), maxPosition(),  other.minPosition(), other.maxPosition());
        }

        public boolean isLowerPosition(final VariantBreakend other)
        {
            if(Position == other.Position)
            {
                if(Coords.Orient == other.Coords.Orient)
                    return false;
                else
                    return Coords.Orient == Orientation.FORWARD;
            }
            else
            {
                return Position < other.Position;
            }
        }

        public String toString() { return format("%d:%d cipos(%d,%d)", Position, Coords.Orient.asByte(), Cipos[0], Cipos[1]); }

        public String coordStr() { return format("%s:%d:%d", Context.getContig(), Position, Coords.Orient.asByte()); }

        public StructuralVariantType svType()
        {
            if(Coords.OtherChromsome.equals(""))
                return SGL;

            return formSvType(
                    Context.getContig(), Coords.OtherChromsome, Position, Coords.OtherPosition,
                    Coords.Orient, Coords.OtherOrient, !Coords.InsertSequence.isEmpty());
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
        BREAKEND,
        VARIANT;
    }

    private enum FieldType
    {
        STRING,
        INTEGER,
        DECIMAL;
    }

    private class VcfCompareField
    {
        public final String VcfTag;
        public final FieldType Type;
        public final GenotypeScope Scope;
        public final double DiffAbs;
        public final double DiffPerc;
        public final VariantTypeScope TypeScope;

        public VcfCompareField(
                final String vcfTag, final FieldType type, final GenotypeScope scope, final VariantTypeScope typeScope,
                final double diffAbs, final double diffPerc)
        {
            VcfTag = vcfTag;
            Type = type;
            Scope = scope;
            TypeScope = typeScope;
            DiffAbs = diffAbs;
            DiffPerc = diffPerc;
        }

        public VcfCompareField(
                final String vcfTag, final FieldType type, final GenotypeScope scope, final VariantTypeScope typeScope)
        {
            this(vcfTag, type, scope, typeScope, 0, 0);
        }

        public String toString() { return format("tag(%s) scope(%s) st(%s)", VcfTag, Scope, TypeScope); }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();

        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addConfigItem(REFERENCE, false, REFERENCE_DESC);
        configBuilder.addPath(ORIGINAL_VCF, true, "Path to the old VCF file");
        configBuilder.addPath(ORIGINAL_UNFILTERED_VCF, false, "Path to the old unfiltered VCF file");
        configBuilder.addPath(NEW_VCF, true, "Path to the new VCF file");
        configBuilder.addPath(NEW_UNFILTERED_VCF, false, "Path to the new unfiltered VCF file");
        configBuilder.addFlag(IGNORE_PON_DIFF, "Ignore diffs if just PON filter");
        configBuilder.addFlag(COMPARE_FILTERS, "Compare filters within same SV caller");

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);
        setLogLevel(configBuilder);

        SvCompareVcfs svVcfCompare = new SvCompareVcfs(configBuilder);
        svVcfCompare.run();
    }
}
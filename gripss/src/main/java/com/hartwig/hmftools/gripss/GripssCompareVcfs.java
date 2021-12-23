package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.LOCAL_LINKED_BY;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PON_FILTER_PON;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REMOTE_LINKED_BY;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.gripss.GripssConfig.GR_LOGGER;
import static com.hartwig.hmftools.gripss.GripssConfig.SAMPLE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.GenotypeIds;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.common.VcfUtils;
import com.hartwig.hmftools.gripss.filters.HotspotCache;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

public class GripssCompareVcfs
{
    private final String mSampleId;
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

    private static final String ORIGINAL_VCF = "original_vcf";
    private static final String NEW_VCF = "new_vcf";
    private static final String IGNORE_PON_DIFF = "ignore_pon_diff";
    private static final String KEY_BY_COORDS = "key_by_coords";

    public GripssCompareVcfs(final CommandLine cmd)
    {
        mSampleId = cmd.getOptionValue(SAMPLE);
        mOriginalVcf = cmd.getOptionValue(ORIGINAL_VCF);
        mNewVcf = cmd.getOptionValue(NEW_VCF);
        mOutputDir = parseOutputDir(cmd);
        mOutputId = cmd.getOptionValue(OUTPUT_ID);

        mOriginalSvData = Maps.newHashMap();
        mOriginalCoordsSvData = Maps.newHashMap();

        mVariantBuilder = new VariantBuilder(null, new HotspotCache(cmd));

        mIgnorePonDiff = cmd.hasOption(IGNORE_PON_DIFF);
        mKeyByCoords = cmd.hasOption(KEY_BY_COORDS);

        mWriter = initialiseWriter();
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
    }

    private void loadOriginalVariants(final String vcfFile)
    {
        mVariantBuilder.clearState();

        final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                vcfFile, new VCFCodec(), false);

        VCFHeader vcfHeader = (VCFHeader)reader.getHeader();
        GenotypeIds genotypeIds = VcfUtils.parseVcfSampleIds(vcfHeader, "", mSampleId);

        if(genotypeIds == null)
        {
            System.exit(1);
        }

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
        catch(IOException e)
        {
            GR_LOGGER.error("error reading vcf({}): {}", vcfFile, e.toString());
        }
    }

    private static String chromosomePair(final SvData sv)
    {
        if(sv.isSgl())
            return sv.chromosomeStart();
        else
            return String.format("%s_%s", sv.chromosomeStart(), sv.chromosomeEnd());
    }

    private void compareVariants(final String newVcfFile)
    {
        mVariantBuilder.clearState();

        GR_LOGGER.info("loading new VCF({})", newVcfFile);

        final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                newVcfFile, new VCFCodec(), false);

        VCFHeader vcfHeader = (VCFHeader)reader.getHeader();
        GenotypeIds genotypeIds = VcfUtils.parseVcfSampleIds(vcfHeader, "", mSampleId);

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

                if(origSv.posStart() != newSv.posStart() || origSv.posEnd() != newSv.posEnd()
                || origSv.orientStart() != newSv.orientStart() || origSv.orientEnd() != newSv.orientEnd())
                {
                    ++diffCount;
                    writeDiffs(origSv, newSv, "COORDS", makeSvCoords(origSv), makeSvCoords(newSv));
                    continue;
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
                        writeDiffs(origSv, newSv, "FILTER_PASS", filtersStr(origFilterDiffs), filtersStr(newFilterDiffs));
                    }
                    else
                    {
                        writeDiffs(origSv, newSv, "FILTER_DIFF", filtersStr(origFilterDiffs), filtersStr(newFilterDiffs));
                    }

                    ++diffCount;
                    continue;
                }

                // check local and remote linked by for assembled links
                boolean origHasStartAssembled = origStart.Context.getAttributeAsString(LOCAL_LINKED_BY, "").contains("asm");
                boolean newHasStartAssembled = newStart.Context.getAttributeAsString(LOCAL_LINKED_BY, "").contains("asm");
                boolean origHasEndAssembled = !origSv.isSgl() && origStart.Context.getAttributeAsString(REMOTE_LINKED_BY, "").contains("asm");
                boolean newHasEndAssembled = !origSv.isSgl() && newStart.Context.getAttributeAsString(REMOTE_LINKED_BY, "").contains("asm");

                if(origHasStartAssembled != newHasStartAssembled || origHasEndAssembled != newHasEndAssembled)
                {
                    writeDiffs(
                            origSv, newSv, "ASSEMBLY",
                            String.format("%s_%s", origHasStartAssembled, origHasEndAssembled),
                            String.format("%s_%s", newHasStartAssembled, newHasEndAssembled));
                    ++diffCount;
                    continue;
                }
            }

            GR_LOGGER.info("loaded {} new SVs", newSvCount);
        }
        catch(IOException e)
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

    private String filtersStr(final Set<String> filters)
    {
        StringJoiner sj = new StringJoiner(";");
        filters.forEach(x -> sj.add(x));
        return sj.toString();
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            String fileName = mOutputDir + mSampleId + ".compare";

            if(mOutputId != null)
                fileName += "." + mOutputId;

            fileName += ".csv";

            GR_LOGGER.info("writing comparison file: {}", fileName);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            if(mKeyByCoords)
                writer.write("OrigId,NewId");
            else
                writer.write("SvId");

            writer.write(",Coords,Type,DiffType,OrigValue,NewValue");
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
                mWriter.write(String.format("%s,%s",
                        origSv != null ? origSv.id() : "", newSv != null ? newSv.id() : ""));
            }
            else
            {
                mWriter.write(String.format("%s", origSv != null ? origSv.id() : newSv.id()));
            }

            mWriter.write(String.format(",%s,%s,%s,%s,%s",
                    coords, type, diffType, origValue, newValue));

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
            return String.format("%s:%d:%d", sv.chromosomeStart(), sv.posStart(), sv.orientStart());
        }
        else
        {
            return String.format("%s:%d:%d-%s:%d:%d",
                    sv.chromosomeStart(), sv.posStart(), sv.orientStart(), sv.chromosomeEnd(), sv.posEnd(), sv.orientEnd());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Name of the tumor sample");
        options.addOption(ORIGINAL_VCF, true, "Optional, name of the reference sample");
        options.addOption(NEW_VCF, true, "Path to the GRIDSS structural variant VCF file");
        options.addOption(IGNORE_PON_DIFF, false, "Ignore diffs if just PON filter");
        options.addOption(KEY_BY_COORDS, false, "Match SVs on coords rather than VcfId");

        addOutputOptions(options);
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        GripssCompareVcfs gripssCompare = new GripssCompareVcfs(cmd);
        gripssCompare.run();

        GR_LOGGER.info("Gripss compare VCFs complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}

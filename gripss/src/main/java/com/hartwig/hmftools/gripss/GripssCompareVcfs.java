package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.gripss.GripssConfig.GR_LOGGER;
import static com.hartwig.hmftools.gripss.GripssConfig.SAMPLE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.GenotypeIds;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.filters.FilterConstants;
import com.hartwig.hmftools.gripss.filters.HotspotCache;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class GripssCompareVcfs
{
    private final String mSampleId;
    private final String mOriginalVcf;
    private final String mNewVcf;
    private final String mOutputDir;
    private final String mOutputId;

    private final VariantBuilder mVariantBuilder;

    private final Map<String,SvData> mOriginalSvData;
    private final Map<String,SvData> mNewSvData;

    private final BufferedWriter mWriter;

    private static final String ORIGINAL_VCF = "original_vcf";
    private static final String NEW_VCF = "new_vcf";

    public GripssCompareVcfs(final CommandLine cmd)
    {
        mSampleId = cmd.getOptionValue(SAMPLE);
        mOriginalVcf = cmd.getOptionValue(ORIGINAL_VCF);
        mNewVcf = cmd.getOptionValue(NEW_VCF);
        mOutputDir = parseOutputDir(cmd);
        mOutputId = cmd.getOptionValue(OUTPUT_ID);

        mOriginalSvData = Maps.newHashMap();
        mNewSvData = Maps.newHashMap();

        mVariantBuilder = new VariantBuilder(FilterConstants.from(cmd), new HotspotCache(cmd));

        mWriter = initialiseWriter();
    }

    public void run()
    {
        if(mOriginalVcf == null || mNewVcf == null)
        {
            GR_LOGGER.error("missing VCFs");
            return;
        }

        GR_LOGGER.info("loading original VCF({})", mOriginalSvData);
        loadVariants(mOriginalVcf, mOriginalSvData);

        GR_LOGGER.info("loading new VCF({})", mNewVcf);
        loadVariants(mNewVcf, mNewSvData);

        GR_LOGGER.info("loaded SVs: original({}) new({})", mOriginalSvData.size(), mNewSvData.size());

        compareVariants();;
        closeBufferedWriter(mWriter);
    }

    private void loadVariants(final String vcfFile, final Map<String,SvData> svDataMap)
    {
        mVariantBuilder.clearState();

        final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                vcfFile, new VCFCodec(), false);

        GenotypeIds genotypeIds = new GenotypeIds(0, 1, mSampleId, mSampleId);

        try
        {
            for(VariantContext variantContext : reader.iterator())
            {
                SvData svData = mVariantBuilder.checkCreateVariant(variantContext, genotypeIds);

                if(svData == null)
                    continue;

                svDataMap.put(svData.id(), svData);
            }
        }
        catch(IOException e)
        {
            GR_LOGGER.error("error reading vcf({}): {}", vcfFile, e.toString());
        }
    }

    private void compareVariants()
    {
        Set<SvData> matchedSvs = Sets.newHashSet();

        int diffCount = 0;


        for(SvData origSv : mOriginalSvData.values())
        {
            SvData newSv = mNewSvData.get(origSv.id());

            if(newSv == null)
            {
                writeDiffs(origSv, null, "NO_NEW", "", "");
                ++diffCount;
                continue;
            }

            matchedSvs.add(newSv);

            /*
            StringJoiner diffTypes = new StringJoiner(";");
            StringJoiner diffDetails = new StringJoiner(";");

            if(!findDiffs(origSv, newsSv, diffTypes, diffDetails))
                continue;

            writeDiffs(origSv, newsSv, diffTypes.toString(), diffDetails.toString());
            */

            boolean hasDiff = false;

            if(origSv.posStart() != newSv.posStart() || origSv.posEnd() != newSv.posEnd()
            || origSv.orientStart() != newSv.orientStart() || origSv.orientEnd() != newSv.orientEnd())
            {
                hasDiff = true;
                writeDiffs(origSv, newSv, "COORDS", makeSvCoords(origSv), makeSvCoords(newSv));
            }

            for(int se = SE_START; se <= SE_END; ++se)
            {
                if(origSv.isSgl() && se == SE_END)
                    continue;

                Breakend origStart = origSv.breakends()[se];
                Breakend newStart = newSv.breakends()[se];

                if(se == SE_START)
                {
                    Set<String> origFilters = origStart.Context.getFilters();
                    Set<String> newFilters = newStart.Context.getFilters();

                    if((!origFilters.isEmpty() || !newFilters.isEmpty())
                    && origFilters.size() != newFilters.size() || origFilters.stream().anyMatch(x -> !newFilters.contains(x)))
                    {
                        writeDiffs(origSv, newSv, "FILTERS", filtersStr(origFilters), filtersStr(newFilters));
                        hasDiff = true;
                    }
                }
            }

            if(hasDiff)
                ++diffCount;
        }

        for(SvData newSv : mNewSvData.values())
        {
            if(matchedSvs.contains(newSv))
                continue;

            ++diffCount;
            writeDiffs(null, newSv, "NO_ORIG", "", "");
        }

        GR_LOGGER.info("diffTotal({})", diffCount);
    }

    private boolean findDiffs(final SvData origSv, final SvData newSv, final StringJoiner diffTypes, final StringJoiner diffDetails)
    {
        if(origSv.posStart() != newSv.posStart() || origSv.posEnd() != newSv.posEnd()
        || origSv.orientStart() != newSv.orientStart() || origSv.orientEnd() != newSv.orientEnd())
        {
            diffTypes.add("COORDS");
            diffDetails.add(String.format("new=%s)", makeSvCoords(newSv)));
        }

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(origSv.isSgl() && se == SE_END)
                continue;

            Breakend origStart = origSv.breakends()[se];
            Breakend newStart = newSv.breakends()[se];

            if(se == SE_START)
            {
                Set<String> origFilters = origStart.Context.getFilters();
                Set<String> newFilters = newStart.Context.getFilters();

                if((!origFilters.isEmpty() || !newFilters.isEmpty())
                && origFilters.size() != newFilters.size() || origFilters.stream().anyMatch(x -> !newFilters.contains(x)))
                {
                    diffTypes.add("FILTERS");
                    diffDetails.add(String.format("orig=%s new=%s", filtersStr(origFilters), filtersStr(newFilters)));
                }
            }
        }

        return diffTypes.length() > 0;
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

            writer.write("SvId,Coords,DiffType,OrigValue,NewValue");
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

            mWriter.write(String.format("%s,%s,%s,%s,%s",
                    origSv != null ? origSv.id() : newSv.id(), coords,  diffType, origValue, newValue));

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

package com.hartwig.hmftools.esvee.output;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.WriteType;
import com.hartwig.hmftools.esvee.processor.VariantCall;

public class ResultsWriter
{
    private final SvConfig mConfig;
    private final BufferedWriter mVariantWriter;
    private final BamWriter mBamWriter;

    public ResultsWriter(final SvConfig config)
    {
        mConfig = config;
        mVariantWriter = initialiseVariantWriter();
        mBamWriter = new BamWriter(config);
    }

    public void close()
    {
        closeBufferedWriter(mVariantWriter);

        mBamWriter.close();
    }

    private BufferedWriter initialiseVariantWriter()
    {
        if(!mConfig.WriteTypes.contains(WriteType.BREAKEND_TSV))
            return null;

        if(mConfig.OutputDir == null)
            return null;

        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.outputFilename(WriteType.BREAKEND_TSV));

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("ChromosomeStart");
            sj.add("PositionStart");
            sj.add("ChromosomeEnd");
            sj.add("PositionEnd");
            sj.add("DescriptorStart");
            sj.add("DescriptorEnd");
            sj.add("MapQualStart");
            sj.add("MapQualEnd");
            sj.add("GermlineSupport");
            sj.add("SomaticSupport");
            sj.add("AssemblyClassification");
            sj.add("Filters");

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise variant writer: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeVariantAssemblyBamRecords(final List<VariantCall> variants)
    {
        mBamWriter.writeVariantAssemblyBamRecords(variants);
    }

    public synchronized void writeVariant(final VariantCall variant)
    {
        if(mVariantWriter == null)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(variant.LeftChromosome);
            sj.add(String.valueOf(variant.LeftPosition));
            // sj.add(String.valueOf(variant.)); // need orientation
            sj.add(variant.RightChromosome);
            sj.add(String.valueOf(variant.RightPosition));

            sj.add(variant.LeftDescriptor);
            sj.add(variant.RightDescriptor);

            sj.add(String.valueOf(variant.LeftMappingQuality));
            sj.add(String.valueOf(variant.RightMappingQuality));

            sj.add(String.valueOf(variant.germlineSupport()));
            sj.add(String.valueOf(variant.somaticSupport()));

            List<String> associatedAssemblies = variant.associatedAssemblies().stream()
                    .map(asm -> asm.Assembly).collect(Collectors.toList());

            // sj.add(associatedAssemblies); // required?
            sj.add(String.valueOf(variant.Classification));

            // List<String> supportingReads = List.of();
            // sj.add(supportingReads); // was never logged

            sj.add(buildFilters(variant));

            mVariantWriter.write(sj.toString());
            mVariantWriter.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise variant writer: {}", e.toString());
        }
    }

    private String buildFilters(final VariantCall variant)
    {
        boolean isLowOverhang = variant.overhang() < SvConstants.VCFLOWOVERHANGTHRESHOLD;
        boolean isLowQuality = variant.quality() < SvConstants.VCFLOWQUALITYTHRESHOLD;
        boolean isLowSupport = variant.supportingFragments().size() < SvConstants.MINREADSTOSUPPORTASSEMBLY;
        boolean isLikelyFalse = isLowSupport || (isLowOverhang && variant.discordantSupport() == 0) || isLowQuality;

        List<String> filters = new ArrayList<>();

        if(variant.isGermline())
            filters.add("GERMLINE");

        if(variant.associatedAssemblies().size() > 1)
            filters.add("MULTIPLE_ASSEMBLIES");

        if(isLowOverhang)
            filters.add("LOW_OVERHANG");

        if(isLowQuality)
            filters.add("LOW_QUALITY");

        if(isLowSupport)
            filters.add("LOW_SUPPORT");

        if(isLikelyFalse)
            filters.add("LIKELY_FALSE");

        return String.join(ITEM_DELIM, filters);
    }
}

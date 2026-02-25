package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.String.format;
import static java.lang.String.valueOf;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ALT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REF;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.FileType.SOMATICS;
import static com.hartwig.hmftools.wisp.purity.FileType.SOMATIC_PEAK;
import static com.hartwig.hmftools.wisp.purity.FileType.SUMMARY;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.addCommonFields;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.addCommonHeaderFields;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.RExecutor;
import com.hartwig.hmftools.wisp.purity.PurityConfig;
import com.hartwig.hmftools.wisp.purity.SampleData;

public final class SomaticWriter
{
    public static BufferedWriter initialiseVariantWriter(final PurityConfig config)
    {
        try
        {
            String fileName = config.formFilename(SOMATICS);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            addCommonHeaderFields(sj, config);

            sj.add(FLD_CHROMOSOME).add(FLD_POSITION).add(FLD_REF).add(FLD_ALT).add("IsProbe");
            sj.add("Filter").add("Tier").add("Type").add("TNC");
            sj.add("Mappability").add("SubclonalPerc").add("RepeatCount");
            sj.add("Gene").add("CodingEffect").add("Hotspot").add("Reported");
            sj.add("VCN").add("CopyNumber");
            sj.add("TumorDP").add("TumorAD");
            sj.add("SampleDP").add("SampleAD");

            sj.add("DualFilter").add("SampleDualDP").add("SampleDualAD").add("SampleQualPerAD");
            sj.add("SeqGcRatio").add("BqrErrorRate").add("DualBqrErrorRate").add("AvgReadDistance");

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to initialise variant output file: {}", e.toString());
            return null;
        }
    }

    protected static synchronized void writeVariant(
            final BufferedWriter writer, final PurityConfig config,
            final SampleData sampleData, final String sampleId, final SomaticVariant variant,
            final GenotypeFragments sampleFragData, final GenotypeFragments tumorData)
    {
        if(writer == null)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            addCommonFields(sj, config, sampleData, sampleId);

            sj.add(variant.Chromosome).add(valueOf(variant.Position)).add(variant.Ref).add(variant.Alt);
            sj.add(valueOf(variant.isProbeVariant()));

            List<FilterReason> filterReasons = Lists.newArrayList(variant.filterReasons());
            filterReasons.addAll(sampleFragData.filterReasons());

            String filtersStr = !filterReasons.isEmpty() ?
                    filterReasons.stream().map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM)) : PASS_FILTER;

            sj.add(filtersStr).add(variant.Tier.toString()).add(variant.Type.toString()).add(variant.TriNucContext);
            sj.add(format("%.2f", variant.Mappability)).add(format("%.2f", variant.SubclonalPerc)).add(valueOf(variant.RepeatCount));
            sj.add(variant.CanonicalGeneName).add(variant.CanonicalCodingEffect).add(valueOf(variant.Hotspot)).add(valueOf(variant.Reported));
            sj.add(format("%.2f", variant.VariantCopyNumber)).add(format("%.2f", variant.CopyNumber));

            sj.add(valueOf(tumorData.Depth)).add(valueOf(tumorData.AlleleCount));
            sj.add(valueOf(sampleFragData.Depth)).add(valueOf(sampleFragData.AlleleCount));

            String dualFiltersStr = !sampleFragData.dualFilterReasons().isEmpty() ?
                    sampleFragData.dualFilterReasons().stream().map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM)) : PASS_FILTER;

            sj.add(dualFiltersStr).add(valueOf(sampleFragData.UmiCounts.TotalDual)).add(valueOf(sampleFragData.UmiCounts.AlleleDual));
            sj.add(format("%.1f", sampleFragData.qualPerAlleleFragment()));
            sj.add(format("%.3f", variant.sequenceGcRatio()));
            sj.add(format("%4.3e", sampleFragData.bqrErrorRate()));
            sj.add(format("%4.3e", sampleFragData.dualBqrErrorRate()));
            sj.add(format("%.2f", sampleFragData.averageReadDistance()));

            writer.write(sj.toString());

            writer.newLine();
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to write output file: {}", e.toString());
            System.exit(1);
        }
    }

    public static boolean plotSomaticVafs(final String patientId, final String sampleId, final PurityConfig config)
    {
        try
        {
            String summaryFile = config.formFilename(SUMMARY);
            String somaticPeaksFile = config.formFilename(SOMATIC_PEAK);

            if(!Files.exists(Paths.get(summaryFile)) || !Files.exists(Paths.get(somaticPeaksFile)))
            {
                CT_LOGGER.warn("plots missing required files: summary({}) somatics({})", summaryFile, somaticPeaksFile);
                return false;
            }

            int runCode = RExecutor.executeFromClasspath(
                    "plots/SomaticVafPlot.R", patientId, sampleId, summaryFile, somaticPeaksFile, config.PlotDir);

            return runCode == 0;
        }
        catch(Exception e)
        {
            CT_LOGGER.error("failed to generate CN plot with R script: {}", e.toString());
            return false;
        }
    }
}

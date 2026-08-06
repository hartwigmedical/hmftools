package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CANDIDATE_PROBES_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CANDIDATE_TARGET_REGIONS_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.COVERED_REGIONS_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.COVERED_TARGET_REGIONS_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.DNA_OUTPUT_PREFIX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_STATS_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PANEL_PROBES_FILE_STEM;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.REJECTED_FEATURES_FILE_STEM;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.RNA_OUTPUT_PREFIX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_VARIANT_INFO_FILE_NAME;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.List;
import java.util.function.Function;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.panelbuilder.samplevariants.SampleVariants;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

// Writes all file output data. Per-panel probe output is delegated to a probe-output writer, one for the DNA panel and, when requested, one
// for the RNA panel. Sample variant info is DNA-panel only.
public class OutputWriter implements AutoCloseable
{
    private final ProbeOutputWriter mDnaOutput;
    // Null when RNA probes are not requested.
    @Nullable
    private final ProbeOutputWriter mRnaOutput;
    private final DelimFileWriter<SampleVariants.VariantInfo> mSampleVariantInfoTsvWriter;

    private enum SampleVariantInfoColumns
    {
        Variant,
        TargetType,
        FilterReason
    }

    private static final Logger LOGGER = LogManager.getLogger(OutputWriter.class);

    public OutputWriter(final String outputDir, @Nullable final String outputId, boolean verboseOutput, boolean rnaOutput) throws IOException
    {
        Function<String, String> outputFilePath = fileName ->
        {
            if(outputId != null)
            {
                fileName = outputId + "." + fileName;
            }
            return Paths.get(outputDir, fileName).toString();
        };

        // DNA and RNA produce the same set of per-panel files, distinguished only by file-name prefix.
        Function<String, String> dnaFilePath = fileName -> outputFilePath.apply(DNA_OUTPUT_PREFIX + fileName);
        Function<String, String> rnaFilePath = fileName -> outputFilePath.apply(RNA_OUTPUT_PREFIX + fileName);

        mDnaOutput = new ProbeOutputWriter(
                dnaFilePath, PANEL_PROBES_FILE_STEM, COVERED_TARGET_REGIONS_FILE_NAME, COVERED_REGIONS_FILE_NAME,
                REJECTED_FEATURES_FILE_STEM, CANDIDATE_TARGET_REGIONS_FILE_NAME, CANDIDATE_PROBES_FILE_NAME, GENE_STATS_FILE_NAME,
                verboseOutput, false);

        mRnaOutput = rnaOutput
                ? new ProbeOutputWriter(
                rnaFilePath, PANEL_PROBES_FILE_STEM, COVERED_TARGET_REGIONS_FILE_NAME, COVERED_REGIONS_FILE_NAME,
                REJECTED_FEATURES_FILE_STEM, CANDIDATE_TARGET_REGIONS_FILE_NAME, CANDIDATE_PROBES_FILE_NAME, GENE_STATS_FILE_NAME,
                verboseOutput, true)
                : null;

        mSampleVariantInfoTsvWriter = new DelimFileWriter<>(
                outputFilePath.apply(SAMPLE_VARIANT_INFO_FILE_NAME), SampleVariantInfoColumns.values(),
                OutputWriter::writeSampleVariantInfoRow);
    }

    public ProbeOutputWriter panelOutput()
    {
        return mDnaOutput;
    }

    // Null when RNA probes are not requested.
    @Nullable
    public ProbeOutputWriter rnaPanelOutput()
    {
        return mRnaOutput;
    }

    public void writeCandidateProbe(final Probe probe)
    {
        mDnaOutput.writeCandidateProbe(probe);
    }

    public void writeRnaCandidateProbe(final Probe probe)
    {
        if(mRnaOutput != null)
        {
            mRnaOutput.writeCandidateProbe(probe);
        }
    }

    public void writeSampleVariantInfos(final List<SampleVariants.VariantInfo> variantInfos)
    {
        variantInfos.forEach(mSampleVariantInfoTsvWriter::writeRow);
    }

    private static void writeSampleVariantInfoRow(final SampleVariants.VariantInfo variantInfo, DelimFileWriter.Row row)
    {
        row.set(SampleVariantInfoColumns.Variant, variantInfo.variant());
        row.set(SampleVariantInfoColumns.TargetType, variantInfo.targetType().name());
        row.set(SampleVariantInfoColumns.FilterReason, variantInfo.filterReason() == null ? "PASS" : variantInfo.filterReason().name());
    }

    @Override
    public void close() throws IOException
    {
        LOGGER.debug("Flushing and closing output files");

        mDnaOutput.close();
        if(mRnaOutput != null)
        {
            mRnaOutput.close();
        }
        mSampleVariantInfoTsvWriter.close();
    }
}

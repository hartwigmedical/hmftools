package com.hartwig.hmftools.panelbuilder;

import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CANDIDATE_PROBES_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.CANDIDATE_TARGET_REGIONS_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.COVERED_REGIONS_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.COVERED_TARGET_REGIONS_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.GENE_STATS_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PANEL_PROBES_FILE_STEM;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.REJECTED_FEATURES_FILE_STEM;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.RNA_CANDIDATE_TARGET_REGIONS_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.RNA_COVERED_REGIONS_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.RNA_COVERED_TARGET_REGIONS_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.RNA_GENE_STATS_FILE_NAME;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.RNA_PANEL_PROBES_FILE_STEM;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.RNA_REJECTED_FEATURES_FILE_STEM;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_VARIANT_INFO_FILE_NAME;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.panelbuilder.samplevariants.SampleVariants;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

// Writes all file output data. The per-panel probe output (probes, regions, rejections, gene stats) is delegated to a ProbeOutputWriter, one
// for the DNA panel and, when requested, one for the RNA panel. The verbose candidate probes and sample variant info are DNA-panel only.
public class OutputWriter implements AutoCloseable
{
    private final ProbeOutputWriter mDnaOutput;
    // Null when RNA probes are not requested.
    @Nullable
    private final ProbeOutputWriter mRnaOutput;
    @Nullable
    private final DelimFileWriter<Probe> mCandidateProbesTsvWriter;
    @Nullable
    private final ArrayList<Probe> mCandidateProbesBuffer;
    private final DelimFileWriter<SampleVariants.VariantInfo> mSampleVariantInfoTsvWriter;

    private enum CandidateProbesColumns
    {
        StartRegion,
        StartRegionOrient,
        MiddleSequence,
        EndRegion,
        EndRegionOrient,
        Sequence,
        TargetType,
        TargetExtra,
        QualityScore,
        GCContent,
        EvalCriteria
    }

    private enum SampleVariantInfoColumns
    {
        Variant,
        TargetType,
        FilterReason
    }

    private static final int CANDIDATE_PROBES_BUFFER_SIZE = 1_000_000;

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

        mDnaOutput = new ProbeOutputWriter(
                outputFilePath, PANEL_PROBES_FILE_STEM, COVERED_TARGET_REGIONS_FILE_NAME, COVERED_REGIONS_FILE_NAME,
                REJECTED_FEATURES_FILE_STEM, CANDIDATE_TARGET_REGIONS_FILE_NAME, GENE_STATS_FILE_NAME, false);

        mRnaOutput = rnaOutput
                ? new ProbeOutputWriter(
                outputFilePath, RNA_PANEL_PROBES_FILE_STEM, RNA_COVERED_TARGET_REGIONS_FILE_NAME, RNA_COVERED_REGIONS_FILE_NAME,
                RNA_REJECTED_FEATURES_FILE_STEM, RNA_CANDIDATE_TARGET_REGIONS_FILE_NAME, RNA_GENE_STATS_FILE_NAME, true)
                : null;

        if(verboseOutput)
        {
            mCandidateProbesTsvWriter = new DelimFileWriter<>(
                    outputFilePath.apply(CANDIDATE_PROBES_FILE_NAME), CandidateProbesColumns.values(),
                    OutputWriter::writeCandidateProbesRow);
            mCandidateProbesBuffer = new ArrayList<>(CANDIDATE_PROBES_BUFFER_SIZE);
        }
        else
        {
            mCandidateProbesTsvWriter = null;
            mCandidateProbesBuffer = null;
        }

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
        if(mCandidateProbesBuffer != null)
        {
            // Buffer probes to improve performance.
            mCandidateProbesBuffer.add(probe);
            checkFlushCandidateProbes(false);
        }
    }

    private void checkFlushCandidateProbes(boolean force)
    {
        if(mCandidateProbesBuffer != null)
        {
            if(mCandidateProbesBuffer.size() >= CANDIDATE_PROBES_BUFFER_SIZE || force)
            {
                writeCandidateProbes(mCandidateProbesBuffer);
                mCandidateProbesBuffer.clear();
            }
        }
    }

    private void writeCandidateProbes(final List<Probe> probes)
    {
        if(mCandidateProbesTsvWriter != null)
        {
            LOGGER.debug("Writing {} candidate probes to file", probes.size());
            probes.forEach(mCandidateProbesTsvWriter::writeRow);
        }
    }

    private static void writeCandidateProbesRow(final Probe probe, DelimFileWriter.Row row)
    {
        BasicProbeLayout layout = BasicProbeLayout.from(probe.definition());
        ChrBaseRegion start = layout.startRegion();
        Orientation startOrientation = layout.startOrientation();
        ChrBaseRegion end = layout.endRegion();
        Orientation endOrientation = layout.endOrientation();
        row.setOrNull(CandidateProbesColumns.StartRegion, start == null ? null : start.toString());
        row.setOrNull(CandidateProbesColumns.StartRegionOrient, startOrientation == null ? null : startOrientation.asChar());
        row.setOrNull(CandidateProbesColumns.MiddleSequence, layout.insertSequence());
        row.setOrNull(CandidateProbesColumns.EndRegion, end == null ? null : end.toString());
        row.setOrNull(CandidateProbesColumns.EndRegionOrient, endOrientation == null ? null : endOrientation.asChar());
        row.setOrNull(CandidateProbesColumns.Sequence, probe.sequence());
        row.set(CandidateProbesColumns.TargetType, probe.metadata().type().name());
        row.set(CandidateProbesColumns.TargetExtra, probe.metadata().extraInfo());
        row.setOrNull(CandidateProbesColumns.QualityScore, probe.qualityScore());
        row.setOrNull(CandidateProbesColumns.GCContent, probe.gcContent());
        row.set(CandidateProbesColumns.EvalCriteria, requireNonNull(probe.evaluationCriteria()).toString());
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

        if(mCandidateProbesTsvWriter != null)
        {
            checkFlushCandidateProbes(true);
            mCandidateProbesTsvWriter.close();
        }

        mSampleVariantInfoTsvWriter.close();
    }
}

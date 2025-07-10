package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class OutputWriter implements AutoCloseable
{
    private final DelimFileWriter<EvaluatedProbe> mPanelProbesWriter;
    private final DelimFileWriter<RejectedRegion> mRejectedRegionsWriter;
    private final DelimFileWriter<EvaluatedProbe> mCandidateProbesWriter;

    private final List<EvaluatedProbe> mCandidateProbesBuffer;

    private static final int CANDIDATE_PROBES_BUFFER_SIZE = 1_000_000;

    private static final String FLD_SEQUENCE = "Sequence";
    private static final String FLD_QUALITY_SCORE = "QualityScore";
    private static final String FLD_GC_CONTENT = "GCContent";
    private static final String FLD_SOURCE_TYPE = "SourceType";
    private static final String FLD_SOURCE_EXTRA_INFO = "SourceExtra";
    private static final String FLD_REJECT_REASON = "RejectReason";
    private static final String FLD_TARGET_START = "TargetPositionStart";
    private static final String FLD_TARGET_END = "TargetPositionEnd";
    private static final String FLD_EVAL_CRITERIA = "EvalCriteria";

    private static final List<String> PANEL_PROBES_COLUMNS = List.of(
            FLD_CHROMOSOME, FLD_POSITION_START, FLD_POSITION_END, FLD_SEQUENCE,
            FLD_SOURCE_TYPE, FLD_SOURCE_EXTRA_INFO,
            FLD_QUALITY_SCORE, FLD_GC_CONTENT);

    private static final List<String> REJECTED_REGIONS_COLUMNS = List.of(
            FLD_CHROMOSOME, FLD_POSITION_START, FLD_POSITION_END,
            FLD_SOURCE_TYPE, FLD_SOURCE_EXTRA_INFO,
            FLD_REJECT_REASON);

    private static final List<String> CANDIDATE_PROBES_COLUMNS = List.of(
            FLD_CHROMOSOME, FLD_POSITION_START, FLD_POSITION_END,
            FLD_TARGET_START, FLD_TARGET_END,
            FLD_SOURCE_TYPE, FLD_SOURCE_EXTRA_INFO,
            FLD_QUALITY_SCORE, FLD_GC_CONTENT,
            FLD_EVAL_CRITERIA, FLD_REJECT_REASON);

    private static final Logger LOGGER = LogManager.getLogger(OutputWriter.class);

    public OutputWriter(final String panelProbesFile, final String rejectedRegionsFile, final String candidateProbesFile)
    {
        mPanelProbesWriter = new DelimFileWriter<>(panelProbesFile, PANEL_PROBES_COLUMNS, OutputWriter::writePanelProbesRow);
        mRejectedRegionsWriter =
                new DelimFileWriter<>(rejectedRegionsFile, REJECTED_REGIONS_COLUMNS, OutputWriter::writeRejectedRegionsRow);
        mCandidateProbesWriter =
                new DelimFileWriter<>(candidateProbesFile, CANDIDATE_PROBES_COLUMNS, OutputWriter::writeCandidateProbesRow);
        mCandidateProbesBuffer = new ArrayList<>(CANDIDATE_PROBES_BUFFER_SIZE);
    }

    public void writePanelProbes(final List<EvaluatedProbe> probes)
    {
        LOGGER.debug("Writing {} panel probes to file", probes.size());
        probes.forEach(mPanelProbesWriter::writeRow);
    }

    private static void writePanelProbesRow(final EvaluatedProbe probe, DelimFileWriter.Row row)
    {
        row.set(FLD_CHROMOSOME, probe.candidate().probeRegion().chromosome());
        row.set(FLD_POSITION_START, probe.candidate().probeRegion().start());
        row.set(FLD_POSITION_END, probe.candidate().probeRegion().end());
        row.set(FLD_SEQUENCE, probe.sequence());
        row.set(FLD_QUALITY_SCORE, probe.qualityScore());
        row.set(FLD_GC_CONTENT, probe.gcContent());
        row.set(FLD_SOURCE_TYPE, probe.candidate().source().type().name());
        row.set(FLD_SOURCE_EXTRA_INFO, probe.candidate().source().extra());
    }

    public void writeRejectedRegions(final List<RejectedRegion> regions)
    {
        LOGGER.debug("Writing {} rejected regions to file", regions.size());
        regions.forEach(mRejectedRegionsWriter::writeRow);
    }

    private static void writeRejectedRegionsRow(final RejectedRegion region, DelimFileWriter.Row row)
    {
        row.set(FLD_CHROMOSOME, region.baseRegion().chromosome());
        row.set(FLD_POSITION_START, region.baseRegion().start());
        row.set(FLD_POSITION_END, region.baseRegion().end());
        row.set(FLD_SOURCE_TYPE, region.source().type().name());
        row.set(FLD_SOURCE_EXTRA_INFO, region.source().extra());
        row.set(FLD_REJECT_REASON, region.reason());
    }

    public void writeCandidateProbe(final EvaluatedProbe probe)
    {
        // Buffer probes to improve performance.
        mCandidateProbesBuffer.add(probe);
        checkFlushCandidateProbes(false);
    }

    private void checkFlushCandidateProbes(boolean force)
    {
        if(mCandidateProbesBuffer.size() >= CANDIDATE_PROBES_BUFFER_SIZE || force)
        {
            writeCandidateProbes(mCandidateProbesBuffer);
            mCandidateProbesBuffer.clear();
        }
    }

    private void writeCandidateProbes(final List<EvaluatedProbe> probes)
    {
        LOGGER.debug("Writing {} candidate probes to file", probes.size());
        probes.forEach(mCandidateProbesWriter::writeRow);
    }

    private static void writeCandidateProbesRow(final EvaluatedProbe probe, DelimFileWriter.Row row)
    {
        // To save disk space, don't write these fields:
        //   - Probe end position (implied from start and probe length)
        //   - Base sequence (implied from probe region)
        row.set(FLD_CHROMOSOME, probe.candidate().probeRegion().chromosome());
        row.set(FLD_POSITION_START, probe.candidate().probeRegion().start());
        row.set(FLD_TARGET_START, probe.candidate().targetRegion().start());
        row.set(FLD_TARGET_END, probe.candidate().targetRegion().end());
        row.set(FLD_SOURCE_TYPE, probe.candidate().source().type().name());
        row.set(FLD_SOURCE_EXTRA_INFO, probe.candidate().source().extra());
        row.setOrNull(FLD_QUALITY_SCORE, probe.qualityScore());
        row.setOrNull(FLD_GC_CONTENT, probe.gcContent());
        row.setOrNull(FLD_EVAL_CRITERIA, probe.criteria().toString());
        row.setOrNull(FLD_REJECT_REASON, probe.rejectionReason());
    }

    @Override
    public void close()
    {
        checkFlushCandidateProbes(true);

        mPanelProbesWriter.close();
        mRejectedRegionsWriter.close();
        mCandidateProbesWriter.close();
    }
}

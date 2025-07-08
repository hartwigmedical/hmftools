package com.hartwig.hmftools.geneutils.paneldesign;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;

import java.util.List;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.logging.log4j.util.BiConsumer;

public class DataWriter
{
    private static final String FLD_BASE_SEQUENCE = "Sequence";
    private static final String FLD_QUALITY_SCORE = "QualityScore";
    private static final String FLD_GC_CONTENT = "GCContent";
    private static final String FLD_SOURCE_TYPE = "SourceType";
    private static final String FLD_SOURCE_EXTRA_INFO = "SourceExtra";
    private static final String FLD_REJECT_REASON = "Reason";

    public static void writePanelProbes(final String filePath, final Stream<EvaluatedProbe> probes)
    {
        List<String> columns = List.of(
                FLD_CHROMOSOME, FLD_POSITION_START, FLD_POSITION_END,
                FLD_BASE_SEQUENCE, FLD_QUALITY_SCORE, FLD_GC_CONTENT,
                FLD_SOURCE_TYPE, FLD_SOURCE_EXTRA_INFO);
        BiConsumer<EvaluatedProbe, DelimFileWriter.Row> rowWriter = (probe, row) ->
        {
            row.set(FLD_CHROMOSOME, probe.candidate().probeRegion().chromosome());
            row.set(FLD_POSITION_START, probe.candidate().probeRegion().start());
            row.set(FLD_POSITION_END, probe.candidate().probeRegion().end());
            row.set(FLD_BASE_SEQUENCE, probe.sequence());
            row.set(FLD_QUALITY_SCORE, probe.qualityScore());
            row.set(FLD_GC_CONTENT, probe.gcContent());
            row.set(FLD_SOURCE_TYPE, probe.candidate().source().type().name());
            row.set(FLD_SOURCE_EXTRA_INFO, probe.candidate().source().extra());
        };
        try(DelimFileWriter<EvaluatedProbe> writer = new DelimFileWriter<>(filePath, columns, rowWriter))
        {
            probes.forEach(writer::writeRow);
        }
    }

    public static void writeRejectedRegions(final String filePath, final Stream<RejectedRegion> regions)
    {
        List<String> columns =
                List.of(FLD_CHROMOSOME, FLD_POSITION_START, FLD_POSITION_END, FLD_SOURCE_TYPE, FLD_SOURCE_EXTRA_INFO, FLD_REJECT_REASON);
        BiConsumer<RejectedRegion, DelimFileWriter.Row> rowWriter = (region, row) ->
        {
            row.set(FLD_CHROMOSOME, region.baseRegion().chromosome());
            row.set(FLD_POSITION_START, region.baseRegion().start());
            row.set(FLD_POSITION_END, region.baseRegion().end());
            row.set(FLD_SOURCE_TYPE, region.source().type().name());
            row.set(FLD_SOURCE_EXTRA_INFO, region.source().extra());
            row.set(FLD_REJECT_REASON, region.reason());
        };
        try(DelimFileWriter<RejectedRegion> writer = new DelimFileWriter<>(filePath, columns, rowWriter))
        {
            regions.forEach(writer::writeRow);
        }
    }
}

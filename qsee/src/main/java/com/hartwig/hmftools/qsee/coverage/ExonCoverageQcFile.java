package com.hartwig.hmftools.qsee.coverage;

import static com.hartwig.hmftools.common.metrics.GeneDepthFile.depthRangeFieldName;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SAMPLE_ID;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.COL_SAMPLE_TYPE;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.DECIMAL_FORMAT;
import static com.hartwig.hmftools.qsee.common.QseeFileCommon.MULTISAMPLE_SAMPLE_ID;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.metrics.DepthRange;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.qsee.QseeConfig;
import com.hartwig.hmftools.qsee.common.QseeFileCommon;
import com.hartwig.hmftools.qsee.common.SampleType;

public class ExonCoverageQcFile
{
    private static final String COL_GENE_NAME = "GeneName";
    private static final String COL_EXON_RANK = "ExonRank";
    private static final String COL_MEDIAN_DEPTH = "MedianDepth";
    private static final String COL_LOD_VAF = "LodVaf";
    private static final String COL_VARIANT_DETECTABILITY_STATUS = "VariantDetectabilityStatus";

    public static String generateFilename(QseeConfig config)
    {
        String sampleId = config.isSinglePatient() ?
                config.getSampleIds(SampleType.TUMOR).get(0) :
                MULTISAMPLE_SAMPLE_ID;

        return QseeFileCommon.generateFilename(config.OutputDir, sampleId, "exon_coverage_qc", config.OutputId, "tsv.gz");
    }

    public static void write(String filename, List<ExonCoverageQc> coverageQcPerExon, boolean writeLodDetails)
    {
        QC_LOGGER.info("Writing exon coverage QC to: {}", filename);

        // All exons share the same depth ranges
        List<DepthRange> depthRanges = writeLodDetails && !coverageQcPerExon.isEmpty()
                ? List.copyOf(coverageQcPerExon.get(0).depthCounts().keySet())
                : List.of();

        List<String> columns = new ArrayList<>(List.of(
                COL_SAMPLE_ID,
                COL_SAMPLE_TYPE,
                COL_GENE_NAME,
                COL_EXON_RANK,
                COL_MEDIAN_DEPTH,
                COL_LOD_VAF,
                COL_VARIANT_DETECTABILITY_STATUS
        ));

        depthRanges.forEach(depthRange -> columns.add(depthRangeFieldName(depthRange)));

        DelimFileWriter.write(filename, columns, coverageQcPerExon, (exonCoverageQc, row) -> {
            row.set(COL_SAMPLE_ID, exonCoverageQc.sampleId());
            row.set(COL_SAMPLE_TYPE, exonCoverageQc.sampleType().name());
            row.set(COL_GENE_NAME, exonCoverageQc.geneName());
            row.set(COL_EXON_RANK, exonCoverageQc.exonRank());
            row.set(COL_MEDIAN_DEPTH, exonCoverageQc.medianDepth(), DECIMAL_FORMAT);
            row.set(COL_LOD_VAF, exonCoverageQc.lodVaf(), DECIMAL_FORMAT);

            row.set(COL_VARIANT_DETECTABILITY_STATUS, exonCoverageQc.variantDetectabilityStatus().name());

            for(DepthRange depthRange : depthRanges)
            {
                row.set(depthRangeFieldName(depthRange), exonCoverageQc.depthCounts().get(depthRange));
            }
        });
    }
}

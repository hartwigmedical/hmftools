package com.hartwig.hmftools.qsee.lod;

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

public class GeneLodFile
{
    private static final String COL_GENE_NAME = "GeneName";
    private static final String COL_EXON_RANK = "ExonRank";
    private static final String COL_MEDIAN_DEPTH = "MedianDepth";
    private static final String COL_SAMPLE_VAF_LOD = "SampleVafLod";
    private static final String COL_TUMOR_VAF_LOD = "TumorVafLod";
    private static final String COL_VARIANT_DETECTABILITY_STATUS = "VariantDetectabilityStatus";

    public static String generateFilename(QseeConfig config)
    {
        String sampleId = config.isSinglePatient() ?
                config.getSampleIds(SampleType.TUMOR).get(0) :
                MULTISAMPLE_SAMPLE_ID;

        return QseeFileCommon.generateFilename(config.OutputDir, sampleId, "gene_lod", config.OutputId, "tsv.gz");
    }

    public static void write(String filename, List<GeneLod> geneLods, boolean writeLodDetails)
    {
        QC_LOGGER.info("Writing gene limit of detection to: {}", filename);

        // All exons share the same depth ranges
        List<DepthRange> depthRanges = writeLodDetails && !geneLods.isEmpty()
                ? List.copyOf(geneLods.get(0).depthCounts().keySet())
                : List.of();

        List<String> columns = new ArrayList<>(List.of(
                COL_SAMPLE_ID,
                COL_SAMPLE_TYPE,
                COL_GENE_NAME,
                COL_EXON_RANK,
                COL_MEDIAN_DEPTH,
                COL_SAMPLE_VAF_LOD,
                COL_TUMOR_VAF_LOD,
                COL_VARIANT_DETECTABILITY_STATUS
        ));

        depthRanges.forEach(depthRange -> columns.add(depthRangeFieldName(depthRange)));

        DelimFileWriter.write(filename, columns, geneLods, (geneLod, row) -> {
            row.set(COL_SAMPLE_ID, geneLod.sampleId());
            row.set(COL_SAMPLE_TYPE, geneLod.sampleType().name());
            row.set(COL_GENE_NAME, geneLod.geneName());
            row.set(COL_EXON_RANK, geneLod.exonRank());
            row.set(COL_MEDIAN_DEPTH, geneLod.medianDepth(), DECIMAL_FORMAT);
            row.set(COL_SAMPLE_VAF_LOD, geneLod.sampleVafLod(), DECIMAL_FORMAT);
            row.set(COL_TUMOR_VAF_LOD, geneLod.tumorVafLod(), DECIMAL_FORMAT);

            row.set(COL_VARIANT_DETECTABILITY_STATUS, geneLod.variantDetectabilityStatus().name());

            for(DepthRange depthRange : depthRanges)
            {
                row.set(depthRangeFieldName(depthRange), geneLod.depthCounts().get(depthRange));
            }
        });
    }
}

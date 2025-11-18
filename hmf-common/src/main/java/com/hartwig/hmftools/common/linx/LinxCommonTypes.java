package com.hartwig.hmftools.common.linx;

import java.io.File;

public final class LinxCommonTypes
{
    // super category for an SV or cluster
    public static final String SUPER_TYPE_SIMPLE = "SIMPLE";
    public static final String SUPER_TYPE_INSERTION = "INSERTION";
    public static final String SUPER_TYPE_RECIPROCAL = "RECIPROCAL";
    public static final String SUPER_TYPE_TEMPLATED_INSERTION = "TEMPLATED_INSERTION";
    public static final String SUPER_TYPE_COMPLEX = "COMPLEX";
    public static final String SUPER_TYPE_DOUBLE_MINUTE = "DOUBLE_MINUTE";
    public static final String SUPER_TYPE_INCOMPLETE = "INCOMPLETE";
    public static final String SUPER_TYPE_ARTIFACT = "ARTIFACT";

    // Visualiser file suffixes
    private static final String SV_FILE_EXTENSION = ".linx.vis_sv_data.tsv";
    private static final String SV_GERMLINE_FILE_EXTENSION = ".linx.germline.vis_sv_data.tsv";

    public static String generateVisSvFilename(final String basePath, final String sample, boolean isGermline)
    {
        return basePath + File.separator + sample + (isGermline ? SV_GERMLINE_FILE_EXTENSION : SV_FILE_EXTENSION);
    }

    private static final String SEGMENT_FILE_EXTENSION = ".linx.vis_segments.tsv";
    private static final String SEGMENT_GERMLINE_FILE_EXTENSION = ".linx.germline.vis_segments.tsv";

    public static String generateVisSegmentFilename(final String basePath, final String sample, boolean isGermline)
    {
        return basePath + File.separator + sample + (isGermline ? SEGMENT_GERMLINE_FILE_EXTENSION : SEGMENT_FILE_EXTENSION);
    }

    private static final String EXON_FILE_EXTENSION = ".linx.vis_gene_exon.tsv";
    private static final String EXON_GERMLINE_FILE_EXTENSION = ".linx.germline.vis_gene_exon.tsv";

    public static String generateVisExonFilename(final String basePath, final String sample, boolean isGermline)
    {
        return basePath + File.separator + sample + (isGermline ? EXON_GERMLINE_FILE_EXTENSION : EXON_FILE_EXTENSION);
    }

    private static final String PROTEIN_FILE_EXTENSION = ".linx.vis_protein_domain.tsv";
    private static final String PROTEIN_GERMLINE_FILE_EXTENSION = ".linx.germline.vis_protein_domain.tsv";

    public static String generateVisProteinFilename(final String basePath, final String sample, boolean isGermline)
    {
        return basePath + File.separator + sample + (isGermline ? PROTEIN_GERMLINE_FILE_EXTENSION : PROTEIN_FILE_EXTENSION);
    }

    private static final String FUSION_FILE_EXTENSION = ".linx.vis_fusion.tsv";
    private static final String FUSION_GERMLINE_FILE_EXTENSION = ".linx.germline.vis_fusion.tsv";

    public static String generateVisFusionFilename(final String basePath, final String sample, boolean isGermline)
    {
        return basePath + File.separator + sample + (isGermline ? FUSION_GERMLINE_FILE_EXTENSION : FUSION_FILE_EXTENSION);
    }

    private static final String CN_FILE_EXTENSION = ".linx.vis_copy_number.tsv";
    private static final String CN_GERMLINE_FILE_EXTENSION = ".linx.germline.vis_copy_number.tsv";

    public static String generateVisCopyNumberFilename(final String basePath, final String sample, boolean isGermline)
    {
        return basePath + File.separator + sample + (isGermline ? CN_GERMLINE_FILE_EXTENSION : CN_FILE_EXTENSION);
    }
}

package com.hartwig.hmftools.compar.resource;

import static com.hartwig.hmftools.common.fusion.KnownFusionData.FLD_CANCER_TYPES;
import static com.hartwig.hmftools.common.fusion.KnownFusionData.FLD_FIVE_GENE;
import static com.hartwig.hmftools.common.fusion.KnownFusionData.FLD_HIGH_IMPACT_PROM;
import static com.hartwig.hmftools.common.fusion.KnownFusionData.FLD_KNOWN_EXON_DOWN_RANGE;
import static com.hartwig.hmftools.common.fusion.KnownFusionData.FLD_KNOWN_EXON_TRANS;
import static com.hartwig.hmftools.common.fusion.KnownFusionData.FLD_KNOWN_EXON_UP_RANGE;
import static com.hartwig.hmftools.common.fusion.KnownFusionData.FLD_OVERRIDES;
import static com.hartwig.hmftools.common.fusion.KnownFusionData.FLD_PUB_MED;
import static com.hartwig.hmftools.common.fusion.KnownFusionData.FLD_THREE_GENE;
import static com.hartwig.hmftools.common.fusion.KnownFusionData.FLD_TYPE;
import static com.hartwig.hmftools.common.pathogenic.PathogenicSummaryFactory.CLNSIG;
import static com.hartwig.hmftools.common.pathogenic.PathogenicSummaryFactory.CLNSIGCONF;
import static com.hartwig.hmftools.common.purple.PurpleCommon.DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO;
import static com.hartwig.hmftools.common.purple.PurpleCommon.DEFAULT_DRIVER_HET_DELETION_THRESHOLD;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ALT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REF;
import static com.hartwig.hmftools.common.variant.pon.PonCache.FLD_MAX_READ_COUNT;
import static com.hartwig.hmftools.common.variant.pon.PonCache.FLD_MULTI_PON_STATUS;
import static com.hartwig.hmftools.common.variant.pon.PonCache.FLD_SAMPLE_COUNT;
import static com.hartwig.hmftools.common.variant.pon.PonCache.FLD_TOTAL_READ_COUNT;
import static com.hartwig.hmftools.compar.resource.ResourceField.ValueType.BOOLEAN;
import static com.hartwig.hmftools.compar.resource.ResourceField.ValueType.DOUBLE;
import static com.hartwig.hmftools.compar.resource.ResourceField.ValueType.GERMLINE_REPORTING;
import static com.hartwig.hmftools.compar.resource.ResourceField.ValueType.STRING;
import static com.hartwig.hmftools.compar.resource.ResourceField.compared;
import static com.hartwig.hmftools.compar.resource.ResourceField.displayOnly;
import static com.hartwig.hmftools.compar.resource.ResourceType.CLINVAR;
import static com.hartwig.hmftools.compar.resource.ResourceType.DRIVER_GENE_PANEL;
import static com.hartwig.hmftools.compar.resource.ResourceType.KNOWN_FUSION_DATA;
import static com.hartwig.hmftools.compar.resource.ResourceType.SAGE_PON;

import java.util.List;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.jetbrains.annotations.Nullable;

public class ResourceDefinition
{
    public final ResourceType Type;
    public final Pattern FileNamePattern;
    public final List<String> KeyColumns;
    @Nullable
    public final String InfoColumn;
    @Nullable
    public final String StripChrPrefixColumn;
    public final List<ResourceField> Fields;

    private static final String VCF_CHROM = "CHROM";
    private static final String VCF_POS = "POS";
    private static final String VCF_REF = "REF";
    private static final String VCF_ALT = "ALT";
    private static final String VCF_INFO = "INFO";

    private static final String BED_NAME = "name";
    private static final String BED_SCORE = "score";

    public static final List<String> BEDPE_COLUMNS = List.of(
            "chrom1", "start1", "end1", "chrom2", "start2", "end2", BED_NAME, BED_SCORE, "strand1", "strand2");

    public static final List<String> BED_COLUMNS = List.of("chrom", "start", "end", BED_NAME, BED_SCORE, "strand");

    private static final List<String> BEDPE_KEY_COLUMNS = List.of(
            "chrom1", "start1", "end1", "chrom2", "start2", "end2", "strand1", "strand2");

    private ResourceDefinition(
            final ResourceType type, final String fileNamePattern, final List<String> keyColumns, @Nullable final String infoColumn,
            @Nullable final String stripChrPrefixColumn, final List<ResourceField> fields)
    {
        Type = type;
        FileNamePattern = Pattern.compile(fileNamePattern);
        KeyColumns = keyColumns;
        InfoColumn = infoColumn;
        StripChrPrefixColumn = stripChrPrefixColumn;
        Fields = fields;
    }

    public static ResourceDefinition create(final ResourceType type)
    {
        return switch(type)
        {
            case DRIVER_GENE_PANEL -> driverGenePanel();
            case KNOWN_FUSION_DATA -> knownFusionData();
            case KNOWN_FUSION_BEDPE -> delimited(
                    type, "known_fusions\\.(37|38)\\.bedpe", BEDPE_KEY_COLUMNS, List.of(displayOnly(BED_NAME), displayOnly(BED_SCORE)));
            case SOMATIC_HOTSPOT_VCF -> hotspotVcf(type, "somatic", somaticHotspotFields());
            case SOMATIC_HOTSPOT_TSV -> hotspotTsv(type, "somatic", somaticHotspotFields());
            case GERMLINE_HOTSPOT_VCF -> hotspotVcf(type, "germline", germlineHotspotFields());
            case GERMLINE_HOTSPOT_TSV -> hotspotTsv(type, "germline", germlineHotspotFields());
            case CLINVAR -> clinvar();
            case SAGE_PON -> sagePon();
            case SV_PON -> delimited(
                    type, "sv_pon.*\\.bedpe\\.gz", BEDPE_KEY_COLUMNS, List.of(compared(BED_SCORE, STRING), displayOnly(BED_NAME)));
            case SGL_PON -> delimited(
                    type, "sgl_pon.*\\.bed\\.gz", List.of("chrom", "start", "end", "strand"),
                    List.of(compared(BED_SCORE, STRING), displayOnly(BED_NAME)));
        };
    }

    private static ResourceDefinition delimited(
            final ResourceType type, final String fileNamePattern, final List<String> keyColumns, final List<ResourceField> fields)
    {
        return new ResourceDefinition(type, fileNamePattern, keyColumns, null, null, fields);
    }

    private static ResourceDefinition vcf(
            final ResourceType type, final String fileNamePattern, final boolean stripChrPrefix, final List<ResourceField> fields)
    {
        return new ResourceDefinition(
                type, fileNamePattern, List.of(VCF_CHROM, VCF_POS, VCF_REF, VCF_ALT), VCF_INFO, stripChrPrefix ? VCF_CHROM : null, fields);
    }

    private static ResourceDefinition driverGenePanel()
    {
        List<ResourceField> fields = List.of(
                compared("reportMissense", BOOLEAN),
                compared("reportNonsense", BOOLEAN),
                compared("reportSplice", BOOLEAN),
                compared("reportDeletion", BOOLEAN),
                compared("reportHetDeletion", BOOLEAN, Boolean.FALSE.toString()),
                compared("reportLoh", BOOLEAN, Boolean.FALSE.toString()),
                compared("hetDeletionThreshold", DOUBLE, String.valueOf(DEFAULT_DRIVER_HET_DELETION_THRESHOLD)),
                compared("reportDisruption", BOOLEAN),
                compared("reportAmplification", BOOLEAN),
                compared("amplificationRatio", DOUBLE, String.valueOf(DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO)),
                compared("reportSomaticHotspot", BOOLEAN),
                compared("likelihoodType", STRING),
                compared("reportGermlineVariant", GERMLINE_REPORTING),
                compared("reportGermlineHotspot", GERMLINE_REPORTING),
                compared("reportGermlineDisruption", GERMLINE_REPORTING),
                compared("reportGermlineDeletion", GERMLINE_REPORTING),
                compared("reportGermlineAmplification", BOOLEAN, Boolean.FALSE.toString()),
                compared("additionalReportedTranscripts", STRING),
                compared("reportPGX", BOOLEAN),
                compared("reportHighExpression", BOOLEAN, Boolean.FALSE.toString()),
                compared("reportLowExpression", BOOLEAN, Boolean.FALSE.toString()),
                compared("reportNovelSpliceJunction", BOOLEAN, Boolean.FALSE.toString()));

        return delimited(DRIVER_GENE_PANEL, "(DriverGenePanel|driver_genes\\..+)\\.(37|38)\\.tsv", List.of("gene"), fields);
    }

    private static ResourceDefinition knownFusionData()
    {
        List<ResourceField> fields = List.of(
                compared(FLD_CANCER_TYPES, STRING),
                compared(FLD_PUB_MED, STRING),
                compared(FLD_KNOWN_EXON_TRANS, STRING),
                compared(FLD_KNOWN_EXON_UP_RANGE, STRING),
                compared(FLD_KNOWN_EXON_DOWN_RANGE, STRING),
                compared(FLD_HIGH_IMPACT_PROM, BOOLEAN),
                compared(FLD_OVERRIDES, STRING));

        return delimited(KNOWN_FUSION_DATA, "known_fusion_data\\.(37|38)\\.csv", List.of(FLD_TYPE, FLD_FIVE_GENE, FLD_THREE_GENE), fields);
    }

    private static List<ResourceField> somaticHotspotFields()
    {
        return List.of(displayOnly("input"), displayOnly("sources"));
    }

    private static List<ResourceField> germlineHotspotFields()
    {
        return List.of(displayOnly(CLNSIG), displayOnly(CLNSIGCONF), displayOnly("GENEINFO"), displayOnly("WHITELIST"));
    }

    private static ResourceDefinition hotspotVcf(final ResourceType type, final String hotspotType, final List<ResourceField> fields)
    {
        return vcf(type, "KnownHotspots\\." + hotspotType + "\\.(37|38)\\.vcf\\.gz", false, fields);
    }

    private static ResourceDefinition hotspotTsv(final ResourceType type, final String hotspotType, final List<ResourceField> fields)
    {
        return new ResourceDefinition(
                type, "KnownHotspots\\." + hotspotType + "\\.(37|38)\\.tsv\\.gz", List.of(FLD_CHROMOSOME, FLD_POSITION, FLD_REF, FLD_ALT),
                "Info", null, fields);
    }

    private static ResourceDefinition clinvar()
    {
        List<String> displayFields = List.of(
                "AF_ESP", "AF_EXAC", "AF_TGP", "ALLELEID", "CLNDN", "CLNDNINCL", "CLNDISDB", "CLNDISDBINCL", "CLNHGVS", "CLNREVSTAT",
                "CLNSIGINCL", "CLNVC", "CLNVCSO", "CLNVI", "DBVARID", "GENEINFO", "MC", "ORIGIN", "RS");

        List<ResourceField> fields = Stream.concat(
                Stream.of(compared(CLNSIG, STRING), compared(CLNSIGCONF, STRING)),
                displayFields.stream().map(ResourceField::displayOnly)).collect(Collectors.toList());

        return vcf(CLINVAR, "clinvar\\.(37|38)\\.vcf\\.gz", true, fields);
    }

    private static ResourceDefinition sagePon()
    {
        List<ResourceField> fields = List.of(
                compared(FLD_SAMPLE_COUNT, STRING),
                compared(FLD_MAX_READ_COUNT, STRING),
                compared(FLD_TOTAL_READ_COUNT, STRING),
                compared(FLD_MULTI_PON_STATUS, STRING),
                compared("SamplesCount", STRING),
                compared("MaxSampleReads", STRING),
                compared("TotalReads", STRING));

        return delimited(
                SAGE_PON, "(hmf_wgs_sage_pon[._].+|SageGermlinePon\\..+|pon_artefacts\\..+)\\.(37|38)\\.tsv\\.gz",
                List.of(FLD_CHROMOSOME, FLD_POSITION, FLD_REF, FLD_ALT), fields);
    }
}

package com.hartwig.hmftools.compar.resource;

import static com.hartwig.hmftools.compar.resource.ResourceType.CLINVAR;
import static com.hartwig.hmftools.compar.resource.ResourceType.DRIVER_GENE_PANEL;
import static com.hartwig.hmftools.compar.resource.ResourceType.GERMLINE_HOTSPOT_TSV;
import static com.hartwig.hmftools.compar.resource.ResourceType.GERMLINE_HOTSPOT_VCF;
import static com.hartwig.hmftools.compar.resource.ResourceType.KNOWN_FUSION_BEDPE;
import static com.hartwig.hmftools.compar.resource.ResourceType.KNOWN_FUSION_DATA;
import static com.hartwig.hmftools.compar.resource.ResourceType.SAGE_PON;
import static com.hartwig.hmftools.compar.resource.ResourceType.SGL_PON;
import static com.hartwig.hmftools.compar.resource.ResourceType.SOMATIC_HOTSPOT_TSV;
import static com.hartwig.hmftools.compar.resource.ResourceType.SOMATIC_HOTSPOT_VCF;
import static com.hartwig.hmftools.compar.resource.ResourceType.SV_PON;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class ResourceDefinitionTest
{
    @Test
    public void matchesReleaseFileNames()
    {
        assertType(DRIVER_GENE_PANEL,
                "DriverGenePanel.37.tsv", "DriverGenePanel.38.tsv", "driver_genes.msk.37.tsv", "driver_genes.omico-avenio.37.tsv",
                "driver_genes.pmh-nhlk.37.tsv", "driver_genes.pmh-panel-v2.37.tsv", "driver_genes.tso500.37.tsv",
                "driver_genes.oncopanel-v3.38.tsv", "driver_genes.tso500.38.tsv");

        assertType(KNOWN_FUSION_DATA, "known_fusion_data.37.csv", "known_fusion_data.38.csv");
        assertType(KNOWN_FUSION_BEDPE, "known_fusions.37.bedpe", "known_fusions.38.bedpe");

        assertType(SOMATIC_HOTSPOT_VCF, "KnownHotspots.somatic.37.vcf.gz", "KnownHotspots.somatic.38.vcf.gz");
        assertType(SOMATIC_HOTSPOT_TSV, "KnownHotspots.somatic.37.tsv.gz", "KnownHotspots.somatic.38.tsv.gz");
        assertType(GERMLINE_HOTSPOT_VCF, "KnownHotspots.germline.37.vcf.gz", "KnownHotspots.germline.38.vcf.gz");
        assertType(GERMLINE_HOTSPOT_TSV, "KnownHotspots.germline.37.tsv.gz", "KnownHotspots.germline.38.tsv.gz");
        assertType(CLINVAR, "clinvar.37.vcf.gz", "clinvar.38.vcf.gz");

        assertType(SAGE_PON,
                "hmf_wgs_sage_pon_1000.38.tsv.gz", "hmf_wgs_sage_pon.ill_sbx.38.tsv.gz", "hmf_wgs_sage_pon.ill_ult.38.tsv.gz",
                "pon_artefacts.msk.37.tsv.gz", "pon_artefacts.omico-avenio.37.tsv.gz", "pon_artefacts.pmh-nhlk.37.tsv.gz",
                "pon_artefacts.pmh-panel-v2.37.tsv.gz", "pon_artefacts.tso500.37.tsv.gz", "pon_artefacts.oncopanel-v3.38.tsv.gz",
                "pon_artefacts.tso500.38.tsv.gz", "SageGermlinePon.1000x.37.tsv.gz");

        assertType(SV_PON,
                "sv_pon.37.bedpe.gz", "sv_pon.38.bedpe.gz", "sv_pon_artefacts.oncopanel-v3.38.bedpe.gz",
                "sv_pon_artefacts.tso500.38.bedpe.gz", "sv_pon_artefacts.pmh-panel-v2.37.bedpe.gz", "sv_pon_hmf_ill_sbx.38.bedpe.gz");

        assertType(SGL_PON,
                "sgl_pon.37.bed.gz", "sgl_pon.38.bed.gz", "sgl_pon_artefacts.oncopanel-v3.38.bed.gz",
                "sgl_pon_artefacts.tso500.38.bed.gz", "sgl_pon_artefacts.pmh-panel-v2.37.bed.gz", "sgl_pon_hmf_ill_sbx.38.bed.gz");

        assertType(null, "DriverGenePanel.38.tsv.dvc", "clinvar.38.vcf.gz.tbi", "known_fusion_data.38.csv.md5", "README.md");
    }

    private static void assertType(@Nullable final ResourceType expectedType, final String... fileNames)
    {
        List<ResourceType> expected = expectedType != null ? List.of(expectedType) : Collections.emptyList();

        for(String fileName : fileNames)
        {
            List<ResourceType> matched = Arrays.stream(ResourceType.values())
                    .filter(x -> ResourceDefinition.create(x).FileNamePattern.matcher(fileName).matches())
                    .collect(Collectors.toList());

            assertEquals(fileName, expected, matched);
        }
    }
}

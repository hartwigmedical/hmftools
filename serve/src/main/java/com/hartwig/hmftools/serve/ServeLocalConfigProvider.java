package com.hartwig.hmftools.serve;

import java.io.IOException;
import java.net.InetAddress;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ServeLocalConfigProvider {

    private static final Logger LOGGER = LogManager.getLogger(ServeLocalConfigProvider.class);

    private ServeLocalConfigProvider() {
    }

    @NotNull
    public static ServeConfig create() throws IOException {
        String hostname = InetAddress.getLocalHost().getHostName();
        LOGGER.info("Generating config for '{}'", hostname);

        ImmutableServeConfig.Builder builder = ImmutableServeConfig.builder()
                .useVicc(false)
                .useIclusion(false)
                .useCkb(false)
                .useActin(false)
                .useDocm(false)
                .docmTsv(Strings.EMPTY)
                .useHartwigCurated(false)
                .hartwigCuratedTsv(Strings.EMPTY)
                .useHartwigCohort(false)
                .hartwigCohortTsv(Strings.EMPTY)
                .skipHotspotResolving(true);

        // Datastore is a shared vm with fixed paths.
        if (hostname.toLowerCase().equals("serve-vm")) {
            builder.viccJson("/data/common/dbs/serve/static_sources/vicc/all.json");
            builder.ckbDir("/data/common/dbs/ckb/210514_flex_dump");
            builder.ckbFilterTsv("/data/common/dbs/serve/curation/ckb_filters.tsv");
            builder.iClusionTrialTsv("/data/common/dbs/iclusion/iclusion_trials_prod.tsv");
            builder.actinTrailTsv("");

            builder.outputDir(System.getProperty("user.home") + "/tmp/serve");
            builder.missingDoidsMappingTsv("/data/common/dbs/serve/curation/missing_doids_mapping.tsv");
            builder.driverGene37Tsv("/data/common/dbs/driver_gene_panel/DriverGenePanel.37.tsv");
            builder.driverGene38Tsv("/data/common/dbs/driver_gene_panel/DriverGenePanel.38.tsv");
            builder.knownFusion37File("/data/common/dbs/fusions/known_fusion_data.37_v3.csv");
            builder.knownFusion38File("/data/common/dbs/fusions/known_fusion_data.38_v3.csv");
            builder.refGenome37FastaFile("/data/common/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta");
            builder.refGenome38FastaFile(
                    "/data/common/refgenomes/Homo_sapiens.GRCh38.no.alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna");
            builder.refGenome37To38Chain("/data/common/refgenomes/liftover/hg19ToHg38.over.chain");
            builder.refGenome38To37Chain("/data/common/refgenomes/liftover/hg38ToHg19.over.chain");
        } else {
            // Assume we run on local machine with fixed paths
            builder.viccJson(System.getProperty("user.home") + "/hmf/serve/static_sources/vicc/all.json");
            builder.ckbDir(System.getProperty("user.home") + "/hmf/serve/ckb");
            builder.ckbFilterTsv(System.getProperty("user.home") + "/hmf/serve/curation/ckb_filters.tsv");
            builder.iClusionTrialTsv(System.getProperty("user.home") + "/hmf/serve/iclusion/iclusion_trials_prod.tsv");
            builder.actinTrailTsv("");

            builder.outputDir(System.getProperty("user.home") + "/hmf/tmp/serve");
            builder.missingDoidsMappingTsv(System.getProperty("user.home") + "/hmf/serve/curation/missing_doids_mapping.tsv");
            builder.driverGene37Tsv(System.getProperty("user.home") + "/hmf/driver_gene_panel/DriverGenePanel.37.tsv");
            builder.driverGene38Tsv(System.getProperty("user.home") + "/hmf/driver_gene_panel/DriverGenePanel.38.tsv");
            builder.knownFusion37File(System.getProperty("user.home") + "/hmf/fusions/known_fusion_data.37.csv");
            builder.knownFusion38File(System.getProperty("user.home") + "/hmf/fusions/known_fusion_data.38.csv");
            builder.refGenome37FastaFile(
                    System.getProperty("user.home") + "/hmf/refgenomes/grch37/Homo_sapiens.GRCh37.GATK.illumina.fasta");
            builder.refGenome38FastaFile(
                    System.getProperty("user.home") + "/hmf/refgenomes/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna");
            builder.refGenome37To38Chain(System.getProperty("user.home") + "/hmf/refgenomes/liftover/hg19ToHg38.over.chain");
            builder.refGenome38To37Chain(System.getProperty("user.home") + "/hmf/refgenomes/liftover/hg38ToHg19.over.chain");
        }

        return builder.build();
    }
}

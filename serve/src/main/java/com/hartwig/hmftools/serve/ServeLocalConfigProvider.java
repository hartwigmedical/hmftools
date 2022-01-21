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

        // SERVE-VM is a shared vm tailored for running SERVE on GCP.
        if (hostname.toLowerCase().contains("serve-vm")) {
            builder.viccJson("/data/resources/crunch/serve/vicc/all.json");
            builder.ckbDir("/data/resources/custom/ckb/latest");
            builder.ckbFilterTsv("/data/resources/crunch/serve/curation/ckb_filters.tsv");
            builder.iClusionTrialTsv("/data/resources/crunch/serve/iclusion/iclusion_trials_prod.tsv");
            builder.iClusionFilterTsv("/data/resources/crunch/serve/curation/iclusion_filters.tsv");
            builder.actinTrialTsv("/data/resources/crunch/serve/actin/actin_knowledgebase.tsv");
            builder.actinFilterTsv("/data/resources/crunch/serve/curation/actin_filters.tsv");

            builder.outputDir(System.getProperty("user.home") + "/tmp/serve");
            builder.missingDoidsMappingTsv("/data/resources/crunch/serve/curation/public_missing_doids_mapping.tsv");
            builder.ensemblDataDir37("/data/resources/public/ensembl_data_cache/37");
            builder.ensemblDataDir38("/data/resources/public/ensembl_data_cache/38");
            builder.driverGene37Tsv("/data/resources/public/gene_panel/37/DriverGenePanel.37.tsv");
            builder.driverGene38Tsv("/data/resources/public/gene_panel/38/DriverGenePanel.38.tsv");
            builder.knownFusion37File("/data/resources/public/fusions/37/known_fusion_data.37.csv");
            builder.knownFusion38File("/data/resources/public/fusions/38/known_fusion_data.38.csv");
            builder.refGenome37FastaFile("/data/resources/bucket/reference_genome/37/Homo_sapiens.GRCh37.GATK.illumina.fasta");
            builder.refGenome38FastaFile("/data/resources/bucket/reference_genome/38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna");
            builder.refGenome37To38Chain("/data/resources/crunch/reference_genome_liftover/hg19ToHg38.over.chain");
            builder.refGenome38To37Chain("/data/resources/crunch/reference_genome_liftover/hg38ToHg19.over.chain");
        } else {
            // Assume we run on local machine with fixed paths
            String baseDir = System.getProperty("user.home") + "/hmf/";
            builder.viccJson(baseDir + "serve/static_sources/vicc/all.json");
            builder.ckbDir(baseDir + "serve/ckb");
            builder.ckbFilterTsv(baseDir + "serve/curation/ckb_filters.tsv");
            builder.iClusionTrialTsv(baseDir + "serve/iclusion/iclusion_trials_prod.tsv");
            builder.iClusionFilterTsv(baseDir + "serve/actin/iclusion_filter.tsv");
            builder.actinTrialTsv(baseDir + "serve/actin/actin_knowledgebase.tsv");
            builder.actinFilterTsv(baseDir + "serve/actin/actin_filter.tsv");

            builder.outputDir(baseDir + "tmp/serve");
            builder.missingDoidsMappingTsv(baseDir + "serve/curation/missing_doids_mapping.tsv");

            builder.ensemblDataDir37(baseDir + "repos/common-resources-public/ensembl_data_cache/37");
            builder.ensemblDataDir38(baseDir + "repos/common-resources-public/ensembl_data_cache/38");
            builder.driverGene37Tsv(baseDir + "repos/common-resources-public/gene_panel/37/DriverGenePanel.37.tsv");
            builder.driverGene38Tsv(baseDir + "repos/common-resources-public/gene_panel/38/DriverGenePanel.38.tsv");
            builder.knownFusion37File(baseDir + "repos/common-resources-public/fusions/37/known_fusion_data.37.csv");
            builder.knownFusion38File(baseDir + "repos/common-resources-public/fusions/38/known_fusion_data.38.csv");
            builder.refGenome37FastaFile(baseDir + "refgenomes/grch37/Homo_sapiens.GRCh37.GATK.illumina.fasta");
            builder.refGenome38FastaFile(baseDir + "refgenomes/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna");
            builder.refGenome37To38Chain(baseDir + "refgenomes/liftover/hg19ToHg38.over.chain");
            builder.refGenome38To37Chain(baseDir + "refgenomes/liftover/hg38ToHg19.over.chain");
        }

        return builder.build();
    }
}

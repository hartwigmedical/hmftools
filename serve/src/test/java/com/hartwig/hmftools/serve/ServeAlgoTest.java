package com.hartwig.hmftools.serve;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.serve.curation.DoidLookupTestFactory;
import com.hartwig.hmftools.serve.refgenome.RefGenomeManagerFactory;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ServeAlgoTest {

    private static final String VICC_JSON = Resources.getResource("vicc/empty.vicc.json").getPath();
    private static final String ACTIN_TRIAL_TSV = Resources.getResource("actin/example.tsv").getPath();
    private static final String ICLUSION_TRIAL_TSV = Resources.getResource("iclusion/empty.iclusion.tsv").getPath();
    private static final String CKB_DIR = Resources.getResource("ckb_data").getPath();
    private static final String CKB_FILTER_TSV = Resources.getResource("ckb_filter/ckb_filters.tsv").getPath();
    private static final String DOCM_TSV = Resources.getResource("docm/example.tsv").getPath();
    private static final String HARTWIG_CURATED_TSV = Resources.getResource("hartwig/example.tsv").getPath();
    private static final String HARTWIG_COHORT_TSV = Resources.getResource("hartwig/example.tsv").getPath();

    private static final String REF_GENOME_37_FASTA_FILE = Resources.getResource("refgenome/v37/ref.fasta").getPath();
    private static final String REF_GENOME_38_FASTA_FILE = Resources.getResource("refgenome/v38/ref.fasta").getPath();
    private static final String REF_GENOME_37_TO_38_CHAIN = Resources.getResource("refgenome/liftover/V37ToV38.over.chain").getPath();
    private static final String REF_GENOME_38_TO_37_CHAIN = Resources.getResource("refgenome/liftover/V38ToV37.over.chain").getPath();

    private static final String DRIVER_GENE_37_TSV = Resources.getResource("driver_gene_panel/driver_gene_panel.37.tsv").getPath();
    private static final String DRIVER_GENE_38_TSV = Resources.getResource("driver_gene_panel/driver_gene_panel.38.tsv").getPath();

    private static final String KNOWN_FUSION_37_FILE = Resources.getResource("known_fusion_data/known_fusion_data.37.csv").getPath();
    private static final String KNOWN_FUSION_38_FILE = Resources.getResource("known_fusion_data/known_fusion_data.38.csv").getPath();

    @Test
    public void canRunServeAlgo() throws IOException {
        ServeConfig config = algoBuilder().useVicc(true)
                .viccJson(VICC_JSON)
                .addViccSources(ViccSource.CIVIC, ViccSource.CGI)
                .useIclusion(true)
                .iClusionTrialTsv(ICLUSION_TRIAL_TSV)
                .useCkb(true)
                .ckbDir(CKB_DIR)
                .ckbFilterTsv(CKB_FILTER_TSV)
                .useActin(true)
                .actinTrialTsv(ACTIN_TRIAL_TSV)
                .useDocm(true)
                .docmTsv(DOCM_TSV)
                .useHartwigCohort(true)
                .hartwigCohortTsv(HARTWIG_COHORT_TSV)
                .useHartwigCurated(true)
                .hartwigCuratedTsv(HARTWIG_CURATED_TSV)
                .refGenome37FastaFile(REF_GENOME_37_FASTA_FILE)
                .refGenome38FastaFile(REF_GENOME_38_FASTA_FILE)
                .refGenome37To38Chain(REF_GENOME_37_TO_38_CHAIN)
                .refGenome38To37Chain((REF_GENOME_38_TO_37_CHAIN))
                .driverGene37Tsv(DRIVER_GENE_37_TSV)
                .driverGene38Tsv(DRIVER_GENE_38_TSV)
                .knownFusion37File(KNOWN_FUSION_37_FILE)
                .knownFusion38File(KNOWN_FUSION_38_FILE)
                .build();

        ServeAlgo algo = new ServeAlgo(RefGenomeManagerFactory.createFromServeConfig(config), DoidLookupTestFactory.dummy());

        assertNotNull(algo.run(config));
    }

    @NotNull
    private static ImmutableServeConfig.Builder algoBuilder() {
        return ImmutableServeConfig.builder().missingDoidsMappingTsv(Strings.EMPTY).outputDir(Strings.EMPTY).skipHotspotResolving(true);
    }
}
package com.hartwig.hmftools.serve;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.serve.curation.DoidLookupTestFactory;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolverFactory;
import com.hartwig.hmftools.serve.util.RefGenomeVersion;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ServeAlgoTest {

    private static final String VICC_JSON = Resources.getResource("vicc/empty.vicc.json").getPath();
    private static final String ICLUSION_TRIAL_TSV = Resources.getResource("iclusion/empty.iclusion.tsv").getPath();
    private static final String DOCM_TSV = Resources.getResource("docm/example.tsv").getPath();
    private static final String HARTWIG_CURATED_TSV = Resources.getResource("hartwig/example.tsv").getPath();
    private static final String HARTWIG_COHORT_TSV = Resources.getResource("hartwig/example.tsv").getPath();
    private static final String REF_GENOME_FASTA_FILE = Resources.getResource("refgenome/ref.fasta").getPath();

    @Test
    public void canRunServeAlgo() throws IOException {
        ServeAlgo algo = new ServeAlgo(Lists.newArrayList(),
                new KnownFusionCache(),
                Maps.newHashMap(),
                ProteinResolverFactory.dummy(),
                DoidLookupTestFactory.dummy());

        ServeConfig config = algoBuilder().useVicc(true)
                .viccJson(VICC_JSON)
                .addViccSources(ViccSource.CIVIC, ViccSource.CGI)
                .useIclusion(true)
                .iClusionTrialTsv(ICLUSION_TRIAL_TSV)
                .useDocm(true)
                .docmTsv(DOCM_TSV)
                .useHartwigCohort(true)
                .hartwigCohortTsv(HARTWIG_COHORT_TSV)
                .useHartwigCurated(true)
                .hartwigCuratedTsv(HARTWIG_CURATED_TSV)
                .refGenomeVersion(RefGenomeVersion.V37)
                .refGenomeFastaFile(REF_GENOME_FASTA_FILE)
                .build();

        assertNotNull(algo.run(config));
    }

    @NotNull
    private static ImmutableServeConfig.Builder algoBuilder() {
        return ImmutableServeConfig.builder()
                .missingDoidsMappingTsv(Strings.EMPTY)
                .driverGeneTsv(Strings.EMPTY)
                .knownFusionFile(Strings.EMPTY)
                .outputDir(Strings.EMPTY)
                .skipHotspotResolving(true);
    }
}
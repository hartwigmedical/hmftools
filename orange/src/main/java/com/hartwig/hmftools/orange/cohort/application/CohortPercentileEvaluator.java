package com.hartwig.hmftools.orange.cohort.application;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidEntry;
import com.hartwig.hmftools.common.doid.DoidParents;
import com.hartwig.hmftools.orange.cohort.datamodel.Evaluation;
import com.hartwig.hmftools.orange.cohort.datamodel.Observation;
import com.hartwig.hmftools.orange.cohort.datamodel.Sample;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMapping;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMappingFile;
import com.hartwig.hmftools.orange.cohort.mapping.DoidCohortMapper;
import com.hartwig.hmftools.orange.cohort.percentile.CohortPercentiles;
import com.hartwig.hmftools.orange.cohort.percentile.CohortPercentilesFile;
import com.hartwig.hmftools.orange.cohort.percentile.CohortPercentilesModel;
import com.hartwig.hmftools.orange.cohort.percentile.PercentileType;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CohortPercentileEvaluator {

    private static final Logger LOGGER = LogManager.getLogger(CohortPercentileEvaluator.class);

    private static final String RESOURCE_DIR = "/data/resources/public";
    private static final String DOID_JSON = RESOURCE_DIR + "/disease_ontology/201015_doid.json";
    private static final String COHORT_MAPPING_TSV = RESOURCE_DIR + "/orange/cohort_mapping.tsv";
    private static final String COHORT_PERCENTILES_TSV = RESOURCE_DIR + "/orange/cohort_percentiles.tsv";

    private static final String OUTPUT_EVALUATION_TSV = "/data/experiments/orange_cohort/sample_sv_tmb_percentile_evaluations.tsv";

    public static void main(String[] args) throws ParseException, IOException {
        LOGGER.info("Running ORANGE Cohort Percentile Evaluator");
        CommandLine cmd = new DefaultParser().parse(createOptions(), args);

        DatabaseAccess database = DatabaseAccess.createDatabaseAccess(cmd);

        LOGGER.info("Querying database");
        List<Sample> samples = SampleQuery.run(database);
        List<Observation> observations = SvTmbQuery.run(database, samples);
        database.close();

        LOGGER.info("Creating percentiles model");
        CohortPercentilesModel model = createModel();
        Map<Observation, Evaluation> evaluations = Maps.newHashMap();
        for (Observation observation : observations) {
            evaluations.put(observation, model.percentile(observation));
        }

        List<String> lines = Lists.newArrayList();
        lines.add(header());
        lines.addAll(toLines(evaluations));

        LOGGER.info("Writing output evaluations to {}", OUTPUT_EVALUATION_TSV);
        Files.write(new File(OUTPUT_EVALUATION_TSV).toPath(), lines);
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        addDatabaseCmdLineArgs(options);
        return options;
    }

    @NotNull
    private static CohortPercentilesModel createModel() throws IOException {
        LOGGER.info(" Reading percentiles from {}", COHORT_PERCENTILES_TSV);
        Multimap<PercentileType, CohortPercentiles> percentilesMap = CohortPercentilesFile.read(COHORT_PERCENTILES_TSV);

        LOGGER.info(" Reading DOIDs from {}", DOID_JSON);
        DoidEntry doidEntry = DiseaseOntology.readDoidOwlEntryFromDoidJson(DOID_JSON);
        DoidParents doidParentModel = DoidParents.fromEdges(doidEntry.edges());

        LOGGER.info(" Reading cohort mappings from {}", COHORT_MAPPING_TSV);
        List<CohortMapping> mappings = CohortMappingFile.read(COHORT_MAPPING_TSV);

        return new CohortPercentilesModel(new DoidCohortMapper(doidParentModel, mappings), percentilesMap);
    }

    @NotNull
    private static String header() {
        return new StringJoiner("\t").add("sampleId")
                .add("type")
                .add("value")
                .add("cancerType")
                .add("panCancerPercentile")
                .add("cancerTypePercentile")
                .toString();
    }

    @NotNull
    private static List<String> toLines(@NotNull Map<Observation, Evaluation> evaluations) {
        List<String> lines = Lists.newArrayList();
        for (Map.Entry<Observation, Evaluation> evaluation : evaluations.entrySet()) {
            lines.add(toLine(evaluation.getKey(), evaluation.getValue()));
        }
        return lines;
    }

    @NotNull
    private static String toLine(@NotNull Observation observation, @NotNull Evaluation evaluation) {
        return new StringJoiner("\t").add(observation.sample().sampleId())
                .add(observation.type().toString())
                .add(String.valueOf(observation.value()))
                .add(evaluation.cancerType())
                .add(String.valueOf(evaluation.panCancerPercentile()))
                .add(String.valueOf(evaluation.cancerTypePercentile()))
                .toString();
    }
}

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
import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidEntry;
import com.hartwig.hmftools.common.doid.DoidParents;
import com.hartwig.hmftools.orange.cohort.datamodel.Sample;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMapper;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMapping;
import com.hartwig.hmftools.orange.cohort.mapping.CohortMappingFile;
import com.hartwig.hmftools.orange.cohort.mapping.DoidCohortMapper;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CohortMapperEvaluator {

    private static final Logger LOGGER = LogManager.getLogger(CohortMapperEvaluator.class);

    private static final String RESOURCE_DIR = "/data/resources/public";
    private static final String DOID_JSON = RESOURCE_DIR + "/disease_ontology/doid.json";
    private static final String COHORT_MAPPING_TSV = RESOURCE_DIR + "/orange/cohort_mapping.tsv";

    private static final String OUTPUT_EVALUATION_TSV = "/data/experiments/orange/all_mapped_samples.tsv";

    public static void main(String[] args) throws IOException, ParseException {
        LOGGER.info("Running ORANGE Cohort Mapper Evaluator");

        CommandLine cmd = new DefaultParser().parse(createOptions(), args);

        DatabaseAccess database = DatabaseAccess.createDatabaseAccess(cmd);

        LOGGER.info("Loading samples from database");
        List<Sample> samples = SampleQuery.run(database);
        LOGGER.info(" Loaded {} samples from database", samples.size());
        database.close();

        CohortMapper mapper = createCohortMapper();

        Map<Sample, String> evaluations = Maps.newHashMap();
        for (Sample sample : samples) {
            evaluations.put(sample, mapper.cancerTypeForSample(sample));
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
    private static CohortMapper createCohortMapper() throws IOException {
        LOGGER.info("Loading cohort mappings from {}", COHORT_MAPPING_TSV);
        List<CohortMapping> mappings = CohortMappingFile.read(COHORT_MAPPING_TSV);
        LOGGER.info(" Loaded {} mappings", mappings.size());

        LOGGER.info("Reading DOID model from {}", DOID_JSON);
        DoidEntry doidEntry = DiseaseOntology.readDoidOwlEntryFromDoidJson(DOID_JSON);
        DoidParents doidParentModel = DoidParents.fromEdges(doidEntry.edges());

        return new DoidCohortMapper(doidParentModel, mappings);
    }

    @NotNull
    private static String header() {
        return new StringJoiner("\t").add("sampleId").add("cancerType").toString();
    }

    @NotNull
    private static List<String> toLines(@NotNull Map<Sample, String> evaluations) {
        List<String> lines = Lists.newArrayList();
        for (Map.Entry<Sample, String> evaluation : evaluations.entrySet()) {
            lines.add(toLine(evaluation.getKey(), evaluation.getValue()));
        }
        return lines;
    }

    @NotNull
    private static String toLine(@NotNull Sample sample, @NotNull String cancerType) {
        return new StringJoiner("\t").add(sample.sampleId()).add(cancerType).toString();
    }
}

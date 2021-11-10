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
import com.hartwig.hmftools.orange.cohort.mapping.CohortMapper;
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
import org.jetbrains.annotations.NotNull;

public class CohortPercentileEvaluator {

    public static void main(String[] args) throws ParseException, IOException {
        CommandLine cmd = new DefaultParser().parse(createOptions(), args);

        DatabaseAccess database = DatabaseAccess.createDatabaseAccess(cmd);

        List<Sample> samples = SampleQuery.run(database);
        List<Observation> observations = SvTmbQuery.run(database, samples);

        String dir = "/data/experiments/orange_cohort";
        String percentileTsv = CohortPercentilesFile.generateOutputTsv(dir);
        Multimap<PercentileType, CohortPercentiles> percentilesMap = CohortPercentilesFile.read(percentileTsv);

        String doidJson = "/data/resources/public/disease_ontology/201015_doid.json";
        DoidEntry doidEntry = DiseaseOntology.readDoidOwlEntryFromDoidJson(doidJson);
        DoidParents doidParentModel = DoidParents.fromEdges(doidEntry.edges());

        String cohortMappingTsv = "/data/experiments/orange_cohort/orange_cohort_mapping.tsv";
        List<CohortMapping> mappings = CohortMappingFile.read(cohortMappingTsv);
        CohortMapper mapper = new DoidCohortMapper(doidParentModel, mappings);

        CohortPercentilesModel model = new CohortPercentilesModel(mapper, percentilesMap);
        Map<Observation, Evaluation> evaluations = Maps.newHashMap();
        for (Observation observation : observations) {
            evaluations.put(observation, model.percentile(observation));
        }

        List<String> lines = com.google.common.collect.Lists.newArrayList();
        lines.add(header());
        lines.addAll(toLines(evaluations));

        String outputTsv = dir + File.separator + "sample_sv_tmb_percentile_evaluations.tsv";
        Files.write(new File(outputTsv).toPath(), lines);
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        addDatabaseCmdLineArgs(options);
        return options;
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

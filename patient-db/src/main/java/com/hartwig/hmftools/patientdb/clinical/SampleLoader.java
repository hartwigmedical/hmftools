package com.hartwig.hmftools.patientdb.clinical;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.lims.Lims;
import com.hartwig.hmftools.common.reportingdb.ReportingDatabase;
import com.hartwig.hmftools.common.reportingdb.ReportingEntry;
import com.hartwig.hmftools.common.runcontext.RunContext;
import com.hartwig.hmftools.patientdb.clinical.datamodel.SampleData;
import com.hartwig.hmftools.patientdb.clinical.readers.LimsSampleReader;
import com.hartwig.hmftools.patientdb.clinical.readers.RunsFolderReader;
import com.hartwig.hmftools.patientdb.clinical.readers.RunsJsonReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SampleLoader {

    private static final Logger LOGGER = LogManager.getLogger(SampleLoader.class);

    @NotNull
    private final Lims lims;

    public SampleLoader(@NotNull final Lims lims) {
        this.lims = lims;
    }

    @NotNull
    public Map<String, List<SampleData>> loadSamplesPerPatient(@NotNull ClinicalAlgoConfig config) throws IOException {
        List<RunContext> runContexts = extractRunContextsFromDirectoryOrJson(config, config.runsDirectory(), config.runsJson());
        Map<String, List<String>> sequencedSamplesPerPatient = extractSequencedSamplesFromRunContexts(runContexts);
        Map<String, String> sampleToSetNameMap = extractSampleToSetNameMap(runContexts);

        LOGGER.info("  Loaded sequence runs for {} patient IDs ({} samples)",
                sequencedSamplesPerPatient.keySet().size(),
                toUniqueSampleIds(sequencedSamplesPerPatient).size());

        LOGGER.info(" Loading reporting database from {}", config.reportingDbTsv());
        List<ReportingEntry> reportingEntries = ReportingDatabase.read(config.reportingDbTsv());
        LOGGER.info("  Loaded {} reporting entries", reportingEntries.size());

        LOGGER.info(" Loading sample data from LIMS in {}", config.limsDirectory());
        Map<String, List<SampleData>> sampleDataPerPatient =
                extractAllSamplesFromLims(lims, reportingEntries, sampleToSetNameMap, sequencedSamplesPerPatient);
        LOGGER.info("  Loaded samples for {} patient IDs ({} samples)",
                sampleDataPerPatient.keySet().size(),
                countValues(sampleDataPerPatient));

        return sampleDataPerPatient;
    }

    @NotNull
    private List<RunContext> extractRunContextsFromDirectoryOrJson(@NotNull ClinicalAlgoConfig config, @Nullable String runsDirectory,
            @Nullable String runsJson) throws IOException {
        if (runsDirectory != null && runsJson == null) {
            String pipelineVersionFile = config.pipelineVersionFile();
            if (pipelineVersionFile != null) {
                LOGGER.info("   Loading sequence runs from {}", runsDirectory);
                return loadRunContextsFromDirectory(runsDirectory, pipelineVersionFile);
            } else {
                throw new IllegalArgumentException("Cannot load from a runs directory if no pipeline version file given");
            }
        } else if (runsJson != null && runsDirectory == null) {
            return loadRunContextsFromJson(runsJson);
        }
        throw new IllegalStateException("Either a runs directory or runs json must be specified, and not both");
    }

    @NotNull
    private static List<RunContext> loadRunContextsFromJson(final String runsJson) {
        List<RunContext> runContexts = RunsJsonReader.extractRunContexts(new File(runsJson));
        LOGGER.info("  Loaded run contexts from {} ({} sets)", runsJson, runContexts.size());
        return runContexts;
    }

    @NotNull
    private static List<RunContext> loadRunContextsFromDirectory(@NotNull String runsDirectory, @NotNull String pipelineVersionFile)
            throws IOException {
        List<RunContext> runContexts = RunsFolderReader.extractRunContexts(new File(runsDirectory), pipelineVersionFile);
        LOGGER.info("  Loaded run contexts from {} ({} sets)", runsDirectory, runContexts.size());

        return runContexts;
    }

    @NotNull
    private static Map<String, String> extractSampleToSetNameMap(@NotNull List<RunContext> runContexts) {
        Map<String, String> sampleToSetNameMap = Maps.newHashMap();
        for (RunContext runContext : runContexts) {
            if (sampleToSetNameMap.containsKey(runContext.tumorSample())) {
                LOGGER.warn(" Duplicate sample ID found in run contexts: {}", runContext.tumorSample());
            }
            sampleToSetNameMap.put(runContext.tumorSample(), runContext.setName());
        }
        return sampleToSetNameMap;
    }

    @NotNull
    private static Map<String, List<String>> extractSequencedSamplesFromRunContexts(@NotNull List<RunContext> runContexts) {
        Map<String, List<String>> sequencedSamplesPerPatient = Maps.newHashMap();
        for (RunContext runContext : runContexts) {
            String patientId = extractPatientIdentifier(runContext.setName());
            List<String> currentSampleIds = sequencedSamplesPerPatient.get(patientId);
            if (currentSampleIds == null) {
                currentSampleIds = Lists.newArrayList(runContext.tumorSample());
            } else {
                currentSampleIds.add(runContext.tumorSample());
            }
            sequencedSamplesPerPatient.put(patientId, currentSampleIds);
        }

        return sequencedSamplesPerPatient;
    }

    @NotNull
    private static Map<String, List<SampleData>> extractAllSamplesFromLims(@NotNull Lims lims,
            @NotNull List<ReportingEntry> reportingEntries, @NotNull Map<String, String> sampleToSetNameMap,
            @NotNull Map<String, List<String>> sequencedSamplesPerPatient) {
        LimsSampleReader sampleReader =
                new LimsSampleReader(lims, reportingEntries, sampleToSetNameMap, toUniqueSampleIds(sequencedSamplesPerPatient));

        Map<String, List<SampleData>> samplesPerPatient = Maps.newHashMap();
        for (String sampleBarcode : lims.sampleBarcodes()) {
            String sampleId = lims.sampleId(sampleBarcode);
            SampleData sampleData = sampleReader.read(sampleBarcode, sampleId);

            if (sampleData != null) {
                String patientId = lims.patientId(sampleBarcode);
                List<SampleData> currentSamples = samplesPerPatient.get(patientId);
                if (currentSamples == null) {
                    currentSamples = Lists.newArrayList(sampleData);
                } else if (!sampleIdExistsInSampleDataList(currentSamples, sampleId)) {
                    // Ideally if a single sample exists in LIMS with multiple barcodes we pick "the most relevant one".
                    // Currently just picking a random one - sampleId has to be unique in this list.
                    currentSamples.add(sampleData);
                }
                samplesPerPatient.put(patientId, currentSamples);
            }
        }

        // Some samples may be missing from LIMS simply because they are old and we did not collect information back in those days.
        for (Map.Entry<String, List<String>> sequencedPatientEntry : sequencedSamplesPerPatient.entrySet()) {
            List<SampleData> samples = samplesPerPatient.get(sequencedPatientEntry.getKey());
            if (samples == null) {
                samples = Lists.newArrayList();
            }
            for (String sampleId : sequencedPatientEntry.getValue()) {
                if (!sampleIdExistsInSampleDataList(samples, sampleId)) {
                    LOGGER.info("  Creating sample data for {}. This sample is not found in LIMS even though it has been sequenced!",
                            sampleId);
                    SampleData sampleData = sampleReader.readSequencedSampleWithoutBarcode(sampleId);
                    if (sampleData != null) {
                        samples.add(sampleData);
                    }
                }
            }
            samplesPerPatient.put(sequencedPatientEntry.getKey(), samples);
        }

        return samplesPerPatient;
    }

    private static boolean sampleIdExistsInSampleDataList(@NotNull List<SampleData> samples, @NotNull String sampleId) {
        for (SampleData sample : samples) {
            if (sample.sampleId().equals(sampleId)) {
                return true;
            }
        }

        return false;
    }

    private static <V, K> int countValues(@NotNull Map<V, List<K>> map) {
        int count = 0;
        for (Map.Entry<V, List<K>> entry : map.entrySet()) {
            count += entry.getValue().size();
        }
        return count;
    }

    @NotNull
    private static Set<String> toUniqueSampleIds(@NotNull Map<String, List<String>> samplesPerPatient) {
        Set<String> uniqueSampleIds = Sets.newHashSet();
        for (Map.Entry<String, List<String>> entry : samplesPerPatient.entrySet()) {
            uniqueSampleIds.addAll(entry.getValue());
        }
        return uniqueSampleIds;
    }

    @NotNull
    private static String extractPatientIdentifier(@NotNull String setName) {
        String[] names = setName.split("_");
        if (names.length < 5) {
            LOGGER.error(" Run name {} had less than 5 parts after splitting on _", setName);
            return Strings.EMPTY;
        }
        return names[4];
    }
}

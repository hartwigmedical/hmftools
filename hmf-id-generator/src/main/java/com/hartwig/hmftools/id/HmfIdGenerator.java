package com.hartwig.hmftools.id;

import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.id.HmfIdConfig.DATA_DELIM;
import static com.hartwig.hmftools.id.HmfIdConfig.ID_LOGGER;
import static com.hartwig.hmftools.id.HmfSample.DELETED;
import static com.hartwig.hmftools.id.HmfSample.PATIENT_ID;
import static com.hartwig.hmftools.id.HmfSample.SAMPLE_HASH;
import static com.hartwig.hmftools.id.HmfSample.SAMPLE_ID;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.amber.AmberAnonymous;
import com.hartwig.hmftools.common.amber.AmberPatient;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class HmfIdGenerator
{
    private final HmfIdConfig mConfig;

    private final DatabaseAccess mDbAccess;


    public HmfIdGenerator(final CommandLine cmd)
    {
        mConfig = new HmfIdConfig(cmd);
        mDbAccess = DatabaseAccess.createDatabaseAccess(cmd);
    }

    public void run() throws IOException
    {
        // ID_LOGGER.info("resource reference directory: {}", mResourceRepoDir);

        if(mConfig.Restore)
        {
            restoreFromFile();
            return;
        }

        ID_LOGGER.info("running HMF Id Generator");

        // Retrieve Data
        List<AmberPatient> amberPatients = mDbAccess.readAmberPatients();
        ID_LOGGER.info("Retrieved ${amberPatients.size} samples from amberPatient table");

        List<AmberAnonymous> currentDbAnonymous = mDbAccess.readAmberAnonymous();
        List<HmfSample> currentDbSamples = currentDbAnonymous.stream().map(x -> HmfSample.fromAmberAnonymous(x)).collect(Collectors.toList());
        ID_LOGGER.info("retrieved {}} sample mappings from amberAnonymous table", currentDbSamples.size());

        List<HmfSample> currentFileSamples = loadSampleFile();

        // validate amberPatient and amberAnonymous tables
        if(!validAmberPatients(currentDbSamples, amberPatients))
        {
            ID_LOGGER.error("amberPatient and amberAnonymous are not synchronized");
            System.exit(1);
        }

        /*
        // validate synchronization between database and file
        val currentDatabaseCsv = currentAnonymous.map { x -> x.toCsv(oldGenerator) }.toSet()
        if (!databaseAndFileInSync(currentDatabaseCsv, currentFileAnonymous))
        {
            ID_LOGGER.error("database and file are not synchronized")
            System.exit(1);
        }
         */

        ID_LOGGER.info("processing samples");

        List<HmfSample> newSamples = anonymize(amberPatients, currentDbSamples);

        // Write to file and database
        ID_LOGGER.info("writing {} samples hashes to {}", newSamples.size(), mConfig.OutputHashFile);
        // CsvWriter.writeCSV(result.map { x -> x.toCsv(newGenerator) }, hashFileOut)

        List<AmberAnonymous> newAnonymous = newSamples.stream().map(x -> x.toAmberAnonymous()).collect(Collectors.toList());
        ID_LOGGER.info("writing {} sample mappings to database", newAnonymous.size());
        mDbAccess.writeAmberAnonymous(newAnonymous);

        ID_LOGGER.info("ID generation complete");
    }

    public static List<HmfSample> anonymize(final List<AmberPatient> amberPatients, final List<HmfSample> existingSamples)
    {
        List<HmfSample> newSamples = Lists.newArrayList();

        int maxPatientId = existingSamples.stream().mapToInt(x -> x.PatientId).max().orElse(0);
        Set<Integer> amberPatientIds = Sets.newHashSet();
        amberPatients.forEach(x -> amberPatientIds.add(x.patientId()));

        for(Integer patientId : amberPatientIds)
        {
            List<String> matchingSampleHashes = amberPatients.stream()
                    .filter(x -> x.patientId() == patientId).map(x -> x.sample()).collect(Collectors.toList());

            List<HmfSample> patientSamples = existingSamples.stream()
                    .filter(x -> matchingSampleHashes.contains(x.SampleHash)).collect(Collectors.toList());

            int hmfPatientId = patientSamples.stream().mapToInt(x -> x.PatientId).max().orElse(-1);
            if(hmfPatientId < 0)
                hmfPatientId = ++maxPatientId;

            int maxHmfSampleId = patientSamples.stream().mapToInt(x -> x.SampleId).max().orElse(0);

            for(String sampleHash : matchingSampleHashes)
            {
                if(existingSamples.stream().noneMatch(x -> x.SampleHash.equals(sampleHash)))
                {
                    newSamples.add(new HmfSample(hmfPatientId, ++maxHmfSampleId, sampleHash, false));
                }
            }
        }

        List<HmfSample> allSamples = Lists.newArrayList(existingSamples);
        allSamples.addAll(newSamples);
        Collections.sort(allSamples, new HmfSample.SampleComparator());
        return allSamples;
    }

    private boolean validAmberPatients(final List<HmfSample> dbSamples, final List<AmberPatient> patients)
    {
        /*
        val actualSamples = patients.map { x -> x.sample() }.toSet()
        val expectedSamples = dbSamples.filter { x -> !x.deleted }.map { x -> x.sample }
        val missingSamples = expectedSamples subtract actualSamples
        for (missingSample in missingSamples) {
            logger.error("Missing sample $missingSample from amberPatient table")
        }

        val deletedSamples = dbSamples.filter { x -> x.deleted }.map { x -> x.sample }
        val unexpectedSamples = actualSamples intersect deletedSamples
        for (unexpectedSample in unexpectedSamples) {
            logger.error("Deleted sample $unexpectedSample in amberPatient table")
        }

        if (missingSamples.isNotEmpty() || unexpectedSamples.isNotEmpty()) {
            return false
        }
        */

        return true;
    }

    private boolean databaseAndFileInSync() // currentDatabaseCsv: Set<HmfSampleCsv>, currentFileCsv: Set<HmfSampleCsv>
    {
        /*
        // Ignore deleted flag for moment
        val adjustedCurrentDatabaseCsv = currentDatabaseCsv.map { x -> x.copy(deleted = false.toString()) }
        val adjustedCurrentFileCsv = currentFileCsv.map { x -> x.copy(deleted = false.toString()) }

        val missingFromDatabase = adjustedCurrentFileCsv subtract adjustedCurrentDatabaseCsv
        val missingFromFile = adjustedCurrentDatabaseCsv subtract adjustedCurrentFileCsv
        if (missingFromDatabase.isNotEmpty()) {
            for (missing in missingFromDatabase) {
                logger.error("${missing.hmfSampleId} missing from database")
            }
        }

        if (missingFromFile.isNotEmpty()) {
            for (missing in missingFromFile) {
                logger.error("${missing.hmfSampleId} missing from file")
            }
        }

        if (missingFromDatabase.isNotEmpty() || missingFromFile.isNotEmpty()) {
            return false
        }
        */

        return true;
    }

    private void restoreFromFile()
    {
        ID_LOGGER.info("precomputing hashes");

        Map<String,String> hashMap = HashGenerator.precomputeHashes();

        List<HmfSample> existingSamples = loadSampleFile();
        List<HmfSample> newSamples = loadSampleFile();

        for(HmfSample sample : existingSamples)
        {
            String sampleFromHash = hashMap.get(sample.SampleHash);

            if(sampleFromHash != null)
            {
                HmfSample newSample = new HmfSample(sample.PatientId, sample.SampleId, sampleFromHash, sample.Deleted);
                newSamples.add(newSample);
            }
            else
            {
                ID_LOGGER.warn("unable to find hash for sample({})", sample);
            }
        }

        if(!newSamples.isEmpty())
        {
            ID_LOGGER.info("writing {} samples to database");
            List<AmberAnonymous> amberAnonymous = newSamples.stream().map(x -> x.toAmberAnonymous()).collect(Collectors.toList());
            mDbAccess.writeAmberAnonymous(amberAnonymous);
        }

        /*
            ID_LOGGER.info("Precomputing hashes")
    val hashes = precomputeHashes(generator)

    val currentFileAnonymous = CsvReader.readCSVByName<HmfSampleCsv>(hashFileIn).toSet()
    ID_LOGGER.info("Retrieved ${currentFileAnonymous.size} sample mappings from $hashFileIn")

    val result = mutableListOf<HmfSample>()
    for (fileEntry in currentFileAnonymous) {
        val sampleFromHash = hashes[fileEntry.sampleHash]
        if (sampleFromHash != null) {
            val record = HmfSample(fileEntry.patientId.toInt(), fileEntry.sampleId.toInt(), fileEntry.deleted.toBoolean(), sampleFromHash)
            result.add(record)
        } else {
            ID_LOGGER.warn("Unable to determine sample of ${fileEntry.hmfSampleId}")
        }
    }

    if (result.isNotEmpty()) {
        ID_LOGGER.info("Writing ${result.size} sample mappings to database")
        mDbAccess.writeAmberAnonymous(result.map { x -> x.toAmberAnonymous() })
    }

         */
    }

    public List<HmfSample> loadSampleFile()
    {
        List<HmfSample> sampleDataList = Lists.newArrayList();

        try
        {
            final List<String> fileData = Files.readAllLines(new File(mConfig.InputHashFile).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DATA_DELIM);

            int patientIdIndex = fieldsIndexMap.get(PATIENT_ID);
            int sampleIdIndex = fieldsIndexMap.get(SAMPLE_ID);
            // int hmfSampleIdIndex = fieldsIndexMap.get(HMF_SAMPLE_ID);
            int sampleHashIndex = fieldsIndexMap.get(SAMPLE_HASH);
            int deletedIndex = fieldsIndexMap.get(DELETED);

            for(final String line : fileData)
            {
                String[] values = line.split(DATA_DELIM);

                HmfSample sample = new HmfSample(
                        Integer.parseInt(values[patientIdIndex]), Integer.parseInt(values[sampleIdIndex]),
                        values[sampleHashIndex], Boolean.parseBoolean(values[deletedIndex]));

                sampleDataList.add(sample);
            }

            ID_LOGGER.info("loaded {} samples from file({})", sampleDataList.size(), mConfig.InputHashFile);
        }
        catch (IOException e)
        {
            ID_LOGGER.error("failed to read sample data file({}): {}", mConfig.InputHashFile, e.toString());
        }

        return sampleDataList;
    }

    public static void main(String[] args) throws IOException, ParseException
    {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        setLogLevel(cmd);

        HmfIdGenerator generator = new HmfIdGenerator(cmd);
        generator.run();
    }

    private static Options createOptions()
    {
        Options options = new Options();

        addOutputDir(options);
        addLoggingOptions(options);

        return options;
    }
}

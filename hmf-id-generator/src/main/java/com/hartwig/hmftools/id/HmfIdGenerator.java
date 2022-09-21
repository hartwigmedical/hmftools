package com.hartwig.hmftools.id;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.id.HmfIdConfig.DATA_DELIM;
import static com.hartwig.hmftools.id.HmfIdConfig.ID_LOGGER;
import static com.hartwig.hmftools.id.HmfIdConfig.addCmdLineArgs;
import static com.hartwig.hmftools.id.SampleData.DELETED;
import static com.hartwig.hmftools.id.SampleData.PATIENT_ID;
import static com.hartwig.hmftools.id.SampleData.SAMPLE_HASH;
import static com.hartwig.hmftools.id.SampleData.SAMPLE_INDEX;

import java.io.BufferedWriter;
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
        if(mConfig.Restore)
        {
            restoreFromFile();
            return;
        }

        ID_LOGGER.info("running HMF Id Generator");

        HashGenerator oldGenerator = new HashGenerator(mConfig.Password, mConfig.MaxPrecomputeCount);

        HashGenerator newGenerator = mConfig.Password.equals(mConfig.NewPassword) ?
                oldGenerator : new HashGenerator(mConfig.NewPassword, mConfig.MaxPrecomputeCount);

        // Retrieve Data
        List<AmberPatient> amberPatients = mDbAccess.readAmberPatients();
        ID_LOGGER.info("Retrieved ${amberPatients.size} samples from amberPatient table");

        List<AmberAnonymous> currentDbAnonymous = mDbAccess.readAmberAnonymous();
        List<SampleData> currentDbSamples = currentDbAnonymous.stream().map(x -> SampleData.fromAmberAnonymous(x)).collect(Collectors.toList());
        ID_LOGGER.info("retrieved {}} sample mappings from amberAnonymous table", currentDbSamples.size());

        // compare entries in the amberPatient and amberAnonymous tables
        if(!validAmberPatients(currentDbSamples, amberPatients))
        {
            ID_LOGGER.error("amberPatient and amberAnonymous are not synchronized");
            System.exit(1);
        }

        if(!compareDatabaseVsFile(currentDbSamples))
        {
            ID_LOGGER.error("sample file and amberAnonymous are not synchronized");
            System.exit(1);
        }

        ID_LOGGER.info("processing samples");

        // create anonymised (HMF) sampleIDs for any new samples
        List<SampleData> newDbSamples = anonymize(amberPatients, currentDbSamples);

        // Write to file and database
        writeSampleFile(newDbSamples, newGenerator);

        List<AmberAnonymous> newAnonymous = newDbSamples.stream().map(x -> x.toAmberAnonymous()).collect(Collectors.toList());
        ID_LOGGER.info("writing {} sample mappings to database", newAnonymous.size());
        mDbAccess.writeAmberAnonymous(newAnonymous);

        ID_LOGGER.info("ID generation complete");
    }

    public static List<SampleData> anonymize(final List<AmberPatient> amberPatients, final List<SampleData> existingDbSamples)
    {
        // searches in amberPatients for samples without an entry in amberAnonymous, and for any new ones an HmfSampleId
        // is generated using the existing patientId and sampleIndex values
        List<SampleData> newSamples = Lists.newArrayList();

        // find the maximum patient ID (including previously deleted samples to avoid reuse of an ID)
        int maxPatientId = existingDbSamples.stream().mapToInt(x -> x.PatientId).max().orElse(0);

        // extract unique patient IDs
        Set<Integer> amberPatientIds = Sets.newHashSet();
        amberPatients.forEach(x -> amberPatientIds.add(x.patientId()));

        for(Integer patientId : amberPatientIds)
        {
            List<String> matchingSampleIds = amberPatients.stream()
                    .filter(x -> x.patientId() == patientId).map(x -> x.sample()).collect(Collectors.toList());

            List<SampleData> patientSamples = existingDbSamples.stream()
                    .filter(x -> matchingSampleIds.contains(x.SampleId)).collect(Collectors.toList());

            int hmfPatientId = patientSamples.stream().mapToInt(x -> x.PatientId).max().orElse(-1);
            if(hmfPatientId < 0)
                hmfPatientId = ++maxPatientId;

            int maxHmfSampleId = patientSamples.stream().mapToInt(x -> x.SampleIndex).max().orElse(0);

            for(String sampleId : matchingSampleIds)
            {
                if(existingDbSamples.stream().noneMatch(x -> x.SampleId.equals(sampleId)))
                {
                    newSamples.add(new SampleData(hmfPatientId, ++maxHmfSampleId, sampleId, "", false));
                }
            }
        }

        List<SampleData> allSamples = Lists.newArrayList(existingDbSamples);
        allSamples.addAll(newSamples);
        Collections.sort(allSamples, new SampleData.SampleComparator());
        return allSamples;
    }

    private boolean validAmberPatients(final List<SampleData> dbSamples, final List<AmberPatient> patients)
    {
        // extract original / non-anonymised sampleIDs from each source
        Set<String> actualSamples = patients.stream().map(x -> x.sample()).collect(Collectors.toSet());
        Set<String> expectedSamples = dbSamples.stream().filter(x -> !x.Deleted).map(x -> x.SampleId).collect(Collectors.toSet());
        Set<String> missingSamples = expectedSamples.stream().filter(x -> !actualSamples.contains(x)).collect(Collectors.toSet());

        for(String sample : missingSamples)
        {
            ID_LOGGER.error("missing sample({}) from amberPatient table", sample);
        }

        Set<String> deleteSamples = dbSamples.stream().filter(x -> x.Deleted).map(x -> x.SampleId).collect(Collectors.toSet());
        Set<String> unexpectedSamples = actualSamples.stream().filter(x -> deleteSamples.contains(x)).collect(Collectors.toSet());

        for(String sample : unexpectedSamples)
        {
            ID_LOGGER.error("deleted sample({}) in amberPatient table", sample);
        }

        if(!missingSamples.isEmpty() || !unexpectedSamples.isEmpty())
            return false;

        return true;
    }

    private boolean compareDatabaseVsFile(final List<SampleData> currentDbSamples)
    {
        List<SampleData> currentFileSamples = loadSampleFile();

        List<SampleData> missingFromFile = Lists.newArrayList();

        for(SampleData dbSample : currentDbSamples)
        {
            boolean matched = false;

            for(int i = 0; i < currentFileSamples.size(); ++i)
            {
                SampleData fileSample = currentFileSamples.get(i);

                if(dbSample.matches(fileSample))
                {
                    currentFileSamples.remove(i);
                    matched = true;
                    break;
                }
            }

            if(!matched)
                missingFromFile.add(dbSample);
        }

        for(SampleData sample : missingFromFile)
        {
            ID_LOGGER.error("sample({}) missing from file");
        }

        for(SampleData sample : currentFileSamples)
        {
            ID_LOGGER.error("sample({}) missing from database");
        }

        return currentFileSamples.isEmpty() && missingFromFile.isEmpty();
    }

    private void restoreFromFile()
    {
        // this method seems flawed since the pre-compute step cannot handle the current variation in actual sampleIDs
        ID_LOGGER.info("precomputing hashes");

        HashGenerator generator = new HashGenerator(mConfig.Password, mConfig.MaxPrecomputeCount);
        Map<String,String> hashMap = generator.precomputeHashes();

        List<SampleData> existingSamples = loadSampleFile();
        List<SampleData> newSamples = Lists.newArrayList();

        for(SampleData sample : existingSamples)
        {
            String sampleId = hashMap.get(sample.sampleHash());

            if(sampleId != null)
            {
                SampleData newSample = new SampleData(sample.PatientId, sample.SampleIndex, sampleId, sample.sampleHash(), sample.Deleted);
                newSamples.add(newSample);
            }
            else
            {
                ID_LOGGER.warn("unable to find hash for sample({})", sample);
            }
        }

        if(!newSamples.isEmpty())
        {
            ID_LOGGER.info("writing {} samples to database", newSamples.size());
            List<AmberAnonymous> amberAnonymous = newSamples.stream().map(x -> x.toAmberAnonymous()).collect(Collectors.toList());
            mDbAccess.writeAmberAnonymous(amberAnonymous);
        }
    }

    public List<SampleData> loadSampleFile()
    {
        List<SampleData> sampleDataList = Lists.newArrayList();

        try
        {
            final List<String> fileData = Files.readAllLines(new File(mConfig.InputHashFile).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DATA_DELIM);

            int patientIdIndex = fieldsIndexMap.get(PATIENT_ID);
            int sampleIdIndex = fieldsIndexMap.get(SAMPLE_INDEX);
            // int hmfSampleIdIndex = fieldsIndexMap.get(HMF_SAMPLE_ID);
            int sampleHashIndex = fieldsIndexMap.get(SAMPLE_HASH);
            int deletedIndex = fieldsIndexMap.get(DELETED);

            for(final String line : fileData)
            {
                String[] values = line.split(DATA_DELIM);

                SampleData sample = new SampleData(
                        Integer.parseInt(values[patientIdIndex]),
                        Integer.parseInt(values[sampleIdIndex]),
                        "", // original sampleID is not stored in this file
                        values[sampleHashIndex],
                        Boolean.parseBoolean(values[deletedIndex]));

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

    private void writeSampleFile(final List<SampleData> samples, final HashGenerator hashGenerator)
    {
        ID_LOGGER.info("regenerating samples hashes for {} samples", samples.size());
        samples.forEach(x -> x.setSampleHash(hashGenerator.hash(x.SampleId)));

        ID_LOGGER.info("writing {} samples hashes to {}", samples.size(), mConfig.OutputHashFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.OutputHashFile, false);

            writer.write(SampleData.header());
            writer.newLine();

            for(SampleData sample : samples)
            {
                writer.write(sample.toCsv());
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            ID_LOGGER.error(" failed to write sample data: {}", e.toString());
        }
    }

    public static void main(String[] args) throws IOException, ParseException
    {
        Options options = new Options();
        addCmdLineArgs(options);

        CommandLine cmd = new DefaultParser().parse(options, args);

        setLogLevel(cmd);

        HmfIdGenerator generator = new HmfIdGenerator(cmd);
        generator.run();
    }
}
